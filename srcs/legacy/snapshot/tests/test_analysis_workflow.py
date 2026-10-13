import copy
import logging
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from PIL import Image
import tifffile
from srcs.utils.config import RunContext
from srcs.prep.resume import StepRunner, check_tiff, check_registration, check_alignment, check_filtered
from srcs.comp.ID_Tools import Edgewise_Contribution

ROOT = Path(__file__).resolve().parents[1]

def config(folder):
    return {
        "paths": {"raw": str(folder / "raw"), "processed": str(folder / "processed"),
                  "output": str(folder / "output"), "mask": str(folder / "mask")},
        "variables": {"dataset": "", "Week_ids": ["S1", "S2", "S3"],
                      "Mouse_ids": [f"M{i}" for i in range(6)], "Run_ids": ["T1", "T2", "T3"],
                      "Atlas": [["A"], ["B"], ["C"]]},
        "ID": {"group_mode": "random", "session_seed": 42, "shuffle": True, "n_iter": 3},
        "Edge": {"seed": 42, "repeat": 2, "threshold": 99, "gap": .25,
                 "cross_only": False, "save_mat": True},
        "analysis": {"statistics_file": "custom_stats.pkl"},
        "supplement": {"structure_atlases": [["A"], ["B"], ["C"]],
                       "structure_group_names": ["A", "B", "C"], "exclude": []},
    }

class PlotStub:
    def __init__(self, *args, **kwargs):
        pass
    def __getattr__(self, name):
        return lambda *args, **kwargs: self

class AnalysisTest(unittest.TestCase):
    def test_output_identity_and_filename(self):
        cfg = config(Path("relative"))
        ctx = RunContext(None, cfg, logging.getLogger("test"))
        first = ctx.output_path
        self.assertTrue(first.startswith(str(Path("relative") / "output")))
        self.assertIn("seed_random_42", first)
        self.assertEqual(Path(ctx.statistics_path).name, "custom_stats.pkl")
        self.assertEqual(Path(ctx.statistics_path).parent.name, "Identification")
        self.assertEqual(first, RunContext(None, copy.deepcopy(cfg), ctx.logger).output_path)
        cfg["ID"]["session_seed"] = 43
        self.assertNotEqual(first, ctx.output_path)
        cfg["ID"]["session_seed"] = 42
        cfg["ID"]["group_mode"] = "sequential"
        self.assertNotEqual(first, ctx.output_path)
        cfg["ID"]["repeat_id"] = "S1"
        self.assertIn("repeatS1", ctx.output_path)
        cfg["analysis"]["statistics_file"] = "../escape.pkl"
        with self.assertRaises(ValueError):
            _ = ctx.statistics_path

    def test_all_axes_modes_and_cover_statistics_pipeline(self):
        from srcs.utils.dataset import analysis_plan
        from srcs.utils.utils import build_half_dataframe
        from srcs.comp.ID_Tools import Identify
        for axis in ('run','week'):
            for mode in ('single','cumulative'):
                with self.subTest(axis=axis,mode=mode), tempfile.TemporaryDirectory() as tmp:
                    cfg=config(Path(tmp))
                    cfg['analysis'].update(split_axis=axis,selection_mode=mode)
                    cfg['ID']['completeness_mode']='cover'
                    cfg['variables']['Week_ids']=['S0','S1','S2','S3']
                    ctx=RunContext(None,cfg,logging.getLogger('test'))
                    Path(ctx.correlation_path).mkdir(parents=True)
                    labels=['SM-r1_L','SM-r2_L','VIS-r3_R','VIS-r4_R']
                    rng=np.random.default_rng(345)
                    original={}
                    for atlas in ctx.atlases:
                        base=rng.uniform(-.3,.5,(6,4,4))
                        for w in ctx.weeks:
                            frame=pd.DataFrame(index=ctx.mice,columns=['T3','T1','T2'],dtype=object)
                            for mi,m in enumerate(ctx.mice):
                                for r in ctx.runs:
                                    mat=base[mi]+rng.normal(0,.2,(4,4))
                                    mat=np.clip((mat+mat.T)/2,-.9,.9)
                                    np.fill_diagonal(mat,1)
                                    frame.at[m,r]=None if w=='S0' else pd.DataFrame(mat,index=labels,columns=labels)
                            frame.to_pickle(Path(ctx.correlation_path)/f'{w}_correlation_matrices_{atlas[0]}_excluded.pkl')
                            original[(atlas[0],w)]=frame
                    from srcs.utils.results import save_json
                    save_json(Path(ctx.output_dir('BatchReports'))/'02_correlation_test.json',{'records':[
                        dict(week=w,mouse=m,run=r,atlas='_'.join(a),status='missing' if w=='S0' else 'success')
                        for w in ctx.weeks for m in ctx.mice for r in ctx.runs for a in ctx.atlases]})
                    with patch('srcs.utils.config.load_context',return_value=ctx), patch('srcs.vis.PlotSeries.PlotSeries',PlotStub), patch('srcs.comp.ID_Tools.PlotSeries',PlotStub), patch('srcs.vis.circularGraph.circularGraph',lambda *a,**k:None):
                        runpy.run_path(str(ROOT/'scripts/03_identification.py'),run_name='__main__')
                        runpy.run_path(str(ROOT/'scripts/04_statistics.py'),run_name='__main__')
                        runpy.run_path(str(ROOT/'scripts/05_supplement.py'),run_name='__main__')
                    weeks=['S1','S2','S3']
                    plan=analysis_plan(ctx,weeks)
                    stats=pd.read_pickle(ctx.statistics_path)
                    self.assertEqual(stats.attrs['cohort_selection']['retained_weeks'],weeks)
                    self.assertEqual(list(stats.columns),list(plan))
                    groups=ctx.split_groups(weeks,ctx.runs)
                    for tag,selection in plan.items():
                        directory=Path(ctx.output_path)/'Identification'/tag
                        numeric=pd.read_pickle(directory/'Numeric/identification_A.pkl')
                        self.assertEqual(numeric['weeks'],selection['weeks'])
                        self.assertEqual(numeric['runs'],selection['runs'])
                        self.assertEqual(numeric['mice'],ctx.mice)
                        self.assertEqual(numeric['split_groups'],groups)
                        # Construct expected Fisher-z halves directly from physical recordings.
                        halves=pd.DataFrame(index=['1st','2nd'],columns=ctx.mice,dtype=object)
                        for half,group in zip(halves.index,groups):
                            for mouse in ctx.mice:
                                wids=group if axis=='week' else selection['weeks']
                                rids=group if axis=='run' else selection['runs']
                                vals=[np.arctanh(np.clip(original[('A',w)].at[mouse,r].to_numpy(),-.999999,.999999)) for w in wids for r in rids]
                                halves.at[half,mouse]=pd.DataFrame(np.mean(vals,axis=0),index=labels,columns=labels)
                        expected=Identify(halves,atlas_exclude=['A']).ID_Test()[1]
                        np.testing.assert_allclose(numeric['similarity'],expected)
                        edge=pd.read_pickle(directory/f'EdgeContri/edges_A_{tag}.pkl')
                        self.assertEqual(edge['status'],'complete')
                        self.assertEqual(edge['split_groups'],groups)
                        shuffle=pd.read_pickle(directory/'Numeric/shuffle_A.pkl')
                        self.assertEqual(len(shuffle['shuffle_acc'][0]),3)
                    tables=Path(ctx.output_path)/'Statistics/Tables'
                    self.assertTrue((tables/'anova_A.csv').is_file())
                    self.assertTrue((tables/f'posthoc_{list(plan)[-1]}.csv').is_file())
                    observations=pd.read_csv(tables/'observations_A.csv')
                    self.assertIn('week' if axis=='run' else 'run',observations)
                    self.assertEqual(set(observations['mouse']),set(ctx.mice))
                    self.assertTrue((tables/'Supplement/manova.pkl').is_file())


class ResumeTest(unittest.TestCase):
    def test_skip_failure_restart_and_downstream_invalidation(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "source"; source.write_text("data")
            out = folder / "out"; later = folder / "later"
            logger = logging.getLogger("test")
            action = lambda: out.write_text("complete")
            valid = lambda: out.exists() and out.read_text() == "complete"
            runner = StepRunner(folder, logger)
            self.assertTrue(runner.run("one", action, [out], valid, [source]))
            runner.run("two", lambda: later.write_text("complete"), [later], later.exists, [out])
            runner = StepRunner(folder, logger)
            self.assertFalse(runner.run("one", action, [out], valid, [source]))
            self.assertFalse(runner.run("two", lambda: None, [later], later.exists, [out]))
            def fail():
                out.write_text("complete")
                raise RuntimeError("interrupted")
            source.write_text("changed input")
            with self.assertRaises(RuntimeError):
                StepRunner(folder, logger).run("one", fail, [out], valid, [source])
            runner = StepRunner(folder, logger)
            self.assertTrue(runner.run("one", action, [out], valid, [source]))
            self.assertTrue(runner.run("two", lambda: later.write_text("complete"), [later], later.exists, [out]))

    def test_file_integrity_and_frame_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "source.tif"
            data = np.ones((12, 3, 4), dtype=np.uint16)
            tifffile.imwrite(source, data, photometric="minisblack")
            self.assertTrue(check_tiff(source, data.shape))
            reg = folder / "reg"; reg.mkdir()
            tifffile.imwrite(reg / "0.tif", data[:6], photometric="minisblack")
            self.assertFalse(check_registration(reg, source))
            tifffile.imwrite(reg / "1.tif", data[6:], photometric="minisblack")
            self.assertTrue(check_registration(reg, source))
            affine = folder / "affine.npy"; np.save(affine, data)
            preview = folder / "preview.png"; Image.fromarray(data[0]).save(preview)
            self.assertTrue(check_alignment(affine, preview, data.shape))
            filtered = folder / "filtered.npy"; np.save(filtered, data[::6])
            self.assertTrue(check_filtered(filtered, affine, (3, 4)))
            np.save(filtered, data[:1])
            self.assertFalse(check_filtered(filtered, affine, (3, 4)))
            with source.open("r+b") as stream:
                stream.truncate(source.stat().st_size // 2)
            with self.assertRaises((ValueError, OSError, tifffile.TiffFileError)):
                check_tiff(source)

class PreprocessingIntegrationTest(unittest.TestCase):
    def test_five_steps_resume_and_changed_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cfg = config(folder)
            cfg["preprocess"] = {"resume": True, "keepna": False}
            ctx = RunContext(None, cfg, logging.getLogger("test"))
            raw = folder / "raw/T1/S1/M0"
            raw.mkdir(parents=True)
            (raw / "info.txt").write_text("test metadata")
            (raw / "img_00000.bin").write_bytes(b"source bytes")
            masks = folder / "mask"; masks.mkdir()
            Image.fromarray(np.full((3, 4), 255, dtype=np.uint8)).save(masks / "Chamber_mask.png")
            data = np.ones((12, 3, 4), dtype=np.uint16)
            called = []
            def bin_to_tiff(**kw):
                called.append(1)
                tifffile.imwrite(kw["output_tiff"], data, photometric="minisblack")
            def suite(**kw):
                called.append(2)
                reg = Path(kw["save_path"]) / "suite2p/plane0/reg_tif"
                reg.mkdir(parents=True)
                tifffile.imwrite(reg / "part.tif", data, photometric="minisblack")
            def merge(input_dir, output):
                called.append(3)
                tifffile.imwrite(output, data, photometric="minisblack")
            def affine(**kw):
                called.append(4)
                out = Path(kw["processed_path"])
                np.save(out / "M0_T1_Affine.npy", data)
                Image.fromarray(data[0]).save(out / "M0_T1_Average_After.png")
            def drr(**kw):
                called.append(5)
                np.save(Path(kw["processed_path"]) / "M0_T1_S1_FiltDRR.npy", data[::6])
            module = runpy.run_path(str(ROOT / "scripts/01_preprocess.py"))
            process = module["process_recording"]
            mocks = {"Bin2Tiff": bin_to_tiff, "run_isolated_process": suite,
                     "merge_tiffs": merge, "Affine_Atlas": affine, "drr_narrow": drr}
            with patch.dict(process.__globals__, mocks):
                process(ctx, "S1", "M0", "T1")
                self.assertEqual(called, [1, 2, 3, 4, 5])
                called.clear()
                process(ctx, "S1", "M0", "T1")
                self.assertEqual(called, [])
                out = folder / "processed/S1/M0/T1"
                (out / "M0_T1_S1_FiltDRR.npy").unlink()
                process(ctx, "S1", "M0", "T1")
                self.assertEqual(called, [5])
                called.clear()
                (out / "M0_T1_Average_After.png").unlink()
                process(ctx, "S1", "M0", "T1")
                self.assertEqual(called, [4, 5])
                called.clear()
                (raw / "img_00000.bin").write_bytes(b"source bytes changed")
                process(ctx, "S1", "M0", "T1")
                self.assertEqual(called, [1, 2, 3, 4, 5])
                self.assertEqual(len(list(out.glob("suite2p.previous-*"))), 1)
                called.clear()
                with patch.dict(process.__globals__, {
                    "suite2p_runtime_signature": lambda: {
                        "suite2p_version": "changed-version", "adapter_sha256": "changed-adapter"}}):
                    process(ctx, "S1", "M0", "T1")
                    self.assertEqual(called, [2, 3, 4, 5])
                    called.clear()
                    process(ctx, "S1", "M0", "T1")
                    self.assertEqual(called, [])
                self.assertEqual(len(list(out.glob("suite2p.previous-*"))), 2)



if __name__ == "__main__":
    unittest.main()
