import copy
import json
import logging
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from PIL import Image

from srcs.comp.correlation import CorrelationCache
from srcs.comp.ID_Tools import Edgewise_Contribution
from srcs.utils.config import RunContext
from srcs.utils.dataset import (axis_frames, analysis_plan, resolve_recording,
                                statistics_metadata)
from srcs.utils.utils import build_half_dataframe, fisher_z

ROOT = Path(__file__).resolve().parents[1]


def config(root):
    return dict(paths={k: str(root/k) for k in ('raw','processed','output','mask')},
                variables=dict(dataset='', Week_ids=['W0','W1','W2'],
                               Mouse_ids=['M1','M2'], Run_ids=['R1','R2'], Atlas=[['AUDv']]),
                analysis=dict(split_axis='run', selection_mode='single'),
                ID=dict(group_mode='random', session_seed=42), Edge=dict(seed=42),
                visualization=dict(defaults=dict(save=False, show=False)))


class LayoutTest(unittest.TestCase):
    def test_fc_layout_invariance_and_reuse(self):
        results=[]
        signal=np.random.default_rng(123).normal(size=(30,2,2))
        for layout, filename in [('week_first','M1_R1_W0'),('run_first','M1_R1_W0')]:
            with self.subTest(layout=layout), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                ctx=RunContext(None,config(root),logging.getLogger('layouts'))
                mask=root/'mask'; mask.mkdir()
                Image.fromarray(np.array([[1,0],[1,0]],dtype=np.uint8)).save(mask/'SM-test.png')
                folder=root/'processed'/('W0/M1/R1' if layout=='week_first' else 'R1/M1/W0')
                folder.mkdir(parents=True)
                np.save(folder/(filename+'_FiltDRR.npy'),signal)
                cache=CorrelationCache(ctx,['AUDv'])
                frame,status=cache.process('W0','M1','R1')
                self.assertEqual(status,'success')
                self.assertEqual(cache.process('W0','M1','R1')[1],'reused')
                results.append(frame.to_numpy())
        for result in results[1:]:
            np.testing.assert_array_equal(results[0],result)

    def test_old_names_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=RunContext(None,config(root),logging.getLogger('layouts'))
            folder=root/'processed/W0/M1/R1';folder.mkdir(parents=True)
            for name in ('M1_R1_FiltDRR.npy','M1_W0_FiltDRR.npy'):
                np.save(folder/name,np.ones((3,2,2)))
            with self.assertRaisesRegex(FileNotFoundError,'M1_R1_W0_FiltDRR.npy'):
                resolve_recording(ctx,'W0','M1','R1')
            correct=folder/'M1_R1_W0_FiltDRR.npy'
            np.save(correct,np.ones((3,2,2)))
            self.assertEqual(resolve_recording(ctx,'W0','M1','R1'),correct.resolve())

    def test_duplicates_missing_and_raw_layouts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=RunContext(None,config(root),logging.getLogger('layouts'))
            with self.assertRaises(FileNotFoundError): resolve_recording(ctx,'W0','M1','R1')
            for relative in ('W0/M1/R1','R1/M1/W0'):
                folder=root/'processed'/relative;folder.mkdir(parents=True)
                np.save(folder/'M1_R1_W0_FiltDRR.npy',np.ones((3,2,2)))
            with self.assertRaisesRegex(ValueError,'Ambiguous'): resolve_recording(ctx,'W0','M1','R1')
            folder=root/'raw/W0/R1/M1';folder.mkdir(parents=True)
            (folder/'info.txt').write_text('metadata')
            self.assertEqual(resolve_recording(ctx,'W0','M1','R1',raw=True),folder.resolve())

    def test_stage02_canonical_reports_tables_and_duplicate_continuation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg=config(root)
            cfg['variables']['Week_ids']=['W0','W1']
            cfg['variables']['Mouse_ids']=['M1']
            ctx=RunContext(None,cfg,logging.getLogger('stage02'))
            mask=root/'mask';mask.mkdir()
            Image.fromarray(np.array([[1,0],[1,0]],dtype=np.uint8)).save(mask/'SM-test.png')
            for rel in ('W0/M1/R1','R2/M1/W0','W1/M1/R1','R1/M1/W1'):
                folder=root/'processed'/rel;folder.mkdir(parents=True)
                run='R2' if rel.startswith('R2') else 'R1'
                week='W0' if 'W0' in rel else 'W1'
                np.save(folder/f'M1_{run}_{week}_FiltDRR.npy',np.random.default_rng(5).normal(size=(30,2,2)))
            with patch('srcs.utils.config.load_context',return_value=ctx):
                runpy.run_path(str(ROOT/'scripts/02_correlation.py'),run_name='__main__')
            frame=pd.read_pickle(Path(ctx.correlation_path)/'W0_correlation_matrices_AUDv_excluded.pkl')
            self.assertEqual(list(frame.index),['M1'])
            self.assertEqual(list(frame.columns),['R1','R2'])
            self.assertEqual(frame.attrs['schema'],'week_mouse_run_v1')
            self.assertTrue(frame.attrs['complete'])
            report=json.loads(next(Path(ctx.output_path).glob('BatchReports/02*.json')).read_text())
            self.assertEqual(report['counts'],dict(success=2,missing=1,failed=1))
            self.assertTrue(all({'week','mouse','run','atlas'}.issubset(row) for row in report['records']))


class AxisTest(unittest.TestCase):
    def test_one_condition_statistics_and_explicit_unavailable_supplement(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=RunContext(None,config(Path(tmp)),logging.getLogger('one-condition'))
            ctx.config['variables']['Mouse_ids']=['M1','M2','M3','M4']
            ctx.config['variables']['Atlas']=[['A'],['B'],['C']]
            rng=np.random.default_rng(7)
            stat=pd.DataFrame(index=['A','B','C'],columns=['W0'],dtype=object)
            for a in stat.index:
                group=rng.normal(.2,.05,4)
                indiv=rng.normal(.5,.05,4)
                stat.at[a,'W0']=[group,indiv,indiv-group,
                    dict(diagonal_dominance=float(rng.random()),
                         normalized_distance_to_diagonal=float(rng.random())),1.]
            stat.attrs['analysis']=dict(split_axis='run',selection_mode='single',
                conditions={'W0':dict(weeks=['W0'],runs=ctx.runs)},mice=ctx.mice)
            Path(ctx.statistics_path).parent.mkdir(parents=True)
            stat.to_pickle(ctx.statistics_path)
            class PlotStub:
                def __init__(self,*args,**kwargs):pass
                def __getattr__(self,name):return lambda *args,**kwargs:self
            with patch('srcs.utils.config.load_context',return_value=ctx),patch('srcs.vis.PlotSeries.PlotSeries',PlotStub):
                runpy.run_path(str(ROOT/'scripts/04_statistics.py'),run_name='__main__')
                runpy.run_path(str(ROOT/'scripts/05_supplement.py'),run_name='__main__')
            folder=Path(ctx.output_path)/'Statistics/Tables'
            self.assertTrue((folder/'anova_A.csv').is_file())
            status=json.loads((folder/'Supplement/inference_status.json').read_text())
            self.assertEqual(status['status'],'unavailable')
            self.assertFalse((folder/'Supplement/manova.pkl').exists())

    def tables(self,ctx):
        result={}
        for wi,w in enumerate(ctx.weeks):
            f=pd.DataFrame(index=ctx.mice,columns=ctx.runs,dtype=object)
            for mi,m in enumerate(ctx.mice):
                for ri,r in enumerate(ctx.runs):
                    v=.1+wi*.03+mi*.01+ri*.02
                    f.at[m,r]=pd.DataFrame([[1,v],[v,1]],index=['L','R'],columns=['L','R'])
            result[w]=f
        return result

    def test_old_run_grouped_identification_is_numerically_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=RunContext(None,config(Path(tmp)),logging.getLogger('axes'))
            data=self.tables(ctx)
            for mode in ('single','cumulative'):
                ctx.config['analysis']['selection_mode']=mode
                for choice in analysis_plan(ctx,ctx.weeks).values():
                    weeks=choice['weeks']
                    # Former layout: session files are Run, columns are week.
                    old=[]
                    for r in ctx.runs:
                        f=pd.DataFrame(index=ctx.mice,columns=ctx.weeks,dtype=object)
                        for w in ctx.weeks:
                            for m in ctx.mice:f.at[m,w]=data[w].at[m,r]
                        old.append(f)
                    expected=build_half_dataframe(old,0,'random',seed=42,trial_ids=weeks)
                    new=axis_frames(data,'run',ctx.mice,weeks,ctx.runs)
                    actual=build_half_dataframe(list(new.values()),0,'random',seed=42)
                    for h in actual.index:
                        for m in ctx.mice:np.testing.assert_array_equal(actual.at[h,m],expected.at[h,m])

    def test_edge_uses_shared_groups_and_single_opposite_axis_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=RunContext(None,config(Path(tmp)),logging.getLogger('axes'))
            data=self.tables(ctx)
            for axis in ('week','run'):
                ctx.config['analysis']['split_axis']=axis
                for choice in analysis_plan(ctx,ctx.weeks).values():
                    weeks,runs=choice['weeks'],choice['runs']
                    frames=axis_frames(data,'week' if axis=='run' else 'run',ctx.mice,weeks,runs)
                    groups=ctx.group_indices(ctx.weeks,ctx.runs)
                    edge=Edgewise_Contribution('',ctx.mice,list(frames),atlas_exclude=['AUDv'],frames=frames,group_indices=groups)
                    phi,dp=edge.run(repeat=50)
                    self.assertTrue(phi and dp)
                    pairs=edge.build_pairs()
                    for fixed,subjects in pairs.items():
                        for mouse,(first,second) in subjects.items():
                            matrices=edge.data[fixed][mouse]
                            np.testing.assert_array_equal(first,np.mean(matrices[groups[0]],axis=0))
                            np.testing.assert_array_equal(second,np.mean(matrices[groups[1]],axis=0))

    def test_config_isolation_and_metadata_mismatch(self):
        ctx=RunContext(None,config(Path('relative')),logging.getLogger('axes'))
        first=ctx.output_path
        ctx.config['analysis']['split_axis']='week'
        self.assertNotEqual(ctx.output_path,first)
        meta=dict(split_axis='run',selection_mode='single',conditions={'W0':{}})
        stat=pd.DataFrame(columns=['W0']);stat.attrs['analysis']=meta
        with self.assertRaisesRegex(ValueError,'axes'): statistics_metadata(ctx,stat)
        ctx.config['analysis']['split_axis']='invalid'
        with self.assertRaises(ValueError):_=ctx.split_axis
        ctx.config['analysis']['selection_mode']='invalid'
        with self.assertRaises(ValueError):_=ctx.selection_mode
