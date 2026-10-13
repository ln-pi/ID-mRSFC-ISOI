import logging
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from PIL import Image
import pandas as pd
import yaml
from srcs.comp.correlation import CorrelationCache
from srcs.utils.config import RunContext

ROOT = Path(__file__).resolve().parents[1]

class CorrelationResumeTest(unittest.TestCase):
    def setup_case(self, root):
        cfg = {"paths": {k: str(root/k) for k in ("raw", "processed", "output", "mask")},
               "variables": {"dataset": "", "Week_ids": ["S1", "S2"], "Mouse_ids": ["M"],
                             "Run_ids": ["T"], "Atlas": [["AUDv"]]},
               "ID": {"group_mode": "random", "session_seed": 42}, "Edge": {"seed": 42},
               "visualization": {"defaults": {"save": True, "show": False}}}
        ctx = RunContext(None, cfg, logging.getLogger("resume-test"))
        mask = root/'mask'; mask.mkdir()
        Image.fromarray(np.array([[1,0],[1,0]], dtype=np.uint8)).save(mask/'SM-test.png')
        source = root/'processed/S1/M/T/M_T_S1_FiltDRR.npy'
        source.parent.mkdir(parents=True)
        np.save(source, np.random.default_rng(1).normal(size=(20,2,2)))
        return ctx, source

    def test_reuse_invalidation_and_output_repair(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx, source = self.setup_case(Path(tmp))
            def plot(data, **kwargs):
                Image.new('RGB',(3,3)).save(Path(kwargs['save_path'])/(kwargs['name']+'.png'))
            with patch('srcs.comp.correlation.save_fc_heatmap', side_effect=plot) as draw:
                cache = CorrelationCache(ctx, ['AUDv'])
                first,status = cache.process('S1','M','T')
                self.assertEqual(status,'success')
                with patch('srcs.comp.correlation.extract_roi', side_effect=AssertionError('must reuse')):
                    self.assertEqual(cache.process('S1','M','T')[1],'reused')
                    self.assertEqual(draw.call_count,1)
                    numeric=next(Path(ctx.correlation_path).glob('FC_Matrix/*.npy'))
                    numeric.write_bytes(b'broken')
                    self.assertEqual(cache.process('S1','M','T')[1],'reused')
                    np.testing.assert_array_equal(np.load(numeric), first.to_numpy())
                np.save(source, np.random.default_rng(4).normal(size=(20,2,2)))
                self.assertEqual(cache.process('S1','M','T')[1],'success')
                Image.fromarray(np.array([[1,0],[0,0]],dtype=np.uint8)).save(Path(ctx.mask_root)/'SM-test.png')
                cache=CorrelationCache(ctx,['AUDv'])
                self.assertEqual(cache.process('S1','M','T')[1],'success')
                ctx.config['correlation']={'resume':False}
                self.assertEqual(CorrelationCache(ctx,['AUDv']).process('S1','M','T')[1],'success')
                source.unlink()
                with self.assertRaises(FileNotFoundError): cache.process('S1','M','T')

    def test_config_isolation_no_legacy(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx, source = self.setup_case(Path(tmp))
            ctx.config['visualization']['defaults']['save']=False
            cache=CorrelationCache(ctx,['AUDv'])
            legacy=pd.DataFrame([[1,.2],[.2,1]],index=cache.names,columns=cache.names)
            legacy.to_pickle(Path(ctx.processed_path)/'S1_correlation_matrices_AUDv_excluded.pkl')
            self.assertEqual(cache.process('S1','M','T')[1],'success')
            first_path=ctx.correlation_path
            self.assertEqual(cache.process('S1','M','T')[1],'reused')
            ctx.config['ID']['session_seed']=43
            self.assertNotEqual(first_path,ctx.correlation_path)
            self.assertEqual(cache.process('S1','M','T')[1],'success')
            ctx.config['ID']['session_seed']=42
            self.assertEqual(first_path,ctx.correlation_path)
            self.assertEqual(cache.process('S1','M','T')[1],'reused')
            cached=next(Path(ctx.correlation_path).glob('FC_Matrix/*.cache.pkl'))
            cached.write_bytes(b'broken')
            self.assertEqual(cache.process('S1','M','T')[1],'success')
            ctx.config['visualization']['fc_heatmap']={'cmap':'viridis'}
            self.assertNotEqual(first_path,ctx.correlation_path)
            self.assertEqual(cache.process('S1','M','T')[1],'success')

    def test_three_atlas_configs(self):
        for name in ('OIS','default'):
            cfg=yaml.safe_load((ROOT/f'configs/{name}.yaml').read_text(encoding='utf-8'))
            atlases=cfg['variables']['Atlas']
            self.assertEqual(atlases, [['SM','VIS','AUD'],['PF','LAT','MED','AUDv'],['AUDv']])
            self.assertEqual(cfg['supplement']['structure_atlases'],atlases)
            self.assertEqual(cfg['supplement']['exclude'],[])
