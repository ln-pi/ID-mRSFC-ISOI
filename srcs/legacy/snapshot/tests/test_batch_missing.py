import logging
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch
import json
import numpy as np
import pandas as pd
from PIL import Image
from srcs.utils.config import RunContext

ROOT = Path(__file__).resolve().parents[1]

class MissingBatchTest(unittest.TestCase):
    def test_both_stages_continue_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = {"paths": {k: str(root/k) for k in ("raw", "processed", "output", "mask")},
                   "variables": {"dataset": "", "Week_ids": ["S1", "S2"],
                   "Mouse_ids": ["M"], "Run_ids": ["missing", "broken", "good"], "Atlas": [["None"]]},
                   "ID": {"group_mode": "random", "session_seed": 42}, "Edge": {"seed": 42}}
            cfg["visualization"] = {"defaults": {"save": False, "show": False}}
            ctx = RunContext(None, cfg, logging.getLogger("batch-test"))
            mod = runpy.run_path(str(ROOT/'scripts/01_preprocess.py'))
            visited = []
            def process(ctx, session, mouse, trial):
                visited.append((session, trial))
                if trial == "missing": raise FileNotFoundError("absent input")
                if trial == "broken": raise ValueError("corrupt input")
            with patch.dict(mod['main'].__globals__, {"load_context": lambda: ctx, "process_recording": process}):
                mod['main']()
            self.assertEqual(len(visited), 6)
            mask = root/'mask'; mask.mkdir()
            Image.fromarray(np.array([[1, 0], [1, 0]], dtype=np.uint8)).save(mask/'SM-test.png')
            for session in ctx.weeks:
                base = root/'processed'/session/'M'
                (base/'broken').mkdir(parents=True)
                (base/f'broken/M_broken_{session}_FiltDRR.npy').write_bytes(b'broken')
                (base/'good').mkdir()
                np.save(base/f'good/M_good_{session}_FiltDRR.npy', np.random.default_rng(1).normal(size=(20, 2, 2)))
            with patch('srcs.utils.config.load_context', return_value=ctx), patch('srcs.vis.PlotSeries.save_fc_heatmap'):
                runpy.run_path(str(ROOT/'scripts/02_correlation.py'), run_name='__main__')
            for stage in ('01_preprocess', '02_correlation'):
                report = next(Path(ctx.output_path).glob(f'BatchReports/{stage}_*.json'))
                data = json.loads(report.read_text(encoding='utf-8'))
                self.assertEqual(data['counts'], {"success": 2, "missing": 2, "failed": 2})
            for session in ctx.weeks:
                frame = pd.read_pickle(Path(ctx.correlation_path)/f'{session}_correlation_matrices_None_excluded.pkl')
                self.assertIsNone(frame.at['M', 'missing'])
                self.assertIsNone(frame.at['M', 'broken'])
                self.assertIsInstance(frame.at['M', 'good'], pd.DataFrame)
                self.assertFalse(frame.attrs['complete'])
