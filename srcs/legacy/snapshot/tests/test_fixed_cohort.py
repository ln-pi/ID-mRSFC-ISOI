import json
import logging
from pathlib import Path
import tempfile
import unittest
import pandas as pd
from srcs.utils.config import RunContext
from srcs.utils.dataset import select_complete_cohort
from srcs.utils.results import save_json

class CohortTest(unittest.TestCase):
    def fixture(self, root, mode, mice, weeks, missing):
        cfg={'paths':{k:str(root/k) for k in ('raw','processed','output','mask')},
             'variables':{'dataset':'','Week_ids':weeks,'Mouse_ids':mice,
                          'Run_ids':['R1','R2'],'Atlas':[['A'],['B']]},
             'ID':{'session_seed':42,'group_mode':'random','completeness_mode':mode},'Edge':{'seed':42}}
        ctx=RunContext(None,cfg,logging.getLogger('cover'))
        Path(ctx.correlation_path).mkdir(parents=True)
        records=[]
        for w in ctx.weeks:
            for a in ('A','B'):
                frame=pd.DataFrame(index=mice,columns=ctx.runs,dtype=object)
                for m in mice:
                    for r in ctx.runs:
                        bad=(w,m,r,a) in missing
                        frame.at[m,r]=None if bad else pd.DataFrame([[1,.2],[.2,1]],index=['L','R'],columns=['L','R'])
                        records.append(dict(week=w,mouse=m,run=r,atlas=a,status='missing' if bad else 'success'))
                frame.to_pickle(Path(ctx.correlation_path)/f'{w}_correlation_matrices_{a}_excluded.pkl')
        save_json(Path(ctx.output_dir('BatchReports'))/'02_correlation_test.json',{'records':records})
        return ctx

    def test_real_pattern_and_config_isolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            mice=[f'A{i:02}' for i in range(1,11)]
            weeks=[f'{i:02}w' for i in range(14)]
            missing={('01w','A01','R1','A'),('13w','A06','R1','B'),('13w','A06','R2','B')}
            strict=self.fixture(Path(tmp),'strict',mice,weeks,missing)
            cover=self.fixture(Path(tmp),'cover',mice,weeks,missing)
            self.assertNotEqual(strict.output_path,cover.output_path)
            before=json.dumps(cover.config,sort_keys=True)
            kept,audit=select_complete_cohort(strict)
            self.assertEqual(kept,[m for m in mice if m not in ('A01','A06')])
            self.assertEqual(audit['retained_weeks'],weeks)
            kept,audit=select_complete_cohort(cover)
            self.assertEqual(kept,[m for m in mice if m!='A01'])
            self.assertEqual(audit['retained_weeks'],weeks[:-1])
            self.assertEqual(audit['selected_score'],117)
            self.assertEqual(len(audit['candidates']),105)
            self.assertEqual(before,json.dumps(cover.config,sort_keys=True))

    def test_ties_and_noninitial_interval(self):
        cases=[(['a','b','c','d'],['W1','W2'],{('W2',m,'R1','A') for m in ['c','d']},['W1']),
               (['a','b'],['W1','W2','W3'],{('W2','a','R2','B')},['W1']),
               (['a','b'],['W1','W2','W3'],{('W1','a','R2','B')},['W2','W3']),
               (['a','b'],['W1'],set(),['W1']),
               (['a','b'],['W1'],{('W1','a','R2','B')},[])]
        for mice,weeks,missing,expected in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as tmp:
                ctx=self.fixture(Path(tmp),'cover',mice,weeks,missing)
                audit=select_complete_cohort(ctx)[1]
                self.assertEqual(audit['retained_weeks'],expected)

    def test_missing_table_and_invalid_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.fixture(Path(tmp),'cover',['a','b'],['W1'],set())
            (Path(ctx.correlation_path)/'W1_correlation_matrices_B_excluded.pkl').unlink()
            with self.assertRaises(FileNotFoundError): select_complete_cohort(ctx)
            ctx.config['ID']['completeness_mode']='invalid'
            with self.assertRaisesRegex(ValueError,'completeness_mode'): select_complete_cohort(ctx)

    def test_report_disagreement_incomplete_and_invalid_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.fixture(Path(tmp),'strict',['a','b','c'],['W1','W2'],{('W1','a','R1','A')})
            report=Path(ctx.output_dir('BatchReports'))/'02_correlation_test.json'
            data=json.loads(report.read_text())
            for row in data['records']: row['status']='success'
            save_json(report,data)
            self.assertEqual(select_complete_cohort(ctx)[0],['b','c'])
            path=Path(ctx.correlation_path)/'W2_correlation_matrices_A_excluded.pkl'
            frame=pd.read_pickle(path)
            frame.at['b','R1'].iloc[0,1]=float('inf')
            frame.to_pickle(path)
            self.assertEqual(select_complete_cohort(ctx)[0],[])
            data['records'].pop()
            save_json(report,data)
            with self.assertRaisesRegex(ValueError,'incomplete'): select_complete_cohort(ctx)
