import json
import logging
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from srcs.utils import results
from srcs.utils.batch import BatchReport

class OutputPermissionTest(unittest.TestCase):
    def test_transient_replace_and_persistent_denial(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"report.json"
            results.save_json(path,{"old":True})
            real=results.os.replace
            count=0
            def replace(src,dst):
                nonlocal count
                count+=1
                if count<3: raise PermissionError("locked")
                return real(src,dst)
            with patch.object(results.os,'replace',side_effect=replace), patch.object(results.time,'sleep'):
                results.save_json(path,{"new":True})
            self.assertEqual(json.loads(path.read_text()),{"new":True})
            with patch.object(results.os,'replace',side_effect=PermissionError('denied')), patch.object(results.time,'sleep'):
                with self.assertRaises(PermissionError): results.save_json(path,{"lost":True})
            self.assertEqual(json.loads(path.read_text()),{"new":True})
            self.assertEqual(len(list(Path(tmp).iterdir())),1)

    def test_report_fallback_and_later_recovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=SimpleNamespace(logger=logging.getLogger('permissions'),output_dir=lambda *a:tmp)
            report=BatchReport(ctx,'test')
            original=report.path
            real=results.save_json
            def save(path,value):
                if path==original: raise PermissionError('destination locked')
                return real(path,value)
            with patch('srcs.utils.batch.save_json',side_effect=save):
                self.assertTrue(report.run(lambda:None,mouse='A01'))
            self.assertNotEqual(report.path,original)
            self.assertEqual(json.loads(report.path.read_text())['counts']['success'],1)
            with patch('srcs.utils.batch.save_json',side_effect=PermissionError('directory denied')):
                self.assertTrue(report.run(lambda:None,mouse='A02'))
                self.assertFalse(report.persisted)
            report.finish()
            self.assertTrue(report.persisted)
            self.assertEqual(len(json.loads(report.path.read_text())['records']),2)
