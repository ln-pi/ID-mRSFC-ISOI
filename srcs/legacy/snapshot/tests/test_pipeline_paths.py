import ast
import importlib.util
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

from srcs.prep.suite2p import run_isolated_process
from srcs.utils.utils import extract_roi, load_chamber_mask


class SharedMasksTest(unittest.TestCase):
    def test_shared_rgba_mask_is_not_a_roi(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rgba = np.full((3, 4, 4), 255, dtype=np.uint8)
            rgba[0, 0, 3] = 0
            Image.fromarray(rgba).save(root / "Chamber_mask.png")
            roi = np.zeros((3, 4), dtype=np.uint8)
            roi[:, 0] = 255
            Image.fromarray(roi).save(root / "AUD-test.png")
            mask = load_chamber_mask(folder, (3, 4))
            self.assertFalse(mask[0, 0])
            self.assertEqual(mask.sum(), 11)
            data = np.arange(24).reshape(2, 3, 4).astype(float)
            series, names = extract_roi(data, folder)
            self.assertEqual(names, ["AUD-test_L", "AUD-test_R"])
            np.testing.assert_allclose(series[:, 0], data[:, :, 0].mean(axis=1))
            np.testing.assert_allclose(series[:, 1], data[:, :, -1].mean(axis=1))
            with self.assertRaisesRegex(ValueError, "does not match"):
                load_chamber_mask(folder, (4, 3))
            with self.assertRaisesRegex(ValueError, "does not match"):
                extract_roi(np.ones((2, 4, 3)), folder)

    def test_missing_empty_and_chamber_only(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(FileNotFoundError):
                load_chamber_mask(folder)
            Image.fromarray(np.zeros((3, 4), dtype=np.uint8)).save(Path(folder) / "Chamber_Mask.png")
            with self.assertRaisesRegex(ValueError, "no valid pixels"):
                load_chamber_mask(folder)
            with self.assertRaises(FileNotFoundError):
                extract_roi(np.ones((2, 3, 4)), folder)


class Suite2pPathsTest(unittest.TestCase):
    def test_session_isolation_and_explicit_input(self):
        captured = []
        temporary_scripts = []
        def fake_run(command, **kwargs):
            script = Path(command[1])
            temporary_scripts.append(script)
            source = script.read_text(encoding="utf-8")
            compile(source, str(script), "exec")
            tree = ast.parse(source)
            db_node = next(node.value for node in tree.body if isinstance(node, ast.Assign)
                           and any(isinstance(t, ast.Name) and t.id == "db" for t in node.targets))
            db = eval(compile(ast.Expression(db_node), "<db>", "eval"), {"os": os})
            captured.append(db)
            self.assertTrue(db["keep_movie_raw"], "Two-step registration needs raw frames")
            if importlib.util.find_spec("suite2p") is not None:
                from suite2p import default_settings
                schema = default_settings()
                for node in tree.body:
                    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
                        continue
                    call = node.value
                    if not isinstance(call.func, ast.Attribute) or call.func.attr != "update":
                        continue
                    target = call.func.value
                    if isinstance(target, ast.Name) and target.id == "settings":
                        known = schema
                    elif isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name) and target.value.id == "settings":
                        known = schema[ast.literal_eval(target.slice)]
                    else:
                        continue
                    keys = [ast.literal_eval(key) for key in call.args[0].keys]
                    self.assertEqual(set(keys) - set(known), set(), f"Unsupported Suite2p keys: {keys}")
            inputs = [Path(db["data_path"][0]) / name for name in db["file_list"]]
            self.assertEqual([p.name for p in inputs], ["A01_00w.tif"])
            self.assertTrue(all(p.is_file() for p in inputs))
            self.assertFalse(db["look_one_level_down"])
            return SimpleNamespace(returncode=0, stderr="")
        with tempfile.TemporaryDirectory() as folder:
            for session in ("Run01", "Run02"):
                recording = Path(folder) / session / "A01" / "00w"
                recording.mkdir(parents=True)
                original = recording / "A01_00w.tif"
                original.touch()
                (recording / "A01_00w_merged.tif").touch()
                with patch("srcs.prep.suite2p.subprocess.run", side_effect=fake_run):
                    run_isolated_process(str(original), "A01", "00w", str(recording))
                self.assertEqual(Path(captured[-1]["save_path0"]), recording)
                self.assertEqual(captured[-1]["save_folder"], "suite2p")
        self.assertNotEqual(captured[0]["save_path0"], captured[1]["save_path0"])
        self.assertTrue(all(not p.exists() for p in temporary_scripts))


if __name__ == "__main__":
    unittest.main()
