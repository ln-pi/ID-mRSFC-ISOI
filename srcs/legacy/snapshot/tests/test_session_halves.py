import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.utils.utils import build_half_dataframe, split_session_indices, fisher_z


class SessionHalvesTest(unittest.TestCase):
    def test_group_invariants(self):
        for n in range(2, 32):
            for mode in ("sequential", "random"):
                a, b = split_session_indices(n, mode, seed=42)
                self.assertEqual(len(a), (n + 1) // 2)
                self.assertEqual(len(a), len(b))
                self.assertEqual(set(a) | set(b), set(range(n)))
                self.assertEqual(len(set(a) & set(b)), n % 2)
                if n % 2:
                    middle = n // 2 if mode == "sequential" else np.random.default_rng(42).permutation(n)[n // 2]
                    self.assertEqual(set(a) & set(b), {middle})

    def test_explicit_repeat(self):
        self.assertEqual(split_session_indices(5), ([0, 1, 2], [2, 3, 4]))
        self.assertEqual(split_session_indices(5, repeat_session=2), ([0, 1, 2], [1, 3, 4]))
        for n in (3, 5, 9):
            for repeat in range(1, n + 1):
                a, b = split_session_indices(n, repeat_session=repeat)
                self.assertEqual(set(a) & set(b), {repeat - 1})
                self.assertEqual(len(a), len(b))
                self.assertEqual(set(a) | set(b), set(range(n)))

    def test_seed_and_legacy(self):
        self.assertEqual(split_session_indices(9, "random", seed=7), split_session_indices(9, "random", seed=7))
        self.assertEqual(split_session_indices(4, 12), ([0, 1], [2, 3]))
        self.assertEqual(split_session_indices(4, 13), ([0, 2], [1, 3]))
        self.assertEqual(split_session_indices(4, 14), ([0, 3], [1, 2]))

    def test_invalid_groups(self):
        for args in [(1,), (2.5,), (True,), (2, "bad"), (3, 14), (4, "sequential", 1), (3, "random", 1), (3, "sequential", 0), (3, "sequential", 4)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                split_session_indices(*args)

    def test_averaging_and_sizes(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = []
            for i in range(6):
                frame = pd.DataFrame(index=["M1", "M2", "M3"], columns=range(8), dtype=object)
                for mouse in frame.index:
                    for trial in frame.columns:
                        r = (i * 8 + trial + 1) / 100
                        frame.at[mouse, trial] = pd.DataFrame([[1, r], [r, 1]], index=["a", "b"], columns=["a", "b"])
                path = Path(folder) / f"{i}.pkl"
                frame.to_pickle(path)
                paths.append(path)
            for n, mode in [(2, "sequential"), (5, "sequential"), (5, "random"), (6, "random"), (4, 12), (4, 13), (4, 14)]:
                for count in (0, 1, 8):
                    result = build_half_dataframe(paths[:n], count, mode, seed=42)
                    self.assertEqual(result.shape, (2, 3))
                    groups = split_session_indices(n, mode, seed=42)
                    for name, group in zip(["1st", "2nd"], groups):
                        expected = np.mean([np.arctanh((i * 8 + t + 1) / 100) for i in group for t in range(count or 8)])
                        self.assertAlmostEqual(result.at[name, "M2"].loc["a", "b"], expected)
                        self.assertEqual(result.at[name, "M2"].index.tolist(), ["a", "b"])
            for count in (-1, 9, 1.5):
                with self.assertRaises(ValueError):
                    build_half_dataframe(paths[:2], count)
            frame = pd.read_pickle(paths[1]).rename(index={"M1": "other"})
            frame.to_pickle(paths[1])
            with self.assertRaisesRegex(ValueError, "identical mouse"):
                build_half_dataframe(paths[:2])


if __name__ == "__main__":
    unittest.main()
