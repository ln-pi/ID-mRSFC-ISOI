"""Check every merged requirement and import core/project/GUI dependencies."""
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import yaml

def main():
    cfg = yaml.safe_load((ROOT / "tests/environment.yaml").read_text(encoding="utf-8"))
    specs = next(d["pip"] for d in cfg["dependencies"] if isinstance(d, dict))
    versions = {}
    for spec in specs:
        name = re.split(r"[<>=!~\[]", spec)[0]
        versions[name] = metadata.version(name)
    modules = [
        "numpy", "scipy", "pandas", "matplotlib", "yaml", "tqdm", "seaborn",
        "sklearn", "statsmodels", "networkx", "igraph", "PIL", "imageio",
        "tifffile", "ipykernel", "pingouin", "pyvis", "cdlib", "dynetx",
        "hmmlearn", "cv2", "numba", "torch", "suite2p", "suite2p.gui.gui2p", "PyQt6.QtWidgets",
        "srcs.prep.preprocess", "srcs.prep.suite2p", "srcs.comp.ID_Tools",
    ]
    imports = {}
    for module in modules:
        print("Import:", module, flush=True)
        try:
            importlib.import_module(module)
            imports[module] = "ok"
        except Exception as exc:
            imports[module] = repr(exc)
    if imports.get("PyQt6.QtWidgets") == "ok":
        from PyQt6.QtWidgets import QApplication
        app = QApplication.instance() or QApplication([])
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        canvas = FigureCanvasQTAgg(Figure())
        canvas.figure.subplots().plot([0, 1], [0, 1])
        canvas.draw()
        imports["Qt6_matplotlib_canvas"] = "ok"
        canvas.close()
    import suite2p
    report = {
        "python": platform.python_version(), "executable": sys.executable,
        "versions": versions, "imports": imports,
        "suite2p_settings_keys": {key: list(value) if isinstance(value, dict) else value
                                  for key, value in suite2p.default_settings().items()},
    }
    target = ROOT / "reports/environment_validation/imports.json"
    target.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    assert platform.python_version() == "3.13.13"
    assert all(value == "ok" for value in imports.values()), imports
    print("PASS: declared packages, project imports and offscreen Qt6 canvas", flush=True)

if __name__ == "__main__":
    main()
