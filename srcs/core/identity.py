import ast
import hashlib
import json
from importlib.util import resolve_name
from pathlib import Path


EXCLUDED = ("srcs.visualization", "srcs.core", "srcs.utils.logger",
            "srcs.utils.provenance", "srcs.utils.config", "srcs.utils.dataset",
            "srcs.utils.batch")


def digest(value):
    serialized = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(serialized.encode()).hexdigest()


def module_path(root, module):
    base = Path(root).joinpath(*module.split("."))
    for path in (base.with_suffix(".py"), base / "__init__.py"):
        if path.is_file():
            return path
    return None


def code_closure(root, module):
    pending = [module]
    hashes = {}
    while pending:
        name = pending.pop()
        if name in hashes or any(name == p or name.startswith(p + ".") for p in EXCLUDED):
            continue
        path = module_path(root, name)
        if path is None:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        package = name if path.name == "__init__.py" else name.rpartition(".")[0]
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                if ast.get_docstring(node, clean=False) is not None:
                    del node.body[0]
            if isinstance(node, ast.Import):
                pending.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = resolve_name("." * node.level + (node.module or ""), package)
                pending.append(base)
                pending.extend(base + "." + alias.name for alias in node.names if alias.name != "*")
        hashes[name] = digest(ast.dump(tree, include_attributes=False))
    if module not in hashes:
        raise ValueError(f"No numerical module found: {module}")
    return hashes
