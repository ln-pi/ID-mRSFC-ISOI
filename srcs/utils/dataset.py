"""Canonical week x mouse x run recording and FC table paths."""
from itertools import permutations
from pathlib import Path


def recording_folders(root, week, mouse, run):
    return list(dict.fromkeys(Path(root).joinpath(*parts)
                             for parts in permutations((week, mouse, run))))


def resolve_recording(ctx, week, mouse, run, raw=False):
    root = Path(ctx.raw_root) / ctx.dataset if raw else Path(ctx.processed_path)
    folders = recording_folders(root, week, mouse, run)
    if raw:
        matches = [p for p in folders if (p / "info.txt").is_file() or any(p.glob("img*.bin"))]
    else:
        name = f"{mouse}_{run}_{week}_FiltDRR.npy"
        matches = [p / name for p in folders if (p / name).is_file()]
    matches = list(dict.fromkeys(p.resolve() for p in matches))
    if len(matches) > 1:
        raise ValueError(f"Ambiguous recording {(week, mouse, run)}: {matches}")
    if not matches:
        raise FileNotFoundError(f"Missing recording {(week, mouse, run)} under {root}; "
                                + ("expected info.txt/img*.bin" if raw else f"expected {mouse}_{run}_{week}_FiltDRR.npy"))
    return matches[0]


def preprocessing_folder(ctx, week, mouse, run):
    folders = [p for p in recording_folders(ctx.processed_path, week, mouse, run) if p.is_dir()]
    if len(folders) > 1:
        raise ValueError(f"Multiple processed recording folders for {(week, mouse, run)}: {folders}")
    return folders[0] if folders else Path(ctx.processed_dir(week, mouse, run))


def table_path(ctx, week, atlas):
    return Path(ctx.correlation_path) / f"{week}_correlation_matrices_{atlas}_excluded.pkl"
