import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm

from srcs.utils.config import load_context
from srcs.utils.batch import BatchReport
from srcs.utils.dataset import resolve_recording
from srcs.utils.results import save_pickle
from srcs.comp.correlation import CorrelationCache


def run(ctx):
    if (Path(ctx.output_path) / 'metadata/fc_import_provenance.json').exists():
        raise ValueError('Imported FC output: run identification, or use a separate config to recompute FC')
    report = BatchReport(ctx, "02_correlation")
    total = len(ctx.atlases) * len(ctx.weeks) * len(ctx.mice) * len(ctx.runs)
    with tqdm(total=total, desc="Correlating", leave=False) as pbar:
        try:
            for atlas in ctx.atlases:
                atlas_tag = "_".join(atlas)
                cache = CorrelationCache(ctx, atlas)
                for week in ctx.weeks:
                    target = Path(ctx.correlation_path) / f"{week}_correlation_matrices_{atlas_tag}_excluded.pkl"
                    df = pd.DataFrame(None, index=ctx.mice, columns=ctx.runs, dtype=object)
                    df[:] = None
                    for mouse in ctx.mice:
                        for run in ctx.runs:
                            pbar.set_postfix(Week=week, Mouse=mouse, Run=run, Atlas=atlas_tag)
                            def correlate():
                                source = resolve_recording(ctx, week, mouse, run)
                                frame, status = cache.process(week, mouse, run)
                                frame.attrs["recording"] = dict(week=week, mouse=mouse, run=run, input_path=str(source))
                                df.at[mouse, run] = frame
                                return dict(status=status, input_path=str(source))
                            report.run(correlate, week=week, mouse=mouse, run=run,
                                atlas=atlas_tag, input_root=ctx.processed_path)
                            pbar.update(1)
                    df.attrs["complete"] = all(isinstance(df.at[m, t], pd.DataFrame)
                                               for m in ctx.mice for t in ctx.runs)
                    df.attrs.update(schema="week_mouse_run_v1", week=week, weeks=ctx.weeks, runs=ctx.runs, mice=ctx.mice)
                    df.attrs["batch_report"] = str(report.path)
                    save_pickle(target, df)
                    ctx.logger.info("FC table saved: %s; complete=%s", target, df.attrs["complete"])
        finally:
            report.finish()


    if not report.persisted:
        raise OSError("FC report could not be saved")
    if any(row["status"] == "failed" for row in report.rows):
        raise RuntimeError(f"FC contains failed records; see {report.path}")
