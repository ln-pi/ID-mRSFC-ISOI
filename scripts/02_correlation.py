import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm

from srcs.utils.config import load_context
from srcs.utils.batch import BatchReport
from srcs.utils.results import save_pickle
from srcs.comp.correlation import CorrelationCache


def main():
    ctx = load_context()
    report = BatchReport(ctx, "02_correlation")
    total = len(ctx.atlases) * len(ctx.sessions) * len(ctx.mice) * len(ctx.trials)
    with tqdm(total=total, desc="Correlating", leave=False) as pbar:
        try:
            for atlas in ctx.atlases:
                atlas_tag = "_".join(atlas)
                cache = CorrelationCache(ctx, atlas)
                for session in ctx.sessions:
                    target = Path(ctx.correlation_path) / f"{session}_correlation_matrices_{atlas_tag}_excluded.pkl"
                    df = pd.DataFrame(None, index=ctx.mice, columns=ctx.trials, dtype=object)
                    df[:] = None
                    for mouse in ctx.mice:
                        for trial in ctx.trials:
                            pbar.set_postfix(Session=session, Mouse=mouse, Trial=trial, Atlas=atlas_tag)
                            def correlate():
                                frame, status = cache.process(session, mouse, trial)
                                df.at[mouse, trial] = frame
                                return status
                            report.run(correlate, session=session, mouse=mouse, trial=trial,
                                atlas=atlas_tag, input_path=os.path.join(ctx.processed_path, session,
                                    mouse, trial, f"{mouse}_{trial}_FiltDRR.npy"))
                            pbar.update(1)
                    df.attrs["complete"] = all(isinstance(df.at[m, t], pd.DataFrame)
                                               for m in ctx.mice for t in ctx.trials)
                    df.attrs["batch_report"] = str(report.path)
                    save_pickle(target, df)
                    ctx.logger.info("FC table saved: %s; complete=%s", target, df.attrs["complete"])
        finally:
            report.finish()


if __name__ == "__main__":
    main()
