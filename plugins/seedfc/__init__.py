"""Read organized SeedFC records and publish labeled native FC week tables."""
from pathlib import Path
import numpy as np
import pandas as pd

from srcs.analysis.connectivity import valid_fc
from srcs.core.artifacts import sha256
from srcs.core.identity import digest
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.results import save_json, save_pickle


def locations(ctx, options):
    root = (PROJECT_ROOT / options["source"] / options["branch"] / options["series"]).resolve()
    for week in ctx.weeks:
        for mouse in ctx.mice:
            for run in ctx.runs:
                yield week, mouse, run, root/week/run/mouse/"FC_matrices.npz"


def input_files(ctx, options):
    return {str(p): sha256(p) if p.is_file() else None for *_,p in locations(ctx,options)}


def run(ctx, dry_run=False):
    options = ctx.config["plugins"]["seedfc"]["options"]
    partitions = options["partitions"]
    if list(partitions) != ["_".join(a) for a in ctx.atlases]:
        raise ValueError("SeedFC partitions must match configured atlas order")
    tables = {(w,part): pd.DataFrame(None,index=ctx.mice,columns=ctx.runs,dtype=object)
              for w in ctx.weeks for part in partitions}
    records, basis = [], None
    for w,m,r,path in locations(ctx,options):
        matrix = None
        if path.exists():
            with np.load(path,allow_pickle=False) as data:
                labels = data["seed_names"].astype(str).tolist()
                if len(set(labels)) != len(labels) or (basis is not None and basis != labels):
                    raise ValueError(f"Inconsistent SeedFC labels: {path}")
                basis = labels
                matrix = pd.DataFrame(data["r"],index=labels,columns=labels)
        for part,rules in partitions.items():
            status, reason = "missing", "missing source file"
            if matrix is not None:
                excluded = rules["exclude"]
                if not isinstance(excluded,list) or not set(excluded).issubset(basis):
                    raise ValueError(f"Unknown seed exclusion in {part}")
                selected = [s for s in basis if s not in excluded]
                if len(selected)<3:
                    raise ValueError(f"Fewer than three retained seeds in {part}")
                frame = matrix.loc[selected,selected]
                if valid_fc(frame,selected):
                    tables[w,part].at[m,r] = frame
                    status,reason = "success",""
                elif np.isfinite(frame).all().all():
                    raise ValueError(f"Invalid finite FC: {path}")
                else:
                    status,reason = "missing","nonfinite FC"
            records.append(dict(week=w,mouse=m,run=r,atlas=part,status=status,reason=reason,input_path=str(path)))
    if basis is None:
        raise FileNotFoundError("No configured SeedFC records found")
    if dry_run:
        return
    products = []
    for (week,part),table in tables.items():
        path = Path(ctx.correlation_path)/f"{week}_correlation_matrices_{part}_excluded.pkl"
        save_pickle(path,table)
        products.append({"name":path.name,"sha256":sha256(path)})
    save_json(Path(ctx.output_dir("BatchReports"))/f"02_correlation_seedfc_{options['branch']}.json",
              {"records":records})
    save_json(Path(ctx.output_dir("Metadata"))/"fc_import_provenance.json",
              {"config_tag":ctx.analysis_tag,"identity_digest":digest(ctx.data_identity),
               "source":"seedfc","tables":products})
