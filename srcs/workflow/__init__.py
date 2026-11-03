"""Shared stage order; no scientific libraries are imported by the launcher."""
STAGES = ("preprocess", "fc", "edges", "geometry", "inference", "pairwise", "descriptive", "stability", "longitudinal", "combined", "statistics", "report")
MODULES = {name: "srcs.workflow." + name for name in STAGES}
MODULES.update(fc="srcs.workflow.correlation")
MODULES["seedfc"] = "plugins.seedfc"
