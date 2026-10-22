"""Shared stage order; no scientific libraries are imported by the launcher."""
STAGES = ("preprocess", "fc", "edges", "pairwise", "descriptive", "combined", "statistics", "report")
MODULES = {name: "srcs.workflow." + name for name in STAGES}
MODULES.update(fc="srcs.workflow.correlation", edges="srcs.utils.edge_store")
