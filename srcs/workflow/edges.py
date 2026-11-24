from srcs.core.runtime import numerical_stage
from srcs.utils.edge_store import run as build_edges

OPTION_KEYS = ()


@numerical_stage("edges", OPTION_KEYS)
def run(ctx):
    build_edges(ctx)
