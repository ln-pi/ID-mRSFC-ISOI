"""Build the ordered workflow from core stages and enabled plugins."""
from dataclasses import dataclass, field

STAGES = ("preprocess", "fc", "edges", "pairwise", "descriptive", "stability",
          "longitudinal", "combined", "statistics", "report")
REPORT_SECTIONS = ("pairwise", "descriptive", "stability", "longitudinal", "combined", "statistics")


@dataclass(frozen=True)
class Stage:
    module: str
    parents: tuple = ()
    products: dict = field(default_factory=dict)
    validate: object = None
    source: bool = False


CORE = {
    "preprocess": Stage("srcs.workflow.preprocess", products={"preprocessed_recordings": ("*",)}),
    "fc": Stage("srcs.workflow.fc", products={"fc_tables": ("*.pkl",), "batch_report": ("*.json",)}),
    "edges": Stage("srcs.workflow.edges", ("fc",), {"edge_store": ("*",)}),
    "pairwise": Stage("srcs.workflow.pairwise", ("edges",), {"pairwise_tensors": ("*.npz",), "pairwise_selection": ("*.json","*.csv")}),
    "descriptive": Stage("srcs.workflow.descriptive", ("edges",), {"descriptive": ("*",)}),
    "stability": Stage("srcs.workflow.stability", ("edges","descriptive"), {"stability": ("*",)}),
    "longitudinal": Stage("srcs.workflow.longitudinal", ("edges","descriptive"), {"longitudinal": ("*",)}),
    "combined": Stage("srcs.workflow.combined", ("edges",), {"combined_identification": ("*",), "identifiability_aggregate": ("@aggregate",)}),
    "statistics": Stage("srcs.workflow.statistics", ("combined",), {"statistics_tables": ("*",)}),
    "report": Stage("srcs.workflow.report", products={"unified_report": ("*",)}),
}


class Workflow:
    def __init__(self, config, root, *, resolve_source=True):
        from srcs.workflow.plugins import enabled_plugins
        from srcs.workflow.sources import active_source

        self.source = active_source(config, root) if resolve_source else None
        self.functions = enabled_plugins(config, "function", root)
        self.stages = dict(CORE)
        self.reports = {}
        self.option_sections = {}
        order = list(STAGES)
        if self.source:
            name = self.source.spec["id"]
            if name in self.stages:
                raise ValueError(f"Duplicate stage {name}")
            self.stages[name] = Stage(self.source.module, products=CORE["fc"].products,
                                      validate=self.source.spec["validate"], source=True)
            self.stages["edges"] = Stage(CORE["edges"].module, (name,), CORE["edges"].products)
            order = [name, *order[order.index("edges"):]]
        for plugin in self.functions:
            for name, spec in plugin.spec.get("stages", {}).items():
                if name in self.stages:
                    raise ValueError(f"Duplicate stage {name}")
                after = spec["after"]
                parents = tuple(spec.get("requires", ()))
                if after not in order or any(parent not in order[:order.index(after)+1] for parent in parents):
                    raise ValueError(f"Unavailable or unordered dependencies for {name}")
                self.stages[name] = Stage(plugin.stage_module(spec), parents, spec["produces"],
                                          spec.get("validate"))
                order.insert(order.index(after) + 1, name)
            for name in plugin.spec.get("option_sections", ()):
                if name in self.option_sections:
                    raise ValueError(f"Duplicate plugin option section {name}")
                self.option_sections[name] = plugin.spec["id"]
            for name, spec in plugin.spec.get("reports", {}).items():
                if name in self.reports or name in REPORT_SECTIONS or name == "all":
                    raise ValueError(f"Duplicate report section {name}")
                self.reports[name] = (plugin, spec)
        self.order = tuple(order)

    def validate_plan(self, stages):
        from importlib.util import find_spec

        seen = set()
        for name in stages:
            stage = self.stages[name]
            if find_spec(stage.module) is None:
                raise ValueError(f"Missing implementation: {name}")
            collision = seen.intersection(stage.products)
            if collision:
                raise ValueError(f"Multiple producers: {sorted(collision)}")
            seen.update(stage.products)
