from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class ProductRef:
    name: str
    schema: str
    path: Path
    sha256: str


@dataclass(frozen=True)
class StageResult:
    products: tuple[ProductRef, ...] = ()


@dataclass(frozen=True)
class StageContext:
    inputs: Mapping[str, ProductRef]
    options: Mapping
    option_keys: tuple[str, ...]
    output_root: Path
    data_tag: str

    def section(self, name):
        if name not in self.option_keys:
            raise KeyError(f"Undeclared configuration section: {name}")
        return self.options.get(name, {})
