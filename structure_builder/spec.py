from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclass
class CompositionSpec:
    metal: Dict[str, float]
    nonmetal: Dict[str, float] = field(default_factory=dict)


@dataclass
class SQSOptions:
    iterations: int = 50000
    backend: str = "auto"


@dataclass
class CompositionalOptions:
    arrangement: str = "random"
    ordered_type: Optional[str] = None
    layer_axis: int = 2
    site_occupancies: Optional[List[Dict[str, Any]]] = None


@dataclass
class OutputSpec:
    formats: List[str] = field(default_factory=lambda: ["vasp", "cif"])
    directory: str = "data/structures/"


@dataclass
class SourceSpec:
    article: str = ""
    description: str = ""
    doi: str = ""


@dataclass
class StructureSpec:
    name: str
    base_structure: str
    type: str
    supercell_matrix: Any
    composition: CompositionSpec
    output: OutputSpec
    source: SourceSpec
    sqs_options: Optional[SQSOptions] = None
    compositional_options: Optional[CompositionalOptions] = None


def load_spec(path: str | Path) -> StructureSpec:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Spec file not found: {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    _validate(raw, path)

    comp_raw = raw["composition"]
    composition = CompositionSpec(
        metal=comp_raw["metal"],
        nonmetal=comp_raw.get("nonmetal", {}),
    )

    output_raw = raw.get("output", {})
    output = OutputSpec(
        formats=output_raw.get("formats", ["vasp", "cif"]),
        directory=output_raw.get("directory", "data/structures/"),
    )

    source_raw = raw.get("source", {})
    source = SourceSpec(
        article=source_raw.get("article", ""),
        description=source_raw.get("description", ""),
        doi=source_raw.get("doi", ""),
    )

    sqs_options = None
    if "sqs_options" in raw:
        o = raw["sqs_options"]
        sqs_options = SQSOptions(
            iterations=o.get("iterations", 50000),
            backend=o.get("backend", "auto"),
        )

    compositional_options = None
    if "compositional_options" in raw:
        o = raw["compositional_options"]
        compositional_options = CompositionalOptions(
            arrangement=o.get("arrangement", "random"),
            ordered_type=o.get("ordered_type"),
            layer_axis=o.get("layer_axis", 2),
            site_occupancies=o.get("site_occupancies"),
        )

    return StructureSpec(
        name=raw["name"],
        base_structure=raw["base_structure"],
        type=raw["type"],
        supercell_matrix=raw["supercell_matrix"],
        composition=composition,
        output=output,
        source=source,
        sqs_options=sqs_options,
        compositional_options=compositional_options,
    )


def _validate(raw: dict, path: Path) -> None:
    required = ["name", "base_structure", "type", "supercell_matrix", "composition"]
    missing = [f for f in required if f not in raw]
    if missing:
        raise ValueError(f"Spec {path} is missing required fields: {missing}")

    if raw["type"] not in ("sqs", "compositional"):
        raise ValueError(
            f"Spec {path}: type must be 'sqs' or 'compositional', got '{raw['type']}'"
        )

    comp = raw.get("composition", {})
    if "metal" not in comp:
        raise ValueError(f"Spec {path}: composition must have a 'metal' key")

    metal_sum = sum(comp["metal"].values())
    if abs(metal_sum - 1.0) > 0.01:
        raise ValueError(
            f"Spec {path}: metal fractions must sum to 1.0, got {metal_sum:.4f}"
        )

    if "nonmetal" in comp and comp["nonmetal"]:
        nonmetal_sum = sum(comp["nonmetal"].values())
        if abs(nonmetal_sum - 1.0) > 0.01:
            raise ValueError(
                f"Spec {path}: nonmetal fractions must sum to 1.0, got {nonmetal_sum:.4f}"
            )
