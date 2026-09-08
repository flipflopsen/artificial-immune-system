from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class PathogenKind(StrEnum):
    BACTERIUM = "bacterium"
    VIRUS = "virus"


@dataclass(frozen=True, slots=True)
class PathogenSpecies:
    """
    Immutable description of a pathogen strain.

    `model` selects a registered PathogenModel implementation. Model-specific
    parameters are held in `parameters`, which keeps the species abstraction
    independent of any particular growth equation.
    """

    name: str
    kind: PathogenKind
    model: str
    virulence: float = 1.0
    immune_evasion: float = 0.0
    parameters: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Pathogen species name must not be empty")
        if self.virulence < 0.0:
            raise ValueError("Virulence must be non-negative")
        if not 0.0 <= self.immune_evasion <= 1.0:
            raise ValueError("Immune evasion must be in [0, 1]")

    def parameter(self, name: str, default: float) -> float:
        return float(self.parameters.get(name, default))


@dataclass(frozen=True, slots=True)
class SimulationEvent:
    time_h: float
    event_type: str
    location: str | None
    details: dict[str, Any]