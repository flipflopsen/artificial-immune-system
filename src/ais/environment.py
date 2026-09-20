from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class CompartmentKind(StrEnum):
    BLOOD = "blood"
    VESSEL = "vessel"
    TISSUE = "tissue"
    ORGAN = "organ"
    LYMPH = "lymph"


@dataclass(slots=True)
class Compartment:
    name: str
    kind: CompartmentKind
    volume_ml: float

    health: float = 1.0
    inflammation: float = 0.0

    # These variables are principally relevant to viral dynamics.
    susceptible_cells: float = 0.0
    infected_cells: float = 0.0

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Compartment name must not be empty")
        if self.volume_ml <= 0.0:
            raise ValueError("Compartment volume must be positive")
        if self.susceptible_cells < 0.0 or self.infected_cells < 0.0:
            raise ValueError("Cell populations must be non-negative")


@dataclass(frozen=True, slots=True)
class Vessel:
    """
    Directed anatomical connection.

    transfer_rate_h:
        Fractional pathogen transport rate per hour.

    cell_migration_weight:
        Relative probability that mobile immune cells use this connection.
    """

    source: str
    target: str
    transfer_rate_h: float
    cell_migration_weight: float = 1.0

    def __post_init__(self) -> None:
        if self.transfer_rate_h < 0.0:
            raise ValueError("Transfer rate must be non-negative")
        if self.cell_migration_weight < 0.0:
            raise ValueError("Cell migration weight must be non-negative")


class Environment:
    """
    Directed graph representing blood vessels, lymphatic connections and organs.
    """

    def __init__(self) -> None:
        self.compartments: dict[str, Compartment] = {}
        self._outgoing: dict[str, list[Vessel]] = {}

    def add_compartment(self, compartment: Compartment) -> None:
        if compartment.name in self.compartments:
            raise ValueError(f"Duplicate compartment: {compartment.name!r}")

        self.compartments[compartment.name] = compartment
        self._outgoing[compartment.name] = []

    def connect(
        self,
        source: str,
        target: str,
        *,
        transfer_rate_h: float,
        cell_migration_weight: float = 1.0,
        bidirectional: bool = False,
    ) -> None:
        self._require_compartment(source)
        self._require_compartment(target)

        self._outgoing[source].append(
            Vessel(
                source=source,
                target=target,
                transfer_rate_h=transfer_rate_h,
                cell_migration_weight=cell_migration_weight,
            )
        )

        if bidirectional:
            self._outgoing[target].append(
                Vessel(
                    source=target,
                    target=source,
                    transfer_rate_h=transfer_rate_h,
                    cell_migration_weight=cell_migration_weight,
                )
            )

    def outgoing(self, source: str) -> tuple[Vessel, ...]:
        self._require_compartment(source)
        return tuple(self._outgoing[source])

    def get(self, name: str) -> Compartment:
        self._require_compartment(name)
        return self.compartments[name]

    def validate(self) -> None:
        if not self.compartments:
            raise ValueError("Environment has no compartments")

        for source, vessels in self._outgoing.items():
            for vessel in vessels:
                if vessel.source != source:
                    raise ValueError(f"Inconsistent vessel source: {vessel}")
                self._require_compartment(vessel.target)

    def _require_compartment(self, name: str) -> None:
        if name not in self.compartments:
            raise KeyError(f"Unknown anatomical compartment: {name!r}")
