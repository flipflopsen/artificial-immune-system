from __future__ import annotations

from typing import ClassVar

from ais.immune.base import ImmuneCell
from ais.models import PathogenKind


def saturation(value: float, half_saturation: float) -> float:
    if value <= 0.0:
        return 0.0
    return value / (value + half_saturation)


class Macrophage(ImmuneCell):
    cell_type: ClassVar[str] = "macrophage"
    lifespan_h: ClassVar[float] = 24.0 * 60.0
    mobility_h: ClassVar[float] = 0.02

    def act(self, simulation: Simulation, dt_h: float) -> None:
        local_burden = simulation.total_pathogen_load(self.location)
        self.activate(saturation(local_burden, 1.0e4), dt_h)

        for species, population in simulation.pathogens_at(self.location):
            if species.kind is not PathogenKind.BACTERIUM:
                continue

            recognition = 1.0 - species.immune_evasion
            capacity = 30.0 * (0.25 + self.activation) * recognition * dt_h
            removed = simulation.remove_pathogens(
                self.location,
                species.name,
                min(population, capacity),
            )

            if removed > 0.0:
                simulation.emit_signal(
                    self.location,
                    "pro-inflammatory",
                    0.01 * removed,
                )
                simulation.record_event(
                    "phagocytosis",
                    self.location,
                    cell=self.cell_type,
                    species=species.name,
                    amount=removed,
                )


class Neutrophil(ImmuneCell):
    cell_type: ClassVar[str] = "neutrophil"
    lifespan_h: ClassVar[float] = 72.0
    mobility_h: ClassVar[float] = 0.35

    def act(self, simulation: Simulation, dt_h: float) -> None:
        local_burden = simulation.total_pathogen_load(self.location)
        inflammatory_signal = simulation.signal(self.location, "pro-inflammatory")

        stimulus = max(
            saturation(local_burden, 5.0e3),
            saturation(inflammatory_signal, 10.0),
        )
        self.activate(stimulus, dt_h)

        for species, population in simulation.pathogens_at(self.location):
            if species.kind is not PathogenKind.BACTERIUM:
                continue

            recognition = 1.0 - species.immune_evasion
            capacity = 75.0 * self.activation * recognition * dt_h
            removed = simulation.remove_pathogens(
                self.location,
                species.name,
                min(population, capacity),
            )

            if removed > 0.0:
                simulation.emit_signal(
                    self.location,
                    "pro-inflammatory",
                    0.02 * removed,
                )


class DendriticCell(ImmuneCell):
    cell_type: ClassVar[str] = "dendritic-cell"
    lifespan_h: ClassVar[float] = 24.0 * 14.0
    mobility_h: ClassVar[float] = 0.08

    def act(self, simulation: Simulation, dt_h: float) -> None:
        burden = simulation.total_pathogen_load(self.location)
        self.activate(saturation(burden, 1.0e3), dt_h)

        if self.activation <= 0.0:
            return

        for species, population in simulation.pathogens_at(self.location):
            presentation = (
                self.activation
                * saturation(population, 1.0e3)
                * (1.0 - species.immune_evasion)
                * dt_h
            )
            simulation.antigen_presentation[species.name] = (
                simulation.antigen_presentation.get(species.name, 0.0) + presentation
            )


class NaturalKillerCell(ImmuneCell):
    cell_type: ClassVar[str] = "natural-killer-cell"
    lifespan_h: ClassVar[float] = 24.0 * 14.0
    mobility_h: ClassVar[float] = 0.15

    def act(self, simulation: Simulation, dt_h: float) -> None:
        compartment = simulation.environment.get(self.location)
        stimulus = saturation(compartment.infected_cells, 100.0)
        self.activate(stimulus, dt_h)

        killed = min(
            compartment.infected_cells,
            8.0 * self.activation * dt_h,
        )
        compartment.infected_cells -= killed

        if killed > 0.0:
            simulation.record_event(
                "infected-cell-killing",
                self.location,
                cell=self.cell_type,
                amount=killed,
            )


class CytotoxicTCell(ImmuneCell):
    cell_type: ClassVar[str] = "cytotoxic-t-cell"
    lifespan_h: ClassVar[float] = 24.0 * 365.0
    mobility_h: ClassVar[float] = 0.18

    def act(self, simulation: Simulation, dt_h: float) -> None:
        local_viral_memory = 0.0

        for species, _ in simulation.pathogens_at(self.location):
            if species.kind is PathogenKind.VIRUS:
                local_viral_memory = max(
                    local_viral_memory,
                    simulation.antigen_presentation.get(species.name, 0.0),
                )

        compartment = simulation.environment.get(self.location)
        stimulus = saturation(local_viral_memory, 5.0)
        self.activate(stimulus, dt_h)

        killed = min(
            compartment.infected_cells,
            20.0 * self.activation * dt_h,
        )
        compartment.infected_cells -= killed


class BCell(ImmuneCell):
    cell_type: ClassVar[str] = "b-cell"
    lifespan_h: ClassVar[float] = 24.0 * 365.0
    mobility_h: ClassVar[float] = 0.08

    def act(self, simulation: Simulation, dt_h: float) -> None:
        strongest_presentation = max(
            simulation.antigen_presentation.values(),
            default=0.0,
        )
        self.activate(saturation(strongest_presentation, 5.0), dt_h)

        for species_name, presentation in simulation.antigen_presentation.items():
            production = 0.25 * self.activation * saturation(presentation, 5.0) * dt_h
            simulation.antibodies[species_name] = (
                simulation.antibodies.get(species_name, 0.0) + production
            )


CELL_TYPES: dict[str, type[ImmuneCell]] = {
    Macrophage.cell_type: Macrophage,
    Neutrophil.cell_type: Neutrophil,
    DendriticCell.cell_type: DendriticCell,
    NaturalKillerCell.cell_type: NaturalKillerCell,
    CytotoxicTCell.cell_type: CytotoxicTCell,
    BCell.cell_type: BCell,
}


# Imported only for static type checking semantics; kept at the bottom to avoid
# a runtime import cycle between simulation.py and cells.py.
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ais.simulation import Simulation
