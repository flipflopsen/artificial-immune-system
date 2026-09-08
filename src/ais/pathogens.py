from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from ais.models import PathogenKind, PathogenSpecies

if TYPE_CHECKING:
    from ais.simulation import Simulation


class PathogenModel(ABC):
    """
    Strategy interface for pathogen-specific local dynamics.

    Implementations mutate only the state belonging to the supplied species
    and compartment. Transport is handled independently by Simulation.
    """

    name: str
    supported_kind: PathogenKind

    @abstractmethod
    def advance(
        self,
        simulation: Simulation,
        species: PathogenSpecies,
        location: str,
        dt_h: float,
    ) -> None:
        raise NotImplementedError


class BacterialGrowthModel(PathogenModel):
    name = "bacterial-logistic"
    supported_kind = PathogenKind.BACTERIUM

    def advance(
        self,
        simulation: Simulation,
        species: PathogenSpecies,
        location: str,
        dt_h: float,
    ) -> None:
        population = simulation.pathogen_load(location, species.name)
        if population <= 0.0:
            return

        compartment = simulation.environment.get(location)

        growth_rate = species.parameter("growth_rate_h", 0.6)
        carrying_capacity_ml = species.parameter(
            "carrying_capacity_per_ml",
            1.0e8,
        )
        natural_clearance = species.parameter("natural_clearance_h", 0.01)

        carrying_capacity = carrying_capacity_ml * compartment.volume_ml
        logistic_factor = max(0.0, 1.0 - population / carrying_capacity)

        growth = growth_rate * population * logistic_factor * dt_h
        clearance = natural_clearance * population * dt_h

        simulation.set_pathogen_load(
            location,
            species.name,
            max(0.0, population + growth - clearance),
        )


class ViralReplicationModel(PathogenModel):
    name = "viral-target-cell"
    supported_kind = PathogenKind.VIRUS

    def advance(
        self,
        simulation: Simulation,
        species: PathogenSpecies,
        location: str,
        dt_h: float,
    ) -> None:
        virions = simulation.pathogen_load(location, species.name)
        compartment = simulation.environment.get(location)

        if virions <= 0.0 and compartment.infected_cells <= 0.0:
            return

        infectivity = species.parameter("infectivity_h", 0.02)
        half_saturation = species.parameter("target_half_saturation", 1.0e6)
        virion_production = species.parameter("virion_production_h", 20.0)
        infected_cell_death = species.parameter("infected_cell_death_h", 0.08)
        virion_clearance = species.parameter("virion_clearance_h", 0.15)

        target_fraction = compartment.susceptible_cells / (
            compartment.susceptible_cells + half_saturation
        )

        newly_infected = min(
            compartment.susceptible_cells,
            infectivity * virions * target_fraction * dt_h,
        )
        infected_deaths = min(
            compartment.infected_cells,
            infected_cell_death * compartment.infected_cells * dt_h,
        )

        produced_virions = virion_production * compartment.infected_cells * dt_h
        cleared_virions = min(
            virions,
            virion_clearance * virions * dt_h,
        )

        compartment.susceptible_cells -= newly_infected
        compartment.infected_cells += newly_infected - infected_deaths

        simulation.set_pathogen_load(
            location,
            species.name,
            max(0.0, virions + produced_virions - cleared_virions),
        )