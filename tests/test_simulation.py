from __future__ import annotations

import pytest

from ais.environment import Compartment, CompartmentKind, Environment
from ais.immune import Macrophage
from ais.models import PathogenKind, PathogenSpecies
from ais.simulation import Simulation


@pytest.fixture
def bacterium() -> PathogenSpecies:
    return PathogenSpecies(
        name="test-bacterium",
        kind=PathogenKind.BACTERIUM,
        model="bacterial-logistic",
        parameters={
            "growth_rate_h": 0.5,
            "carrying_capacity_per_ml": 1.0e8,
            "natural_clearance_h": 0.0,
        },
    )


def simple_environment() -> Environment:
    environment = Environment()
    environment.add_compartment(
        Compartment(
            name="blood",
            kind=CompartmentKind.BLOOD,
            volume_ml=1_000.0,
        )
    )
    environment.add_compartment(
        Compartment(
            name="organ",
            kind=CompartmentKind.ORGAN,
            volume_ml=100.0,
        )
    )
    environment.connect(
        "blood",
        "organ",
        transfer_rate_h=0.2,
        bidirectional=True,
    )
    return environment


def test_injection_registers_species_and_population(
    bacterium: PathogenSpecies,
) -> None:
    simulation = Simulation(simple_environment(), seed=1)

    simulation.inject(bacterium, 1_000.0, location="blood")

    assert simulation.pathogen_load("blood", bacterium.name) == pytest.approx(
        1_000.0
    )


def test_bacteria_grow_without_immune_response(
    bacterium: PathogenSpecies,
) -> None:
    simulation = Simulation(simple_environment(), seed=1)
    simulation.register_species(bacterium)
    simulation.inject(bacterium.name, 1_000.0, location="blood")

    before = simulation.total_pathogen_load()
    simulation.run(1.0, dt_h=0.05)
    after = simulation.total_pathogen_load()

    assert after > before


def test_macrophages_reduce_bacterial_load(
    bacterium: PathogenSpecies,
) -> None:
    control = Simulation(simple_environment(), seed=1)
    control.register_species(bacterium)
    control.inject(bacterium.name, 10_000.0, location="blood")

    treated = Simulation(simple_environment(), seed=1)
    treated.register_species(bacterium)
    treated.inject(bacterium.name, 10_000.0, location="blood")
    treated.add_cells(Macrophage, 500, location="blood")

    control.run(2.0, dt_h=0.05)
    treated.run(2.0, dt_h=0.05)

    assert treated.total_pathogen_load() < control.total_pathogen_load()


def test_pathogens_are_transported_between_compartments(
    bacterium: PathogenSpecies,
) -> None:
    simulation = Simulation(simple_environment(), seed=1)
    simulation.register_species(bacterium)
    simulation.inject(bacterium.name, 1_000.0, location="blood")

    simulation.step(0.5)

    assert simulation.pathogen_load("organ", bacterium.name) > 0.0


def test_invalid_injection_is_rejected(
    bacterium: PathogenSpecies,
) -> None:
    simulation = Simulation(simple_environment())

    with pytest.raises(ValueError):
        simulation.inject(bacterium, -1.0)

    with pytest.raises(KeyError):
        simulation.inject(bacterium, 10.0, location="nonexistent")


def test_simulation_is_deterministic_for_fixed_seed(
    bacterium: PathogenSpecies,
) -> None:
    first = Simulation(simple_environment(), seed=17)
    second = Simulation(simple_environment(), seed=17)

    for simulation in (first, second):
        simulation.register_species(bacterium)
        simulation.inject(bacterium.name, 10_000.0)
        simulation.add_cells(Macrophage, 20, location="blood")
        simulation.run(5.0, dt_h=0.1)

    assert first.snapshot() == second.snapshot()