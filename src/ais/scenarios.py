from __future__ import annotations

from ais.environment import Compartment, CompartmentKind, Environment
from ais.immune import (
    BCell,
    CytotoxicTCell,
    DendriticCell,
    Macrophage,
    NaturalKillerCell,
    Neutrophil,
)
from ais.models import PathogenKind, PathogenSpecies
from ais.simulation import Simulation


def create_human_environment() -> Environment:
    environment = Environment()

    environment.add_compartment(
        Compartment(
            name="blood",
            kind=CompartmentKind.BLOOD,
            volume_ml=5_000.0,
            susceptible_cells=2.0e8,
        )
    )
    environment.add_compartment(
        Compartment(
            name="lungs",
            kind=CompartmentKind.ORGAN,
            volume_ml=1_000.0,
            susceptible_cells=5.0e8,
        )
    )
    environment.add_compartment(
        Compartment(
            name="liver",
            kind=CompartmentKind.ORGAN,
            volume_ml=1_500.0,
            susceptible_cells=3.0e8,
        )
    )
    environment.add_compartment(
        Compartment(
            name="spleen",
            kind=CompartmentKind.ORGAN,
            volume_ml=200.0,
            susceptible_cells=5.0e7,
        )
    )
    environment.add_compartment(
        Compartment(
            name="lymph",
            kind=CompartmentKind.LYMPH,
            volume_ml=2_000.0,
            susceptible_cells=1.0e8,
        )
    )
    environment.add_compartment(
        Compartment(
            name="peripheral-tissue",
            kind=CompartmentKind.TISSUE,
            volume_ml=8_000.0,
            susceptible_cells=1.0e9,
        )
    )

    environment.connect(
        "blood",
        "lungs",
        transfer_rate_h=0.12,
        cell_migration_weight=1.3,
        bidirectional=True,
    )
    environment.connect(
        "blood",
        "liver",
        transfer_rate_h=0.08,
        cell_migration_weight=1.0,
        bidirectional=True,
    )
    environment.connect(
        "blood",
        "spleen",
        transfer_rate_h=0.06,
        cell_migration_weight=1.5,
        bidirectional=True,
    )
    environment.connect(
        "blood",
        "peripheral-tissue",
        transfer_rate_h=0.10,
        cell_migration_weight=1.0,
        bidirectional=True,
    )
    environment.connect(
        "peripheral-tissue",
        "lymph",
        transfer_rate_h=0.03,
        cell_migration_weight=1.2,
    )
    environment.connect(
        "lymph",
        "blood",
        transfer_rate_h=0.04,
        cell_migration_weight=1.0,
    )
    environment.connect(
        "spleen",
        "lymph",
        transfer_rate_h=0.04,
        cell_migration_weight=1.2,
        bidirectional=True,
    )

    return environment


def escherichia_coli() -> PathogenSpecies:
    return PathogenSpecies(
        name="E. coli",
        kind=PathogenKind.BACTERIUM,
        model="bacterial-logistic",
        virulence=0.7,
        immune_evasion=0.1,
        parameters={
            "growth_rate_h": 0.75,
            "carrying_capacity_per_ml": 1.0e8,
            "natural_clearance_h": 0.01,
        },
    )


def influenza_a() -> PathogenSpecies:
    return PathogenSpecies(
        name="Influenza A",
        kind=PathogenKind.VIRUS,
        model="viral-target-cell",
        virulence=0.8,
        immune_evasion=0.2,
        parameters={
            "infectivity_h": 0.03,
            "target_half_saturation": 1.0e6,
            "virion_production_h": 30.0,
            "infected_cell_death_h": 0.08,
            "virion_clearance_h": 0.18,
        },
    )


def create_baseline_simulation(seed: int | None = 42) -> Simulation:
    simulation = Simulation(create_human_environment(), seed=seed)

    simulation.register_species(escherichia_coli())
    simulation.register_species(influenza_a())

    # These counts are simulation agents, not literal physiological counts.
    # Each agent may be interpreted as a representative cellular cohort.
    simulation.add_cells(Macrophage, 100, location="peripheral-tissue")
    simulation.add_cells(Macrophage, 40, location="lungs")
    simulation.add_cells(Neutrophil, 250, location="blood")
    simulation.add_cells(DendriticCell, 30, location="peripheral-tissue")
    simulation.add_cells(DendriticCell, 20, location="lungs")
    simulation.add_cells(NaturalKillerCell, 80, location="blood")
    simulation.add_cells(CytotoxicTCell, 60, location="lymph")
    simulation.add_cells(BCell, 60, location="spleen")

    return simulation
