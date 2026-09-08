from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
from math import exp
from typing import Any, Iterable

import numpy as np

from ais.environment import Environment
from ais.immune.base import ImmuneCell
from ais.models import PathogenSpecies, SimulationEvent
from ais.pathogens import (
    BacterialGrowthModel,
    PathogenModel,
    ViralReplicationModel,
)


class Simulation:
    """
    Hybrid artificial immune-system simulation.

    Pathogens are represented as aggregate populations because an agent for
    every bacterium or virion would be computationally prohibitive. Immune
    cells remain individual agents because their heterogeneous state and
    behavior are biologically relevant.
    """

    def __init__(
        self,
        environment: Environment,
        *,
        seed: int | None = None,
        event_history_limit: int = 100_000,
    ) -> None:
        environment.validate()

        self.environment = environment
        self.rng = np.random.default_rng(seed)
        self.time_h = 0.0

        self.species: dict[str, PathogenSpecies] = {}
        self._pathogens: dict[tuple[str, str], float] = {}
        self.cells: dict[str, ImmuneCell] = {}

        self.signals: dict[tuple[str, str], float] = {}
        self.antigen_presentation: dict[str, float] = {}
        self.antibodies: dict[str, float] = {}

        self.events: list[SimulationEvent] = []
        self.event_history_limit = event_history_limit

        self.pathogen_models: dict[str, PathogenModel] = {}
        self.register_pathogen_model(BacterialGrowthModel())
        self.register_pathogen_model(ViralReplicationModel())

    def register_species(self, species: PathogenSpecies) -> None:
        if species.name in self.species:
            raise ValueError(f"Duplicate pathogen species: {species.name!r}")

        model = self.pathogen_models.get(species.model)
        if model is None:
            raise ValueError(f"Unknown pathogen model: {species.model!r}")
        if model.supported_kind is not species.kind:
            raise ValueError(
                f"Model {species.model!r} does not support {species.kind.value!r}"
            )

        self.species[species.name] = species

    def register_pathogen_model(self, model: PathogenModel) -> None:
        if model.name in self.pathogen_models:
            raise ValueError(f"Duplicate pathogen model: {model.name!r}")
        self.pathogen_models[model.name] = model

    def inject(
        self,
        species: PathogenSpecies | str,
        amount: float,
        *,
        location: str = "blood",
    ) -> None:
        if amount <= 0.0:
            raise ValueError("Injection amount must be positive")

        self.environment.get(location)

        if isinstance(species, PathogenSpecies):
            if species.name not in self.species:
                self.register_species(species)
            species_name = species.name
        else:
            species_name = species
            self._require_species(species_name)

        key = (location, species_name)
        self._pathogens[key] = self._pathogens.get(key, 0.0) + amount

        self.record_event(
            "pathogen-injection",
            location,
            species=species_name,
            amount=amount,
        )

    def add_cell(self, cell: ImmuneCell) -> None:
        self.environment.get(cell.location)

        if cell.identifier in self.cells:
            raise ValueError(f"Duplicate cell identifier: {cell.identifier!r}")

        self.cells[cell.identifier] = cell

    def add_cells(
        self,
        cell_type: type[ImmuneCell],
        count: int,
        *,
        location: str,
    ) -> None:
        if count < 0:
            raise ValueError("Cell count must be non-negative")

        for _ in range(count):
            self.add_cell(cell_type(location=location))

    def step(self, dt_h: float = 0.1) -> None:
        if dt_h <= 0.0:
            raise ValueError("Time step must be positive")

        self._transport_pathogens(dt_h)
        self._advance_pathogens(dt_h)
        self._advance_immune_cells(dt_h)
        self._apply_antibodies(dt_h)
        self._update_tissue_state(dt_h)
        self._decay_systemic_state(dt_h)

        self.time_h += dt_h

    def run(self, duration_h: float, *, dt_h: float = 0.1) -> None:
        if duration_h < 0.0:
            raise ValueError("Duration must be non-negative")

        full_steps = int(duration_h // dt_h)
        remainder = duration_h - full_steps * dt_h

        for _ in range(full_steps):
            self.step(dt_h)

        if remainder > 1.0e-12:
            self.step(remainder)

    def pathogen_load(self, location: str, species_name: str) -> float:
        return self._pathogens.get((location, species_name), 0.0)

    def set_pathogen_load(
        self,
        location: str,
        species_name: str,
        amount: float,
    ) -> None:
        self.environment.get(location)
        self._require_species(species_name)

        key = (location, species_name)

        if amount <= 1.0e-12:
            self._pathogens.pop(key, None)
        else:
            self._pathogens[key] = amount

    def remove_pathogens(
        self,
        location: str,
        species_name: str,
        amount: float,
    ) -> float:
        if amount < 0.0:
            raise ValueError("Removal amount must be non-negative")

        current = self.pathogen_load(location, species_name)
        removed = min(current, amount)
        self.set_pathogen_load(location, species_name, current - removed)
        return removed

    def total_pathogen_load(self, location: str | None = None) -> float:
        if location is None:
            return sum(self._pathogens.values())

        return sum(
            amount
            for (pathogen_location, _), amount in self._pathogens.items()
            if pathogen_location == location
        )

    def pathogens_at(
        self,
        location: str,
    ) -> Iterable[tuple[PathogenSpecies, float]]:
        snapshot = [
            (self.species[species_name], amount)
            for (pathogen_location, species_name), amount in self._pathogens.items()
            if pathogen_location == location and amount > 0.0
        ]
        return tuple(snapshot)

    def emit_signal(self, location: str, signal_name: str, amount: float) -> None:
        if amount <= 0.0:
            return

        key = (location, signal_name)
        self.signals[key] = self.signals.get(key, 0.0) + amount

    def signal(self, location: str, signal_name: str) -> float:
        return self.signals.get((location, signal_name), 0.0)

    def record_event(
        self,
        event_type: str,
        location: str | None,
        **details: Any,
    ) -> None:
        self.events.append(
            SimulationEvent(
                time_h=self.time_h,
                event_type=event_type,
                location=location,
                details=details,
            )
        )

        excess = len(self.events) - self.event_history_limit
        if excess > 0:
            del self.events[:excess]

    def snapshot(self) -> dict[str, Any]:
        return {
            "time_h": self.time_h,
            "total_pathogen_load": self.total_pathogen_load(),
            "pathogens": {
                f"{location}:{species}": amount
                for (location, species), amount in sorted(self._pathogens.items())
            },
            "compartments": {
                name: asdict(compartment)
                for name, compartment in self.environment.compartments.items()
            },
            "immune_cells": {
                cell_type: sum(
                    1 for cell in self.cells.values() if cell.cell_type == cell_type
                )
                for cell_type in sorted({cell.cell_type for cell in self.cells.values()})
            },
            "antigen_presentation": dict(self.antigen_presentation),
            "antibodies": dict(self.antibodies),
        }

    def _transport_pathogens(self, dt_h: float) -> None:
        """
        Apply transport simultaneously.

        Using a delta buffer prevents graph iteration order from influencing
        the biological result within one simulation step.
        """

        original = dict(self._pathogens)
        deltas: defaultdict[tuple[str, str], float] = defaultdict(float)

        for (source, species_name), amount in original.items():
            vessels = self.environment.outgoing(source)
            if not vessels or amount <= 0.0:
                continue

            requested: list[tuple[str, float]] = []

            for vessel in vessels:
                fraction = 1.0 - exp(-vessel.transfer_rate_h * dt_h)
                requested.append((vessel.target, amount * fraction))

            requested_total = sum(transfer for _, transfer in requested)
            scaling = min(1.0, amount / requested_total) if requested_total > 0.0 else 1.0

            for target, transfer in requested:
                actual_transfer = transfer * scaling
                deltas[(source, species_name)] -= actual_transfer
                deltas[(target, species_name)] += actual_transfer

        for key, delta in deltas.items():
            self.set_pathogen_load(
                key[0],
                key[1],
                max(0.0, self._pathogens.get(key, 0.0) + delta),
            )

    def _advance_pathogens(self, dt_h: float) -> None:
        active_pairs = tuple(self._pathogens.keys())

        for location, species_name in active_pairs:
            species = self.species[species_name]
            model = self.pathogen_models[species.model]
            model.advance(self, species, location, dt_h)

    def _advance_immune_cells(self, dt_h: float) -> None:
        dead_cells: list[str] = []

        for identifier, cell in tuple(self.cells.items()):
            cell.act(self, dt_h)
            cell.age_h += dt_h

            if cell.age_h >= cell.lifespan_h:
                dead_cells.append(identifier)
                continue

            self._migrate_cell(cell, dt_h)

        for identifier in dead_cells:
            cell = self.cells.pop(identifier)
            self.record_event(
                "immune-cell-death",
                cell.location,
                cell=cell.cell_type,
                identifier=identifier,
            )

    def _migrate_cell(self, cell: ImmuneCell, dt_h: float) -> None:
        migration_probability = 1.0 - exp(-cell.mobility_h * dt_h)
        if self.rng.random() >= migration_probability:
            return

        vessels = self.environment.outgoing(cell.location)
        if not vessels:
            return

        weights = np.asarray(
            [
                vessel.cell_migration_weight
                * (1.0 + 5.0 * self.environment.get(vessel.target).inflammation)
                for vessel in vessels
            ],
            dtype=float,
        )

        weight_sum = float(weights.sum())
        if weight_sum <= 0.0:
            return

        probabilities = weights / weight_sum
        selected = int(self.rng.choice(len(vessels), p=probabilities))
        cell.location = vessels[selected].target

    def _apply_antibodies(self, dt_h: float) -> None:
        for species_name, antibody_level in tuple(self.antibodies.items()):
            if antibody_level <= 0.0:
                continue

            species = self.species.get(species_name)
            if species is None:
                continue

            efficacy = (1.0 - species.immune_evasion) * antibody_level
            clearance_fraction = 1.0 - exp(-0.05 * efficacy * dt_h)

            for location in self.environment.compartments:
                load = self.pathogen_load(location, species_name)
                self.remove_pathogens(
                    location,
                    species_name,
                    load * clearance_fraction,
                )

    def _update_tissue_state(self, dt_h: float) -> None:
        for location, compartment in self.environment.compartments.items():
            weighted_burden = 0.0

            for species, population in self.pathogens_at(location):
                weighted_burden += species.virulence * population

            target_inflammation = weighted_burden / (weighted_burden + 1.0e5)
            inflammatory_signal = self.signal(location, "pro-inflammatory")
            target_inflammation = max(
                target_inflammation,
                inflammatory_signal / (inflammatory_signal + 100.0),
            )

            compartment.inflammation += (
                target_inflammation - compartment.inflammation
            ) * min(1.0, 0.5 * dt_h)
            compartment.inflammation = min(
                1.0,
                max(0.0, compartment.inflammation),
            )

            damage_rate = 0.002 * compartment.inflammation
            recovery_rate = 0.0005 * (1.0 - compartment.inflammation)

            compartment.health += (
                recovery_rate * (1.0 - compartment.health) - damage_rate
            ) * dt_h
            compartment.health = min(1.0, max(0.0, compartment.health))

    def _decay_systemic_state(self, dt_h: float) -> None:
        self._decay_mapping(self.signals, rate_h=0.7, dt_h=dt_h)
        self._decay_mapping(
            self.antigen_presentation,
            rate_h=0.02,
            dt_h=dt_h,
        )
        self._decay_mapping(self.antibodies, rate_h=0.005, dt_h=dt_h)

    @staticmethod
    def _decay_mapping(
        values: dict[Any, float],
        *,
        rate_h: float,
        dt_h: float,
    ) -> None:
        factor = exp(-rate_h * dt_h)

        for key in tuple(values):
            values[key] *= factor
            if values[key] <= 1.0e-12:
                del values[key]

    def _require_species(self, species_name: str) -> None:
        if species_name not in self.species:
            raise KeyError(f"Unknown pathogen species: {species_name!r}")