# Artificial Immune System: Examples and Extension Guide

This document describes how to run simulations and extend the environment,
immune-cell system and pathogen models.

Commands assume that the project environment has already been synchronized:

```bash
uv sync --extra dev
```

Use `uv run` to execute commands without manually activating `.venv`.

## Run a bacterial infection

Inject `100,000` bacterial units into the bloodstream and simulate 72 hours
using a 0.1-hour integration step:

```bash
uv run ais-simulate \
    --pathogen bacteria \
    --amount 100000 \
    --location blood \
    --duration 72 \
    --dt 0.1
```

If the virtual environment has been activated manually, the command can also be
executed directly:

```bash
ais-simulate \
    --pathogen bacteria \
    --amount 100000 \
    --location blood \
    --duration 72 \
    --dt 0.1
```

## Run a viral infection

Inject `50,000` viral units into the lungs and simulate 168 hours:

```bash
uv run ais-simulate \
    --pathogen virus \
    --amount 50000 \
    --location lungs \
    --duration 168 \
    --dt 0.1
```

The baseline command-line interface maps:

- `bacteria` to `E. coli`;
- `virus` to `Influenza A`.

More specific species selection can be performed through the programmatic API
or by extending the command-line interface.

## Programmatic API

```python
from ais.scenarios import create_baseline_simulation


simulation = create_baseline_simulation(seed=42)

simulation.inject(
    "E. coli",
    100_000,
    location="blood",
)

simulation.run(
    72.0,
    dt_h=0.1,
)

print(simulation.snapshot())
```

A fixed seed makes stochastic operations reproducible for the same source code,
dependency versions and initial conditions.

### Advance the simulation incrementally

Use `step()` when measurements or interventions must be applied during a run:

```python
from ais.scenarios import create_baseline_simulation


simulation = create_baseline_simulation(seed=42)
simulation.inject("Influenza A", 50_000, location="lungs")

for _ in range(24 * 10):
    simulation.step(dt_h=0.1)

    if simulation.time_h >= 12.0:
        lung_load = simulation.total_pathogen_load("lungs")
        print(f"{simulation.time_h:.1f} h: lung load = {lung_load:.3f}")
```

## Extension points

The principal extension interfaces are:

- `ImmuneCell` for immune-cell phenotypes;
- `PathogenModel` for local pathogen dynamics;
- `PathogenSpecies` for strain definitions;
- `Environment`, `Compartment` and `Vessel` for anatomical topology;
- scenario factory functions for reproducible experimental configurations.

## Add an immune-cell phenotype

Subclass `ImmuneCell` and implement `act()`.

The following regulatory T-cell model suppresses a local pro-inflammatory
signal:

```python
from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from ais.immune.base import ImmuneCell

if TYPE_CHECKING:
    from ais.simulation import Simulation


class RegulatoryTCell(ImmuneCell):
    cell_type: ClassVar[str] = "regulatory-t-cell"
    lifespan_h: ClassVar[float] = 24.0 * 365.0
    mobility_h: ClassVar[float] = 0.1

    def act(self, simulation: Simulation, dt_h: float) -> None:
        signal_key = (self.location, "pro-inflammatory")
        local_signal = simulation.signal(
            self.location,
            "pro-inflammatory",
        )

        suppressed = min(local_signal, 2.0 * dt_h)
        remaining_signal = local_signal - suppressed

        if remaining_signal <= 1.0e-12:
            simulation.signals.pop(signal_key, None)
        else:
            simulation.signals[signal_key] = remaining_signal
```

Add instances to a simulation:

```python
from ais.scenarios import create_baseline_simulation


simulation = create_baseline_simulation(seed=42)
simulation.add_cells(
    RegulatoryTCell,
    count=25,
    location="lymph",
)
```

Optionally register the class in `CELL_TYPES` when string-based cell
construction or serialization is required:

```python
from ais.immune.cells import CELL_TYPES


CELL_TYPES[RegulatoryTCell.cell_type] = RegulatoryTCell
```

For a larger extension, place the implementation in a separate module, for
example:

```text
src/ais/immune/regulatory.py
```

Then export it from `src/ais/immune/__init__.py`.

## Add a pathogen model

Subclass `PathogenModel`, register an instance with the simulation and reference
its name from a `PathogenSpecies`.

```python
from __future__ import annotations

from typing import TYPE_CHECKING

from ais.models import PathogenKind
from ais.pathogens import PathogenModel

if TYPE_CHECKING:
    from ais.models import PathogenSpecies
    from ais.simulation import Simulation


class IntracellularBacterialModel(PathogenModel):
    name = "intracellular-bacterial"
    supported_kind = PathogenKind.BACTERIUM

    def advance(
        self,
        simulation: Simulation,
        species: PathogenSpecies,
        location: str,
        dt_h: float,
    ) -> None:
        current = simulation.pathogen_load(
            location,
            species.name,
        )

        growth_rate_h = species.parameter(
            "growth_rate_h",
            0.2,
        )
        natural_clearance_h = species.parameter(
            "natural_clearance_h",
            0.01,
        )

        net_rate_h = growth_rate_h - natural_clearance_h
        updated = max(
            0.0,
            current * (1.0 + net_rate_h * dt_h),
        )

        simulation.set_pathogen_load(
            location,
            species.name,
            updated,
        )
```

Register the model before registering a species that uses it:

```python
from ais.models import PathogenKind, PathogenSpecies
from ais.scenarios import create_human_environment
from ais.simulation import Simulation


simulation = Simulation(
    create_human_environment(),
    seed=42,
)

simulation.register_pathogen_model(
    IntracellularBacterialModel()
)

simulation.register_species(
    PathogenSpecies(
        name="Intracellular bacterium",
        kind=PathogenKind.BACTERIUM,
        model="intracellular-bacterial",
        virulence=0.8,
        immune_evasion=0.4,
        parameters={
            "growth_rate_h": 0.25,
            "natural_clearance_h": 0.015,
        },
    )
)

simulation.inject(
    "Intracellular bacterium",
    10_000,
    location="blood",
)

simulation.run(48.0, dt_h=0.1)
```

The model must be registered before its corresponding species because
`register_species()` validates that the selected model exists and supports the
species' pathogen kind.

### Numerical consideration

The example uses explicit Euler integration:

```text
N(t + Δt) = N(t) × (1 + r × Δt)
```

For stiff, strongly nonlinear or multiscale pathogen dynamics, implement a more
appropriate numerical method or reduce `dt_h`. Model modularity does not itself
guarantee numerical stability.

## Add a pathogen species

A species combines a pathogen kind, a registered model and model-specific
parameters:

```python
from ais.models import PathogenKind, PathogenSpecies


pseudomonas = PathogenSpecies(
    name="Pseudomonas aeruginosa",
    kind=PathogenKind.BACTERIUM,
    model="bacterial-logistic",
    virulence=0.9,
    immune_evasion=0.35,
    parameters={
        "growth_rate_h": 0.55,
        "carrying_capacity_per_ml": 2.0e8,
        "natural_clearance_h": 0.005,
    },
)
```

Register and inject it:

```python
from ais.scenarios import create_baseline_simulation


simulation = create_baseline_simulation(seed=42)
simulation.register_species(pseudomonas)

simulation.inject(
    pseudomonas.name,
    25_000,
    location="lungs",
)
simulation.run(96.0, dt_h=0.1)
```

Species names must be unique within a simulation.

## Add anatomical structures

Add `Compartment` objects and connect them with directed `Vessel` edges through
`Environment.connect()`.

A pair of directed edges represents bidirectional circulation. The convenience
parameter `bidirectional=True` creates both directions using the same
coefficients.

```python
from ais.environment import (
    Compartment,
    CompartmentKind,
    Environment,
)


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
        name="kidney",
        kind=CompartmentKind.ORGAN,
        volume_ml=300.0,
        susceptible_cells=1.5e8,
    )
)

environment.add_compartment(
    Compartment(
        name="renal-lymph",
        kind=CompartmentKind.LYMPH,
        volume_ml=100.0,
        susceptible_cells=1.0e7,
    )
)

environment.connect(
    "blood",
    "kidney",
    transfer_rate_h=0.08,
    cell_migration_weight=1.0,
    bidirectional=True,
)

environment.connect(
    "kidney",
    "renal-lymph",
    transfer_rate_h=0.02,
    cell_migration_weight=1.3,
)

environment.connect(
    "renal-lymph",
    "blood",
    transfer_rate_h=0.03,
    cell_migration_weight=1.0,
)

environment.validate()
```

Construct a simulation from the custom environment:

```python
from ais.simulation import Simulation


simulation = Simulation(
    environment,
    seed=42,
)
```

### Connection semantics

`transfer_rate_h` defines the first-order pathogen transport rate per hour. For
a step of length `dt_h`, the transported fraction is calculated as:

```text
fraction = 1 - exp(-transfer_rate_h × dt_h)
```

`cell_migration_weight` is a relative migration weight, not a probability.
During migration, outgoing weights are normalized across available
connections and modified by target-compartment inflammation.

## Create a custom scenario

A scenario factory should construct and return a fully configured simulation.
This makes experimental conditions reproducible and testable.

```python
from ais.immune import Macrophage, Neutrophil
from ais.models import PathogenKind, PathogenSpecies
from ais.scenarios import create_human_environment
from ais.simulation import Simulation


def create_sepsis_scenario(seed: int | None = 42) -> Simulation:
    simulation = Simulation(
        create_human_environment(),
        seed=seed,
    )

    bacterium = PathogenSpecies(
        name="Sepsis bacterium",
        kind=PathogenKind.BACTERIUM,
        model="bacterial-logistic",
        virulence=1.0,
        immune_evasion=0.25,
        parameters={
            "growth_rate_h": 0.9,
            "carrying_capacity_per_ml": 2.0e8,
            "natural_clearance_h": 0.005,
        },
    )

    simulation.register_species(bacterium)
    simulation.add_cells(
        Macrophage,
        count=100,
        location="peripheral-tissue",
    )
    simulation.add_cells(
        Neutrophil,
        count=250,
        location="blood",
    )

    simulation.inject(
        bacterium.name,
        1_000_000,
        location="blood",
    )

    return simulation
```

Execute the scenario:

```python
simulation = create_sepsis_scenario(seed=42)
simulation.run(24.0, dt_h=0.05)

print(simulation.snapshot())
```

## Performance characteristics

Let:

- `L` be the number of anatomical compartments;
- `S` be the number of active pathogen strains;
- `C` be the number of immune-cell agents;
- `E` be the number of anatomical edges;
- `d̄` be the average outgoing graph degree.

The approximate computational cost per simulation step is:

```text
O(L × S + active pathogen edges + C × d̄)
```

More explicitly:

- local pathogen dynamics scale with active compartment-strain pairs;
- pathogen transport scales with outgoing edges of infected compartments;
- immune-cell behavior scales approximately linearly with the number of agents;
- immune-cell migration additionally depends on the local graph degree;
- snapshots and full-state serialization scale with the amount of retained
  simulation state.

Memory use is approximately:

```text
O(L + E + active compartment-strain pairs + C + event history)
```

For simulations requiring millions of immune-cell equivalents, individual
agents should be replaced by cohort agents carrying an integer or floating-point
multiplicity.

A cohort design might represent one object as:

```python
from dataclasses import dataclass

from ais.immune.base import ImmuneCell


@dataclass(slots=True)
class ImmuneCellCohort(ImmuneCell):
    multiplicity: float = 1.0
```

Cell actions must then scale production, clearance and death by
`multiplicity`. Cohort splitting or redistribution may be required during
migration, differentiation or stochastic mortality.

## Testing extensions

Every new model should include tests for:

1. parameter validation;
2. non-negative population invariants;
3. deterministic behavior under a fixed seed;
4. expected behavior in the absence of immune pressure;
5. expected behavior under immune pressure;
6. numerical stability at supported time steps;
7. compatibility with pathogen transport;
8. behavior at zero population and zero signal;
9. conservation properties where biologically applicable.

Run the complete test suite with:

```bash
uv run pytest
```

Run a specific test module with:

```bash
uv run pytest tests/test_simulation.py
```

Run with coverage reporting:

```bash
uv run pytest \
    --cov=ais \
    --cov-report=term-missing
```

## Scientific interpretation

The framework provides a mechanism-oriented simulation architecture, not a
validated physiological digital twin. New mechanisms should document:

- biological assumptions;
- units;
- parameter sources;
- integration method;
- validity range;
- conservation assumptions;
- known omissions;
- calibration and validation procedures.

A modular implementation improves maintainability and experimental isolation,
but it does not establish biological validity.