# Artificial Immune System

A hybrid agent-based and population-based artificial immune-system simulator
written in Python.

## Project objective

The project aims to provide a modular simulation framework that can be extended
without modifying its core execution engine. Biological behavior is separated
into independently replaceable environment, pathogen, immune-cell and simulation
components.

The principal extension mechanisms are:

- **Environment topology:** Add organs, tissues, blood compartments, lymphatic
  compartments and directed vascular connections.
- **Pathogen models:** Add new bacterial, viral, fungal, parasitic or synthetic
  pathogen dynamics by implementing the `PathogenModel` interface.
- **Immune-cell models:** Add new innate or adaptive immune-cell phenotypes by
  subclassing `ImmuneCell`.
- **Pathogen species:** Configure distinct strains using immutable
  `PathogenSpecies` definitions and model-specific parameters.
- **Signals and interactions:** Introduce additional cytokines, chemokines,
  antibodies and intercellular signaling mechanisms.
- **Scenarios:** Construct reproducible anatomical configurations, baseline
  immune populations and infection protocols independently of the core model.

The architecture is intended to support both experimental model development and
comparative computational studies. New biological mechanisms should generally
be implemented through the public extension interfaces rather than by modifying
`Simulation`.

## Modeling strategy

Pathogens are represented as aggregate populations indexed by anatomical
compartment and strain. Immune cells are represented as explicit agents.

This hybrid representation avoids allocating one Python object per bacterium or
virion, which would be computationally prohibitive, while retaining heterogeneous
state and behavior for immune-cell agents.

The model is qualitative and mechanistic. Its default coefficients are
illustrative and have not been calibrated as clinical parameters. Simulation
output must therefore not be interpreted as medical prediction.

## Project structure

```text
artificial-immune-system/
├── pyproject.toml
├── README.md
├── README_EXAMPLES.md
├── src/
│   └── ais/
│       ├── __init__.py
│       ├── environment.py
│       ├── models.py
│       ├── pathogens.py
│       ├── simulation.py
│       ├── scenarios.py
│       ├── cli.py
│       └── immune/
│           ├── __init__.py
│           ├── base.py
│           └── cells.py
└── tests/
    └── test_simulation.py
```

## Requirements

- Python 3.11 or newer
- [`uv`](https://docs.astral.sh/uv/)

## Environment setup with `uv`

Install a `uv`-managed Python 3.11 interpreter and pin it for the project:

```bash
uv python install 3.11
uv python pin 3.11
```

Create the virtual environment and synchronize all runtime and development
dependencies:

```bash
uv venv --python 3.11 .venv
uv sync --extra dev
```

The explicit `uv venv` step is optional because `uv sync` creates the project
environment when necessary. The minimal setup is therefore:

```bash
uv python install 3.11
uv python pin 3.11
uv sync --extra dev
```

Verify the interpreter and installation:

```bash
uv run python --version
uv run pytest
```

Manual activation is normally unnecessary because `uv run` executes commands
inside the project environment.

On Linux or macOS:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

## Testing and static analysis

Run the test suite:

```bash
uv run pytest
```

Run linting and formatting checks:

```bash
uv run ruff check .
uv run ruff format --check .
```

Run static type analysis:

```bash
uv run mypy
```

Apply automatic formatting:

```bash
uv run ruff format .
```

## Usage and extension documentation

See [`README_EXAMPLES.md`](README_EXAMPLES.md) for:

- bacterial and viral infection examples;
- programmatic simulation control;
- custom immune-cell phenotypes;
- custom pathogen dynamics;
- pathogen-species definitions;
- custom anatomical environments;
- performance characteristics.

## Architectural components

### `Environment`

`Environment` is a directed anatomical graph. Its nodes are `Compartment`
objects and its edges are `Vessel` objects.

This representation permits arbitrary anatomical topologies without coupling
the environment to pathogen or immune-cell implementations.

### `PathogenModel`

`PathogenModel` defines local pathogen dynamics. Implementations may represent,
for example:

- extracellular bacterial growth;
- target-cell-limited viral replication;
- intracellular bacterial replication;
- fungal growth;
- parasitic life cycles;
- antimicrobial resistance;
- latency and reactivation.

Transport between compartments remains an independent simulation concern.

### `ImmuneCell`

`ImmuneCell` is the base abstraction for immune-cell agents. Subclasses define
phenotype-specific behavior through `act()`, while the simulation engine handles
generic concerns such as:

- aging;
- mortality;
- anatomical location;
- graph-based migration;
- deterministic random-number generation.

### `Simulation`

`Simulation` coordinates pathogen transport, local pathogen dynamics,
immune-cell actions, antibody activity, tissue state and systemic signal decay.

Domain-specific behavior should be implemented through registered models and
cell types rather than embedded directly into the simulation loop.

## Reproducibility

Supply a fixed seed when constructing a simulation:

```python
simulation = Simulation(environment, seed=42)
```

The seeded NumPy random-number generator controls stochastic migration and other
future stochastic mechanisms. Reproducibility additionally requires identical:

- source code;
- dependency lock file;
- Python version;
- initial conditions;
- time-step configuration.

Commit the following files:

```text
.python-version
pyproject.toml
uv.lock
```

Do not commit the virtual environment:

```gitignore
.venv/
__pycache__/
.pytest_cache/
.ruff_cache/
.mypy_cache/
*.py[cod]
```

## Scientific limitations

The current implementation intentionally simplifies several biological
mechanisms:

1. Pathogens are compartment-level populations rather than individual agents.
2. Viral infected-cell populations are pooled per compartment.
3. Cytokines are scalar compartment-level signals rather than
   reaction-diffusion fields.
4. Antigen presentation and antibodies are globally indexed by strain.
5. Hemodynamics are approximated using first-order transfer rates on directed
   graph edges.
6. Immune agents may represent cellular cohorts rather than literal individual
   cells.
7. Default coefficients are illustrative and are not experimentally calibrated.
8. The simulation is not a diagnostic, therapeutic or clinical prediction
   system.

These are model assumptions rather than architectural constraints. The modular
interfaces are intended to permit progressively more detailed implementations.