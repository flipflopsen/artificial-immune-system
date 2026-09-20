from __future__ import annotations

import argparse
import json
from pathlib import Path

from ais.scenarios import create_baseline_simulation
from utils.config_parser import load_config


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run an artificial immune-system simulation")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("./config/ais.yaml"),
        help=("Path to the YAML configuration file. Defaults to <project-root>/config/ais.yaml."),
    )
    return parser.parse_args()


def main() -> None:
    """Load the configuration and execute the simulation."""
    arguments = parse_arguments()
    config = load_config(arguments.config)

    simulation = create_baseline_simulation(
        seed=config.simulation.seed,
    )

    simulation.inject(
        config.pathogen.type.species_name,
        config.pathogen.amount,
        location=config.pathogen.location,
    )

    simulation.run(
        config.simulation.duration_h,
        dt_h=config.simulation.dt_h,
    )

    print(json.dumps(simulation.snapshot(), indent=2))


if __name__ == "__main__":
    main()
