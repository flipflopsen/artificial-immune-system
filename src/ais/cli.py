from __future__ import annotations

import argparse
import json

from ais.scenarios import create_baseline_simulation


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an artificial immune-system simulation"
    )
    parser.add_argument(
        "--pathogen",
        choices=("bacteria", "virus"),
        default="bacteria",
    )
    parser.add_argument("--amount", type=float, default=100_000.0)
    parser.add_argument("--location", default="blood")
    parser.add_argument("--duration", type=float, default=72.0)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    simulation = create_baseline_simulation(seed=arguments.seed)

    species_name = (
        "E. coli" if arguments.pathogen == "bacteria" else "Influenza A"
    )

    simulation.inject(
        species_name,
        arguments.amount,
        location=arguments.location,
    )
    simulation.run(arguments.duration, dt_h=arguments.dt)

    print(json.dumps(simulation.snapshot(), indent=2))


if __name__ == "__main__":
    main()