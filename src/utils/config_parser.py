from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

import yaml


class ConfigurationError(ValueError):
    """Raised when an AIS configuration is syntactically or semantically invalid."""


class PathogenType(StrEnum):
    """Pathogen categories supported by the baseline simulation."""

    BACTERIA = "bacteria"
    VIRUS = "virus"

    @property
    def species_name(self) -> str:
        """Return the simulation species associated with this pathogen type."""
        species_by_pathogen = {
            PathogenType.BACTERIA: "E. coli",
            PathogenType.VIRUS: "Influenza A",
        }
        return species_by_pathogen[self]


@dataclass(frozen=True, slots=True)
class SimulationConfig:
    """Temporal and stochastic simulation parameters."""

    seed: int = 42
    duration_h: float = 72.0
    dt_h: float = 0.1

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> SimulationConfig:
        _reject_unknown_keys(
            values,
            allowed={"seed", "duration_h", "dt_h"},
            section="simulation",
        )

        seed = _read_int(values, "seed", default=42)
        duration_h = _read_float(values, "duration_h", default=72.0)
        dt_h = _read_float(values, "dt_h", default=0.1)

        if duration_h < 0.0:
            raise ConfigurationError(
                "'simulation.duration_h' must be greater than or equal to zero"
            )

        if dt_h <= 0.0:
            raise ConfigurationError("'simulation.dt_h' must be greater than zero")

        return cls(
            seed=seed,
            duration_h=duration_h,
            dt_h=dt_h,
        )


@dataclass(frozen=True, slots=True)
class PathogenConfig:
    """Parameters describing the initial pathogen injection."""

    type: PathogenType = PathogenType.BACTERIA
    amount: float = 100_000.0
    location: str = "blood"

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> PathogenConfig:
        _reject_unknown_keys(
            values,
            allowed={"type", "amount", "location"},
            section="pathogen",
        )

        pathogen_type_value = _read_string(
            values,
            "type",
            default=PathogenType.BACTERIA.value,
        )

        try:
            pathogen_type = PathogenType(pathogen_type_value)
        except ValueError as error:
            supported = ", ".join(item.value for item in PathogenType)
            raise ConfigurationError(f"'pathogen.type' must be one of: {supported}") from error

        amount = _read_float(values, "amount", default=100_000.0)
        location = _read_string(values, "location", default="blood")

        if amount < 0.0:
            raise ConfigurationError("'pathogen.amount' must be greater than or equal to zero")

        if not location.strip():
            raise ConfigurationError("'pathogen.location' must not be empty")

        return cls(
            type=pathogen_type,
            amount=amount,
            location=location,
        )


@dataclass(frozen=True, slots=True)
class AISConfig:
    """Root configuration for an artificial immune-system simulation."""

    simulation: SimulationConfig
    pathogen: PathogenConfig

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> AISConfig:
        _reject_unknown_keys(
            values,
            allowed={"simulation", "pathogen"},
            section="root",
        )

        return cls(
            simulation=SimulationConfig.from_mapping(
                _read_mapping(values, "simulation", default={})
            ),
            pathogen=PathogenConfig.from_mapping(_read_mapping(values, "pathogen", default={})),
        )


def find_project_root(start: str | Path | None = None) -> Path:
    """
    Locate the project root by searching upward for ``pyproject.toml``.

    The search defaults to the location of this module rather than the current
    working directory. Consequently, changing the process working directory
    does not affect configuration resolution.

    Args:
        start:
            Optional file or directory from which to begin searching.

    Returns:
        The directory containing ``pyproject.toml``.

    Raises:
        FileNotFoundError:
            If no project root can be identified.
    """
    candidate = Path(start).resolve() if start is not None else Path(__file__).resolve()

    if candidate.is_file():
        candidate = candidate.parent

    for directory in (candidate, *candidate.parents):
        if (directory / "pyproject.toml").is_file():
            return directory

    raise FileNotFoundError(
        f"Could not locate a project root from '{candidate}': "
        "no parent directory contains pyproject.toml"
    )


def resolve_config_path(path: str | Path | None = None) -> Path:
    """
    Resolve an explicit or default AIS configuration path.

    Relative explicit paths are interpreted relative to the current working
    directory, following normal command-line conventions. If no path is given,
    ``config/ais.yaml`` is resolved relative to the project root.
    """
    if path is None:
        return find_project_root() / "config" / "ais.yaml"

    return Path(path).expanduser().resolve()


def load_config(path: str | Path | None = None) -> AISConfig:
    """
    Read and validate an AIS YAML configuration.

    Args:
        path:
            Optional configuration path. If omitted, the parser loads
            ``config/ais.yaml`` from the project root.

    Raises:
        FileNotFoundError:
            If the resolved configuration file does not exist.
        ConfigurationError:
            If the file cannot be read or its contents are invalid.
    """
    config_path = resolve_config_path(path)

    if not config_path.is_file():
        raise FileNotFoundError(f"AIS configuration file not found: {config_path}")

    try:
        with config_path.open("r", encoding="utf-8") as file:
            document = yaml.safe_load(file)
    except yaml.YAMLError as error:
        raise ConfigurationError(
            f"Could not parse YAML configuration '{config_path}': {error}"
        ) from error
    except OSError as error:
        raise ConfigurationError(
            f"Could not read configuration file '{config_path}': {error}"
        ) from error

    if document is None:
        document = {}

    if not isinstance(document, Mapping):
        raise ConfigurationError("The YAML document root must be a mapping")

    return AISConfig.from_mapping(document)


def _read_mapping(
    values: Mapping[str, Any],
    key: str,
    *,
    default: Mapping[str, Any],
) -> Mapping[str, Any]:
    value = values.get(key, default)

    if not isinstance(value, Mapping):
        raise ConfigurationError(f"'{key}' must be a mapping")

    return value


def _read_string(
    values: Mapping[str, Any],
    key: str,
    *,
    default: str,
) -> str:
    value = values.get(key, default)

    if not isinstance(value, str):
        raise ConfigurationError(f"'{key}' must be a string")

    return value


def _read_int(
    values: Mapping[str, Any],
    key: str,
    *,
    default: int,
) -> int:
    value = values.get(key, default)

    # bool is an int subclass and must therefore be excluded explicitly.
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"'{key}' must be an integer")

    return value


def _read_float(
    values: Mapping[str, Any],
    key: str,
    *,
    default: float,
) -> float:
    value = values.get(key, default)

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"'{key}' must be numeric")

    return float(value)


def _reject_unknown_keys(
    values: Mapping[str, Any],
    *,
    allowed: set[str],
    section: str,
) -> None:
    unknown = set(values) - allowed

    if unknown:
        formatted_keys = ", ".join(sorted(map(str, unknown)))
        raise ConfigurationError(
            f"Unknown key(s) in configuration section '{section}': {formatted_keys}"
        )
