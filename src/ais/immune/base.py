from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar
from uuid import uuid4

if TYPE_CHECKING:
    from ais.simulation import Simulation


@dataclass(slots=True)
class ImmuneCell(ABC):
    """
    Base class for independently simulated immune-cell agents.

    Subclasses should contain phenotype-specific behavior in `act()`. Generic
    aging and migration are handled by Simulation.
    """

    location: str
    identifier: str = field(default_factory=lambda: uuid4().hex)
    age_h: float = 0.0
    activation: float = 0.0

    cell_type: ClassVar[str] = "immune-cell"
    lifespan_h: ClassVar[float] = 24.0 * 7.0
    mobility_h: ClassVar[float] = 0.05

    @abstractmethod
    def act(self, simulation: Simulation, dt_h: float) -> None:
        raise NotImplementedError

    def activate(self, stimulus: float, dt_h: float) -> None:
        stimulus = min(1.0, max(0.0, stimulus))
        relaxation_rate = 0.5
        self.activation += relaxation_rate * (stimulus - self.activation) * dt_h
        self.activation = min(1.0, max(0.0, self.activation))
