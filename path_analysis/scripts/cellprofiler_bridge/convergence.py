"""Per-species convergence tracking via SEM%.

Ported from the active `cellprofiler_test` workflow so the eventual merge
keeps one in-repo implementation.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


class ConvergenceTracker:
    """Track per-species measurement convergence via SEM%."""

    def __init__(self, min_cells: int = 30, max_cells: int = 200, sem_pct: float = 3.0):
        self.min_cells = min_cells
        self.max_cells = max_cells
        self.sem_pct = sem_pct
        self.species_values: dict[str, list[float]] = defaultdict(list)
        self.exempt_species: set[str] = set()

    def add(self, species: str, values: list[float]) -> None:
        self.species_values[species].extend(values)

    def add_cells(self, species: str, measurements: list[dict], key: str = "area_um2") -> None:
        for measurement in measurements:
            self.species_values[species].append(float(measurement[key]))

    def exempt(self, species: str) -> None:
        self.exempt_species.add(species)

    def _sem_pct_value(self, species: str) -> float:
        values = self.species_values.get(species, [])
        n_values = len(values)
        if n_values < 2:
            return float("inf")
        arr = np.array(values)
        mean = arr.mean()
        if mean <= 0:
            return float("inf")
        sem = arr.std(ddof=1) / np.sqrt(n_values)
        return 100 * sem / mean

    def truly_converged(self, species: str) -> bool:
        n_values = len(self.species_values.get(species, []))
        if n_values < self.min_cells:
            return False
        return self._sem_pct_value(species) < self.sem_pct

    def is_capped(self, species: str) -> bool:
        if species in self.exempt_species:
            return False
        n_values = len(self.species_values.get(species, []))
        return n_values >= self.max_cells and not self.truly_converged(species)

    def is_done(self, species: str) -> bool:
        if self.truly_converged(species):
            return True
        if species in self.exempt_species:
            return False
        n_values = len(self.species_values.get(species, []))
        return n_values >= self.max_cells

    def cells_needed(self, species: str) -> int:
        values = self.species_values.get(species, [])
        n_values = len(values)
        cap = None if species in self.exempt_species else self.max_cells
        if cap is not None and n_values >= cap:
            return 0
        if n_values < 2:
            return (cap or 10000) - n_values
        arr = np.array(values)
        mean = arr.mean()
        if mean <= 0:
            return (cap or 10000) - n_values
        if self._sem_pct_value(species) < self.sem_pct:
            return 0
        sd = arr.std(ddof=1)
        target_n = int(np.ceil((100 * sd / (mean * self.sem_pct)) ** 2))
        remaining = max(0, target_n - n_values)
        if cap is not None:
            remaining = min(remaining, cap - n_values)
        return remaining

    def cell_count(self, species: str) -> int:
        return len(self.species_values.get(species, []))

    def status(self, species: str) -> str:
        values = self.species_values.get(species, [])
        n_values = len(values)
        if n_values < 2:
            return f"n={n_values}"
        sem_pct = self._sem_pct_value(species)
        if self.truly_converged(species):
            tag = " CONVERGED"
        elif self.is_capped(species):
            tag = " CAPPED (not converged)"
        elif species in self.exempt_species:
            tag = " EXEMPT (no cap)"
        else:
            tag = ""
        arr = np.array(values)
        return f"n={n_values}, mean={arr.mean():.1f}, SEM%={sem_pct:.1f}%{tag}"

    def summary(self) -> str:
        return "\n".join(f"  {species}: {self.status(species)}" for species in sorted(self.species_values))

    def unconverged_species(self) -> list[str]:
        return [species for species in self.species_values if self.is_capped(species)]


def interleave_by_species(jobs: list[dict], species_key: str = "species") -> list[dict]:
    """Round-robin jobs across species for even sampling."""

    species_queues: dict[str, list[dict]] = defaultdict(list)
    for job in jobs:
        species_queues[job[species_key]].append(job)
    interleaved: list[dict] = []
    while any(species_queues.values()):
        for species in sorted(species_queues):
            if species_queues[species]:
                interleaved.append(species_queues[species].pop(0))
        species_queues = {key: value for key, value in species_queues.items() if value}
    return interleaved
