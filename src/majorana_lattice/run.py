"""Running a Program: the only place anything is executed.

An executor is any function (stim.Circuit, shots) -> bool array of shape
(shots, num_measurements), one row of measurement records per shot. Local
Stim is the default. A program can carry its own executor (prepare(...,
executor=...) or p.with_executor(...)), and run can override it.

    res = run(p, shots=1000)
    res["v"]                  # a named readout, +1/-1 per shot
    res.pairs()               # every pair's value per shot, without measuring it
    res.summary(at=p_mid)     # fermions, e and m per shot, as they were at p_mid
"""
from __future__ import annotations

import numpy as np
import stim

from . import queries as Q
from .lattice import Vertex
from .program import Executor, Pair, Program, present


def stim_executor(seed: int | None = None) -> Executor:
    """Sample with Stim's own sampler, locally."""
    def execute(circuit: stim.Circuit, shots: int) -> np.ndarray:
        return circuit.compile_sampler(seed=seed).sample(shots)
    return execute


def run(p: Program, shots: int = 1000, executor: Executor | None = None,
        seed: int | None = None) -> Result:
    """Execute the program and return its results. Uses `executor` if given,
    else the program's own, else local Stim (seeded with `seed`)."""
    executor = executor or p.executor or stim_executor(seed)
    bits = np.asarray(executor(p.circuit, shots), dtype=bool)
    if bits.shape[1] != p.num_records:
        raise ValueError(f"executor returned {bits.shape[1]} records per shot, expected {p.num_records}")
    return Result(p, bits)


class Result:
    """Measurement records of every shot, and the values they imply."""

    def __init__(self, program: Program, bits: np.ndarray):
        self.program = program
        self.bits = bits

    @property
    def shots(self) -> int:
        return self.bits.shape[0]

    def __getitem__(self, name: str) -> np.ndarray:
        """A named readout (or named link outcome), +1/-1 per shot."""
        return 1 - 2 * self.bits[:, self.program.readouts[name]].astype(np.int8)

    def outcomes(self) -> np.ndarray:
        """Every measurement record as +1/-1, shape (shots, records)."""
        return 1 - 2 * self.bits.astype(np.int8)

    def value(self, pair: Pair) -> np.ndarray | None:
        """A tracked pair's value per shot (None if undetermined)."""
        if pair.records is None:
            return None
        par = self.bits[:, sorted(pair.records)].sum(axis=1) % 2 if pair.records else np.zeros(self.shots, int)
        return (1 - 2 * (par ^ int(pair.flip))).astype(np.int8)

    def pairs(self, at: Program | None = None) -> dict[tuple, np.ndarray | None]:
        """Every pair's value per shot, keyed by its two ends. `at` is an
        earlier Program this one extends: values as they were then."""
        at = self._at(at)
        return {tuple(present(at.lattice, Vertex, q) for q in sorted((P.a, P.b))): self.value(P)
                for P in at.pairs}

    def plaquettes(self, at: Program | None = None) -> np.ndarray:
        """Plaquette values, shape (shots, faces); 0 where undetermined."""
        v = np.array(Q.plaquette_values(self._at(at)), dtype=np.int8)
        return np.broadcast_to(v, (self.shots, len(v)))

    def summary(self, at: Program | None = None) -> dict:
        """Fermions, e and m anyons per shot, plus the Majoranas."""
        at = self._at(at)
        vals = [v for v in self.pairs(at).values() if v is not None]
        col, ok = Q.colouring(at)
        plaq = Q.plaquette_values(at)
        return {
            "fermions": np.sum(np.array(vals) == -1, axis=0) if vals else np.zeros(self.shots, int),
            "e": sum(1 for f, x in enumerate(plaq) if x == -1 and col[f] == 0),
            "m": sum(1 for f, x in enumerate(plaq) if x == -1 and col[f] == 1),
            "majoranas": Q.majoranas(at),
            "colouring_consistent": ok,
        }

    def _at(self, at: Program | None) -> Program:
        if at is None:
            return self.program
        if not self.program.extends(at):
            raise ValueError("`at` must be an earlier stage of the program that was run")
        return at
