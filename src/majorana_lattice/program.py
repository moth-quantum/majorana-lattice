"""Programs: matching-code experiments as Stim circuits, built without running them.

A Program is an immutable value: the circuit so far, plus the classical
bookkeeping needed to decide what can be done next (who is paired with whom,
which plaquettes hold anyons, and what each measurement record means).
Operations in `ops` take a Program and return a new one; nothing executes
until `run`.

Nothing in the bookkeeping depends on measurement outcomes. A pair's value is
instead stored as an expression: the parity of some measurement records,
times a fixed sign. Measuring a link (i, j) between pairs A and B gives

  * the dimer (i, j), whose value is the new record k;
  * the pair (a, b), whose value is value(A) * value(B) * record k.

Applying a Pauli flips the sign of every pair and plaquette it anticommutes
with. The paper's fix for a -1 outcome is a Pauli conditioned on record k,
so it toggles k in the expressions it touches. Every tracked pair stays
determined by the records, so `run` can report every pair's value in every
shot without measuring it.

The vacuum (every dimer and every plaquette reading +1) is prepared the way
hardware would: measure every plaquette, then undo each -1 with a Pauli
conditioned on that record (the plaquette's destabilizer, from Stim). No
postselection is needed, and the plaquette records are kept as
("prepare", f) at the start of the measurement record.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import stim

from .lattice import Lattice, Vertex


def link_pauli(lat: Lattice, e: int) -> stim.PauliString:
    """The link operator on edge e: its label's Pauli on both ends."""
    E = lat.edges[e]
    p = stim.PauliString(lat.n)
    p[E.u] = E.label.upper()
    p[E.v] = E.label.upper()
    return p


def simplify_path(path: list[int]) -> list[int]:
    """Cancel back-and-forth steps (an edge used twice in a row)."""
    out: list[int] = []
    for e in path:
        if out and out[-1] == e:
            out.pop()
        else:
            out.append(e)
    return out


@dataclass(frozen=True, eq=False)
class Pair:
    """Two vertices whose parity operator is tracked: a dimer, a longer
    string, or a released pair of Majoranas."""
    a: int
    b: int
    path: tuple[int, ...]          # edge ids from a to b
    op: stim.PauliString           # signed, Hermitian
    kind: str = "string"           # "string" (stabilizer) or "majorana"
    flip: bool = False             # value = (-1)^(flip + parity of records)
    records: frozenset[int] | None = frozenset()   # None: not determined

    def other(self, q: int) -> int:
        return self.b if q == self.a else self.a

    def path_from(self, q: int) -> tuple[int, ...]:
        return self.path if q == self.a else self.path[::-1]

    @property
    def is_dimer(self) -> bool:
        return self.kind == "string" and len(self.path) == 1


@dataclass(frozen=True, eq=False)
class Program:
    lattice: Lattice
    background: str
    prep: stim.Circuit                          # shared between Programs; never mutated
    instructions: tuple[stim.CircuitInstruction, ...]
    pairs: tuple[Pair, ...]
    owner: tuple[int, ...]                      # vertex -> index into pairs
    plaq_ops: tuple[stim.PauliString, ...]
    plaq_flips: tuple[bool | None, ...]         # None: not determined
    measurements: tuple[tuple, ...]             # per record: ("link", e) or ("readout", name)
    readouts: dict[str, int] = field(default_factory=dict)

    @property
    def circuit(self) -> stim.Circuit:
        """The whole program as a stim.Circuit (a fresh copy each time)."""
        c = self.prep.copy()
        for inst in self.instructions:
            c.append(inst)
        return c

    @property
    def num_records(self) -> int:
        return len(self.measurements)

    def pair_at(self, v: Vertex | int) -> Pair:
        """The pair holding vertex v."""
        return self.pairs[self.owner[resolve(self.lattice, v, Vertex)]]

    def extends(self, other: Program) -> bool:
        """True if this program starts with everything `other` did."""
        n = len(other.instructions)
        return (self.lattice is other.lattice and self.prep is other.prep
                and len(self.instructions) >= n and self.instructions[:n] == other.instructions)

    @classmethod
    def from_circuit(cls, lattice: Lattice, circuit: stim.Circuit) -> Program:
        """Rebuild a Program from its circuit (as made by Program.circuit)."""
        from .ops import from_circuit
        return from_circuit(lattice, circuit)

    def _evolve(self, **kw) -> Program:
        return replace(self, **kw)


def prepare(lattice: Lattice, background: str = "z") -> Program:
    """The vacuum: every `background` link a dimer, every dimer and plaquette +1."""
    prep, plaq_ops = _vacuum(lattice, background)
    pairs, owner, seen = [], [-1] * lattice.n, set()
    for q in range(lattice.n):
        e = lattice.edge_with_label(q, background)
        if e in seen:
            continue
        seen.add(e)
        E = lattice.edges[e]
        owner[E.u] = owner[E.v] = len(pairs)
        pairs.append(Pair(E.u, E.v, (e,), link_pauli(lattice, e)))
    return Program(lattice, background, prep, (), tuple(pairs), tuple(owner),
                   plaq_ops, (False,) * len(plaq_ops),
                   tuple(("prepare", f) for f in range(len(plaq_ops))))


def _vacuum(lattice: Lattice, background: str):
    """Prep circuit and signed plaquette operators, cached on the lattice.

    Reset every qubit, rotate so every `background` link reads +1 (nothing to
    do for z), measure every plaquette, and apply each measured plaquette's
    destabilizer if it read -1. The destabilizers come from a Stim tableau
    whose stabilizers are the dimers, the plaquettes, and the products of
    single-qubit background Paulis that commute with every plaquette (the
    logicals the start state already fixes). Those must be pure products of
    the background Pauli: a correction that commuted only with a logical times
    a plaquette would flip the logical."""
    cache = lattice.__dict__.setdefault("_vacuum_cache", {})
    if background in cache:
        return cache[background]
    n, B = lattice.n, background.upper()
    plaqs = []
    for f in lattice.faces:
        p = stim.PauliString(n)
        for q, P in f.pauli.items():
            p[q] = P
        plaqs.append(p)
    dimers = [link_pauli(lattice, e) for e, E in enumerate(lattice.edges) if E.label == background]
    # products of B that commute with every plaquette: measure, from |0..0>,
    # an X on each qubit where the plaquette anticommutes with B; the pure-Z
    # stabilizers left over are exactly those products, with Z standing for B
    sim = stim.TableauSimulator()
    for p in plaqs:
        xp = stim.PauliString(n)
        for q in p.pauli_indices():
            if "_XYZ"[p[q]] != B:
                xp[q] = "X"
        sim.measure_observable(xp)
    logicals = []
    for g in sim.canonical_stabilizers():
        if not g.pauli_indices("XY"):
            b = stim.PauliString(n)
            for q in g.pauli_indices():
                b[q] = B
            logicals.append(b)
    T = stim.Tableau.from_stabilizers(dimers + plaqs + logicals, allow_redundant=True)
    row = {str(T.z_output(k)): k for k in range(n)}

    # written as text: Stim parses it far faster than per-target Python calls
    lines = [f"R[prepare:{background}] " + " ".join(map(str, range(n)))]
    rotate = {"x": "H", "y": "H_YZ"}.get(background)
    if rotate:
        lines.append(f"{rotate}[prepare:{background}] " + " ".join(map(str, range(n))))
    for f, p in enumerate(plaqs):       # distinct tags stop Stim merging these
        lines.append(f"MPP[prepare:{background}:{f}] " + "*".join(
            f"{'_XYZ'[p[q]]}{q}" for q in p.pauli_indices()))
    F = len(plaqs)
    for f, p in enumerate(plaqs):
        k = row.get(str(p))
        if k is None:                   # fixed by the other plaquettes
            continue
        C = T.x_output(k)
        for P in "XYZ":
            qs = C.pauli_indices(P)
            if qs:
                lines.append(f"C{P}[prepare:{background}:{f}] "
                             + " ".join(f"rec[{f - F}] {q}" for q in qs))
    prep = stim.Circuit("\n".join(lines))
    # Every tableau stabilizer ends at +1. A plaquette fixed by the others is
    # +/- a product of them, so it may end at -1: flip its sign so vacuum
    # still reads +1.
    signed = []
    for p in plaqs:
        if str(p) in row:
            signed.append(p)
            continue
        prod = stim.PauliString(n)
        for k in range(n):
            if not p.commutes(T.x_output(k)):
                prod *= T.z_output(k)
        assert prod in (p, -p)
        signed.append(prod)
    cache[background] = (prep, tuple(signed))
    return cache[background]


def _mpp_targets(op: stim.PauliString) -> list:
    """MPP targets for a signed Pauli product."""
    targets = []
    for q in op.pauli_indices():
        if targets:
            targets.append(stim.target_combiner())
        P = "_XYZ"[op[q]]
        make = {"X": stim.target_x, "Y": stim.target_y, "Z": stim.target_z}[P]
        targets.append(make(q, invert=(len(targets) == 0 and op.sign == -1)))
    return targets


def resolve(lattice: Lattice, x, kind: type) -> int:
    """Index of an address of the given kind; an int is passed through."""
    if isinstance(x, bool):
        raise TypeError(f"expected a {kind.__name__} or an int, got {x!r}")
    if isinstance(x, int):
        return x
    if not isinstance(x, kind):
        raise TypeError(f"expected a {kind.__name__}, got {x!r}")
    return lattice.index(x)


def present(lattice: Lattice, kind: type, k: int):
    """An address when the lattice has them, otherwise the index."""
    return lattice.address(kind, k) if lattice.cells is not None else k


__all__ = ["Pair", "Program", "link_pauli", "prepare"]
