"""Matching codes on trivalent lattices, simulated with Stim.

Every vertex holds one Majorana mode. Vertices are paired up; each pair is
joined by a path of links, and the product of link operators along the path is
that pair's parity operator (Wootton 2015, arXiv:1501.07779).

  * A pair that is a stabilizer is a "string" (a "dimer" when it is one link).
    Value -1 means it holds a fermion.
  * A pair that has been released is a pair of computational Majoranas; its
    parity is an observable, not a stabilizer.
  * Plaquette value -1 means an e or m anyon, depending on the bicolouring.

How to read the board:
  * a dimer is a single-link pair that is a stabilizer;
  * every vertex not in a dimer is a Majorana (the ends of longer strings,
    and the ends of split dimers);
  * vacuum means every dimer and every plaquette reads +1. Plaquettes are
    prepared at +1 (not at a random first outcome), so a Majorana carried
    around an empty plaquette picks up no sign.

Stim holds the quantum state and is the source of truth for every value. This
module only decides which operators to track, and keeps each one's sign
fixed by history so that vacuum always reads +1.

Actions (all JSON-friendly, replayable from the log):
  measure_link(e)          measure the link operator; rewires the pairing
                           (outcome=+1/-1 postselects a random outcome)
  apply_link(e)            apply the link operator as a unitary
  apply_pauli(q, P)        apply X, Y or Z to one qubit
  release(q)               split a dimer into two Majoranas, or join them back
                           (bookkeeping only: no operation on the qubits)
  measure_label(L)         measure every link with label L (a Floquet round)
  move_majorana(q, e, fix) move the Majorana at q across link e, fixing a -1
                           with probability fix (default: the fix_rate field)
  move_fermion(q, e)       move the fermion on q's pair across link e
  move_anyon(f, e)         move the e or m anyon on plaquette f across link e
  route_majorana(a, b)     route a Majorana from a to b by moves

Each particle has a matching list of the moves on offer: majorana_moves(q),
fermion_moves(q) and anyon_moves(f), all as [(link, destination)].
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import stim

from .lattice import Lattice


def link_pauli(lat: Lattice, e: int) -> stim.PauliString:
    E = lat.edges[e]
    p = stim.PauliString(lat.n)
    p[E.u] = E.label.upper()
    p[E.v] = E.label.upper()
    return p


@dataclass
class Pair:
    a: int
    b: int
    path: list[int]            # edge ids from a to b
    op: stim.PauliString       # signed, Hermitian
    kind: str = "string"       # "string" (stabilizer) or "majorana"

    def other(self, q: int) -> int:
        return self.b if q == self.a else self.a

    def path_from(self, q: int) -> list[int]:
        return self.path if q == self.a else self.path[::-1]


@dataclass
class MatchingCode:
    lat: Lattice
    background: str = "z"
    seed: int | None = None
    fix_rate: float = 1.0      # default chance that a Majorana move fixes a -1 outcome
    plaquette_outcomes: list[int] | None = None  # replay: initial raw plaquette results
    sim: stim.TableauSimulator = field(init=False)
    pairs: dict[int, Pair] = field(init=False)
    owner: list[int] = field(init=False)     # vertex -> pair id
    plaq: list[stim.PauliString] = field(init=False)
    log: list[dict] = field(init=False)

    def __post_init__(self):
        self.sim = stim.TableauSimulator(seed=self.seed)
        self._fix_rng = random.Random(self.seed)   # separate from Stim's, so fixing never changes outcomes
        self.pairs, self.owner, self._next = {}, [-1] * self.lat.n, 0
        self.log = []
        # |0...0> is already +1 for every ZZ link. For another background
        # label, rotate each qubit so that label's link operators start at +1.
        for q in range(self.lat.n):
            if self.background == "x":
                self.sim.h(q)
            elif self.background == "y":
                self.sim.h_yz(q)
        seen = set()
        for q in range(self.lat.n):
            e = self.lat.edge_with_label(q, self.background)
            if e in seen:
                continue
            seen.add(e)
            E = self.lat.edges[e]
            self._add(Pair(E.u, E.v, [e], link_pauli(self.lat, e)))
        for q in range(self.lat.n):
            assert self.owner[q] >= 0
        # Plaquettes commute with every link operator. Vacuum means every
        # plaquette W_p = +(product of Paulis) reads +1: a Majorana carried
        # around an empty plaquette then picks up no sign. We postselect onto
        # that sector (on hardware this would be measurement plus correction).
        # Adopting random first outcomes instead leaves hidden flux that makes
        # Majorana fusion depend on the route taken.
        if self.plaquette_outcomes is None:
            self.plaquette_outcomes = [1] * len(self.lat.faces)
        self.plaq, self.plaquette_raw = [], []
        for k, f in enumerate(self.lat.faces):
            p = stim.PauliString(self.lat.n)
            for q, P in f.pauli.items():
                p[q] = P
            if self.plaquette_outcomes is not None and self.sim.peek_observable_expectation(p) == 0:
                self.sim.postselect_observable(p, desired_value=self.plaquette_outcomes[k] == -1)
            raw = -1 if self.sim.measure_observable(p) else 1
            self.plaquette_raw.append(raw)
            if raw == -1:
                p = -p
            self.plaq.append(p)
        assert all(v == 1 for v in self.values()["plaquettes"])

    # ------------------------------------------------------------ bookkeeping
    def _add(self, pair: Pair) -> int:
        pid = self._next
        self._next += 1
        self.pairs[pid] = pair
        self.owner[pair.a] = self.owner[pair.b] = pid
        return pid

    def _remove(self, pid: int) -> Pair:
        return self.pairs.pop(pid)

    def _measure(self, op: stim.PauliString, outcome: int | None = None) -> int:
        """Measure a Pauli observable. If the outcome is random and `outcome` is
        given, postselect on it; a deterministic outcome is never changed."""
        exp = self.sim.peek_observable_expectation(op)
        if exp != 0 or outcome is None:
            return -1 if self.sim.measure_observable(op) else 1
        self.sim.postselect_observable(op, desired_value=(outcome == -1))
        return outcome

    # ---------------------------------------------------------------- actions
    def measure_link(self, e: int, outcome: int | None = None) -> int:
        """Measure link e = (i, j). The pairs holding i and j are rewired:
        (i, a) and (j, b) become the dimer (i, j) and the pair (a, b), whose
        path runs a -> i -> j -> b. Returns the outcome (+1 or -1). Passing
        `outcome` postselects a random outcome on that value (replay uses this);
        a deterministic outcome is never changed."""
        E = self.lat.edges[e]
        i, j = E.u, E.v
        K = link_pauli(self.lat, e)
        pi, pj = self.owner[i], self.owner[j]
        if pi == pj:
            # i and j already paired: the link commutes with everything.
            old = self._remove(pi)
            m = self._measure(K, outcome)
            self._add(Pair(i, j, [e], K, "string"))
            self._record("measure_link", edge=e, outcome=m)
            return m
        A, B = self._remove(pi), self._remove(pj)
        a, b = A.other(i), B.other(j)
        new_op = A.op * B.op * K          # commutes with K; real sign
        assert new_op.sign in (1, -1)
        m = self._measure(K, outcome)
        kind = "string" if (A.kind == "string" and B.kind == "string") else "majorana"
        path = A.path_from(a) + [e] + B.path_from(j)
        self._add(Pair(i, j, [e], K, "string"))
        self._add(Pair(a, b, _simplify(path), new_op, kind))
        self._record("measure_link", edge=e, outcome=m)
        return m

    def apply_link(self, e: int):
        self.sim.do(link_pauli(self.lat, e))
        self._record("apply_link", edge=e)

    def apply_pauli(self, q: int, P: str):
        p = stim.PauliString(self.lat.n)
        p[q] = P.upper()
        self.sim.do(p)
        self._record("apply_pauli", qubit=q, pauli=P.upper())

    def release(self, q: int):
        """Stop treating the pair holding q as a stabilizer: its two ends
        become computational Majoranas (or re-pair them if already released)."""
        P = self.pairs[self.owner[q]]
        P.kind = "majorana" if P.kind == "string" else "string"
        self._record("release", qubit=q)

    def measure_label(self, label: str) -> list[int]:
        """Measure every link with this label: one round of a Floquet schedule."""
        out = []
        for e, E in enumerate(self.lat.edges):
            if E.label == label:
                out.append(self.measure_link(e))
        return out

    @classmethod
    def replay(cls, lat: Lattice, record: dict) -> "MatchingCode":
        """Rebuild a run from a log, forcing every random outcome to match.
        Raises if a deterministic outcome disagrees."""
        c = cls(lat, plaquette_outcomes=record["plaquette_raw"])
        for step in record["log"]:
            a = step["action"]
            if a == "measure_link":
                m = c.measure_link(step["edge"], outcome=step["outcome"])
                if m != step["outcome"]:
                    raise AssertionError(f"outcome mismatch at {step}")
            elif a == "apply_link":
                c.apply_link(step["edge"])
            elif a == "apply_pauli":
                c.apply_pauli(step["qubit"], step["pauli"])
            elif a == "release":
                c.release(step["qubit"])
        return c

    # ------------------------------------------------------- Majorana helpers
    @staticmethod
    def is_dimer(P: Pair) -> bool:
        return P.kind == "string" and len(P.path) == 1

    def majoranas(self) -> list[int]:
        """Vertices not in a dimer. Ends of a longer string are a Majorana
        pair with a definite parity; ends of a released pair are too."""
        return sorted(q for p in self.pairs.values() if not self.is_dimer(p) for q in (p.a, p.b))

    def majorana_moves(self, q: int) -> list[tuple[int, int]]:
        """Moves for the Majorana at q: [(link to measure, destination)]."""
        out = []
        for e in self.lat.incident[q]:
            j = self.lat.edges[e].other(q)
            Pj = self.pairs[self.owner[j]]
            if self.is_dimer(Pj) and self.owner[j] != self.owner[q]:
                out.append((e, Pj.other(j)))
        return out

    def move_majorana(self, q: int, e: int, fix: float | None = None) -> int:
        """Move the Majorana at q across link e; returns where it lands. On a
        -1 outcome, apply the old dimer's link operator to remove the fermion
        (the paper's fix) with probability `fix`, or `self.fix_rate` if not
        given (1 always fixes, 0 never does; True and False work too). Raises
        ValueError unless (e, destination) is in majorana_moves(q)."""
        fix = self.fix_rate if fix is None else fix
        k = next((k for ee, k in self.majorana_moves(q) if ee == e), None)
        if k is None:
            raise ValueError(f"the Majorana at {q} cannot move across link {e}")
        Pj = self.pairs[self.owner[self.lat.edges[e].other(q)]]
        old_edge = Pj.path_from(self.lat.edges[e].other(q))
        m = self.measure_link(e)
        if m == -1 and (fix >= 1 or (fix > 0 and self._fix_rng.random() < fix)):
            # the old pair j..k's operator anticommutes with both new pairs;
            # applying its links (= the operator up to a phase) removes both fermions
            for f in old_edge:
                self.apply_link(f)
        return k

    def _plan(self, src: int, dst: int, avoid) -> list[tuple[int, int]] | None:
        """Shortest route of moves on the current pairing: [(from, link)].

        Each move rewires the vertex it leaves and the vertex it moves across,
        so a valid route must never reuse a vertex, as a stop or as a vertex
        moved across. Dimers away from the route stay as they are now, so
        the current moves are valid along such a route."""
        from collections import deque
        # breadth-first over positions; each route carries its own set of
        # rewired vertices so it can check it never crosses its own wake
        best = {src: ((), frozenset((src,)))}
        dq = deque([src])
        while dq:
            q = dq.popleft()
            steps, used = best[q]
            if q == dst:
                return list(steps)
            for e, k in self.majorana_moves(q):
                j = self.lat.edges[e].other(q)
                if j in used or k in used or k in avoid or j in avoid or k in best:
                    continue
                best[k] = (steps + ((q, e),), used | {j, k})
                dq.append(k)
        return None

    def route_majorana(self, src: int, dst: int, avoid: set[int] = frozenset(), fix: float | None = None) -> list[int]:
        """Walk the Majorana at src to dst by moves, avoiding given vertices.
        Plans a route that never reuses a vertex, since every move rewires the
        dimers it passes. Returns the vertices visited."""
        if self.is_dimer(self.pairs[self.owner[src]]):
            raise ValueError(f"vertex {src} is not a Majorana")
        plan = self._plan(src, dst, avoid)
        if plan is None:
            raise ValueError(f"no route from {src} to {dst}")
        route = [src]
        for q0, e in plan:
            assert q0 == route[-1]
            route.append(self.move_majorana(q0, e, fix))
        assert route[-1] == dst
        return route

    # ----------------------------------------------------------- fermion helpers
    def fermion_moves(self, q: int) -> list[tuple[int, int]]:
        """Moves for the fermion on the pair holding q: [(link, vertex reached)].

        Applying a link operator toggles a fermion on each of the two pairs at
        its ends, so a fermion on one side moves across. Every link leaving
        the pair (from either end) to a different pair is listed."""
        P = self.pairs[self.owner[q]]
        out = []
        for end in (P.a, P.b):
            for e in self.lat.incident[end]:
                j = self.lat.edges[e].other(end)
                if self.owner[j] != self.owner[q]:
                    out.append((e, j))
        return out

    def move_fermion(self, q: int, e: int) -> int:
        """Move the fermion on q's pair across link e; returns the vertex reached.
        Raises ValueError unless (e, vertex) is in fermion_moves(q)."""
        j = next((j for ee, j in self.fermion_moves(q) if ee == e), None)
        if j is None:
            raise ValueError(f"the pair holding {q} has no move across link {e}")
        self.apply_link(e)
        return j

    # ------------------------------------------------------------ anyon helpers
    def anyon_moves(self, f: int) -> list[tuple[int, int]]:
        """Moves for the anyon on plaquette f: [(link to cross, destination face)].

        Applying a link's own Pauli to one endpoint flips exactly the two
        plaquettes either side of that link, so the anyon moves across it. Only
        moves that keep the anyon's type (both faces the same colour) onto a
        plaquette with no anyon are listed, so the number of e and of m anyons
        never changes and no fermion is made. Empty plaquettes have no moves."""
        v = self.values()["plaquettes"]
        if v[f] != -1:
            return []
        col, _ = self.colouring()
        return [(e, g) for g, e in self.lat.face_neighbours()[f] if col[g] == col[f] and v[g] == 1]

    def move_anyon(self, f: int, e: int) -> int:
        """Move the anyon on plaquette f across link e; returns the new face.
        Raises ValueError unless (e, face) is in anyon_moves(f)."""
        dest = next((g for ee, g in self.anyon_moves(f) if ee == e), None)
        if dest is None:
            raise ValueError(f"the anyon on plaquette {f} cannot move across link {e}")
        E = self.lat.edges[e]
        self.apply_pauli(E.u, E.label.upper())
        return dest

    # ----------------------------------------------------------------- reading
    def values(self) -> dict:
        peek = self.sim.peek_observable_expectation
        return {
            "plaquettes": [peek(p) for p in self.plaq],
            "pairs": {pid: peek(P.op) for pid, P in self.pairs.items()},
        }

    def colouring(self) -> tuple[list[int], bool]:
        """Bicolour faces (0 = white/e, 1 = black/m). Faces sharing a link used
        by an odd number of paths get the same colour, otherwise different.
        Returns (colours, consistent); inconsistent means an even-link loop
        winds around the torus, so e and m are only defined locally."""
        count = [0] * len(self.lat.edges)
        for P in self.pairs.values():
            for e in P.path:
                count[e] += 1
        nb = self.lat.face_neighbours()
        col = [-1] * len(self.lat.faces)
        ok = True
        for start in range(len(col)):
            if col[start] >= 0:
                continue
            col[start] = 0
            stack = [start]
            while stack:
                f = stack.pop()
                for g, e in nb[f]:
                    want = col[f] if count[e] % 2 else 1 - col[f]
                    if col[g] < 0:
                        col[g] = want
                        stack.append(g)
                    elif col[g] != want:
                        ok = False
        return col, ok

    def summary(self) -> dict:
        v = self.values()
        col, ok = self.colouring()
        return {
            "fermions": sum(1 for x in v["pairs"].values() if x == -1),
            "e": sum(1 for f, x in enumerate(v["plaquettes"]) if x == -1 and col[f] == 0),
            "m": sum(1 for f, x in enumerate(v["plaquettes"]) if x == -1 and col[f] == 1),
            "majoranas": self.majoranas(),
            "colouring_consistent": ok,
        }

    def record(self) -> dict:
        """The run as JSON-friendly data; MatchingCode.replay rebuilds it."""
        return {"lattice": self.lat.name, "seed": self.seed,
                "plaquette_raw": list(self.plaquette_raw), "log": list(self.log)}

    def _record(self, action, **kw):
        self.log.append({"action": action, **kw})


def _simplify(path: list[int]) -> list[int]:
    """Cancel back-and-forth steps (an edge used twice in a row)."""
    out = []
    for e in path:
        if out and out[-1] == e:
            out.pop()
        else:
            out.append(e)
    return out
