"""Operations on a Program. Each takes a Program and returns a new one.

Targets are addresses (Vertex, Link, Plaquette) or plain indices. Nothing is
executed: each operation appends Stim instructions and updates the
bookkeeping. Instructions carry tags (MPP[link:12], X[pauli], ...) so that
`from_circuit` can rebuild the Program from its circuit alone.

  measure_link(p, link)        measure a link operator; rewires the pairing
  apply_link(p, link)          apply a link operator as a unitary
  apply_pauli(p, v, P)         apply X, Y or Z to one qubit
  release(p, v)                split the dimer holding v into two Majoranas,
                               or join a released pair back (bookkeeping only)
  split_dimer(p, link)         release, addressed by the dimer's link
  measure_label(p, L)          measure every link with label L
  create_majoranas(p, link)    measure a link between two dimers, fixing -1
  fuse(p, link, name)          measure the link between two Majoranas
  move_majorana(p, v, link)    move a Majorana across a link, fixing -1
  route_majorana(p, a, b)      walk a Majorana from a to b by moves
  move_fermion(p, v, link)     move the fermion on v's pair across a link
  move_anyon(p, f, link)       move the anyon on plaquette f across a link
  create_anyons(p, link)       flip the two plaquettes either side of a link
  readout(p, name, op)         measure any operator, for `run` to report
"""
from __future__ import annotations

import random

import stim

from . import queries as Q
from .lattice import Lattice, Link, Plaquette, Vertex
from .program import Pair, Program, _mpp_targets, link_pauli, prepare, resolve, simplify_path

_CTRL = {"X": "CX", "Y": "CY", "Z": "CZ"}


# ------------------------------------------------------------- primitives
def measure_link(p: Program, link, fix: bool | float = False, rng: random.Random | None = None,
                 name: str | None = None) -> Program:
    """Measure link (i, j). The pairs holding i and j are rewired: (i, a) and
    (j, b) become the dimer (i, j) and the pair (a, b). With `fix`, a -1
    outcome is undone by applying the links of the old dimer on one side,
    conditioned on the outcome (a float is the chance of doing so, drawn
    from `rng` now, not at run time). `name` makes the outcome a readout."""
    if name is not None and name in p.readouts:
        raise ValueError(f"readout {name!r} already exists")
    lat = p.lattice
    e = resolve(lat, link, Link)
    E = lat.edges[e]
    i, j = E.u, E.v
    K = link_pauli(lat, e)
    do_fix = _decide(fix, rng)
    k = p.num_records
    pairs, owner = list(p.pairs), list(p.owner)
    pi, pj = owner[i], owner[j]
    fixer = None
    if pi == pj:
        # i and j already paired: the link commutes with everything
        pairs[pi] = Pair(i, j, (e,), K, "string", False, frozenset({k}))
    else:
        A, B = pairs[pi], pairs[pj]
        a, b = A.other(i), B.other(j)
        new_op = A.op * B.op * K          # commutes with K; real sign
        assert new_op.sign in (1, -1)
        kind = "string" if (A.kind == "string" and B.kind == "string") else "majorana"
        path = tuple(simplify_path(list(A.path_from(a)) + [e] + list(B.path_from(j))))
        known = A.records is not None and B.records is not None
        pairs[pi] = Pair(i, j, (e,), K, "string", False, frozenset({k}))
        pairs[pj] = Pair(a, b, path, new_op, kind, A.flip ^ B.flip,
                         A.records ^ B.records ^ {k} if known else None)
        owner[i] = owner[j] = pi
        owner[a] = owner[b] = pj
        fixer = B if B.is_dimer else A if A.is_dimer else None
    if do_fix and fixer is None:
        raise ValueError(f"no dimer next to link {e} to fix a -1 outcome with")
    tag = f"link:{e}" if name is None else f"link:{e}:{name}"
    readouts = p.readouts if name is None else {**p.readouts, name: k}
    q = p._evolve(instructions=p.instructions + (_mpp(K, tag),), pairs=tuple(pairs),
                  owner=tuple(owner), measurements=p.measurements + (("link", e),),
                  readouts=readouts)
    if do_fix:
        for f in fixer.path:
            q = _conditional_link(q, f, k)
    return q


def apply_link(p: Program, link) -> Program:
    """Apply a link operator as a unitary: toggles a fermion on each of the
    two pairs it joins."""
    lat = p.lattice
    e = resolve(lat, link, Link)
    E = lat.edges[e]
    inst = stim.CircuitInstruction(E.label.upper(), [E.u, E.v], tag=f"link:{e}")
    return _flip_by(p, link_pauli(lat, e), inst)


def apply_pauli(p: Program, v, P: str) -> Program:
    """Apply X, Y or Z to one qubit."""
    q = resolve(p.lattice, v, Vertex)
    P = P.upper()
    op = stim.PauliString(p.lattice.n)
    op[q] = P
    return _flip_by(p, op, stim.CircuitInstruction(P, [q], tag="pauli"))


def release(p: Program, v) -> Program:
    """Stop treating the pair holding v as a stabilizer: its two ends become
    Majoranas (or join back up if already released). No gates; a tagged I
    marks it in the circuit."""
    q = resolve(p.lattice, v, Vertex)
    pairs = list(p.pairs)
    P = pairs[p.owner[q]]
    pairs[p.owner[q]] = Pair(P.a, P.b, P.path, P.op, "majorana" if P.kind == "string" else "string",
                             P.flip, P.records)
    inst = stim.CircuitInstruction("I", [q], tag=f"release:{q}")
    return p._evolve(instructions=p.instructions + (inst,), pairs=tuple(pairs))


def readout(p: Program, name: str, target) -> Program:
    """Measure an operator so `run` reports it under `name`. The target is a
    stim.PauliString, a Pair (e.g. from an earlier Program) or a Plaquette.
    Tracked pairs and plaquettes it anticommutes with become undetermined."""
    if name in p.readouts:
        raise ValueError(f"readout {name!r} already exists")
    if isinstance(target, Pair):
        op = target.op
    elif isinstance(target, (Plaquette, int)) and not isinstance(target, bool):
        op = p.plaq_ops[resolve(p.lattice, target, Plaquette)]
    elif isinstance(target, stim.PauliString):
        op = target
    else:
        raise TypeError(f"cannot read out {target!r}")
    pairs = tuple(P if P.op.commutes(op) else Pair(P.a, P.b, P.path, P.op, P.kind, False, None)
                  for P in p.pairs)
    flips = tuple(f if o.commutes(op) else None for o, f in zip(p.plaq_ops, p.plaq_flips))
    return p._evolve(instructions=p.instructions + (_mpp(op, f"readout:{name}"),), pairs=pairs,
                     plaq_flips=flips, measurements=p.measurements + (("readout", name),),
                     readouts={**p.readouts, name: p.num_records})


# --------------------------------------------------------------- compound
def split_dimer(p: Program, link) -> Program:
    """Release the dimer on this link into two Majoranas."""
    if not Q.is_dimer(p, link):
        raise ValueError(f"link {link} is not a dimer")
    return release(p, p.lattice.edges[resolve(p.lattice, link, Link)].u)


def measure_label(p: Program, label: str) -> Program:
    """Measure every link with this label: one round of a Floquet schedule."""
    for e, E in enumerate(p.lattice.edges):
        if E.label == label:
            p = measure_link(p, e)
    return p


def create_majoranas(p: Program, link, fix: bool | float = True, rng=None) -> Program:
    """Measure a link between two dimers: a new dimer forms and a Majorana
    pair appears at the two far ends. A -1 outcome is fixed by default."""
    lat = p.lattice
    E = lat.edges[resolve(lat, link, Link)]
    A, B = p.pairs[p.owner[E.u]], p.pairs[p.owner[E.v]]
    if p.owner[E.u] == p.owner[E.v] or not (A.is_dimer and B.is_dimer):
        raise ValueError(f"link {link} does not join two dimers")
    return measure_link(p, link, fix=fix, rng=rng)


def fuse(p: Program, link, name: str | None = None) -> Program:
    """Measure the link between two Majoranas; the outcome is their fusion
    channel (+1 vacuum, -1 fermion), reported as `name`."""
    lat = p.lattice
    E = lat.edges[resolve(lat, link, Link)]
    if p.pairs[p.owner[E.u]].is_dimer or p.pairs[p.owner[E.v]].is_dimer:
        raise ValueError(f"link {link} does not join two Majoranas")
    return measure_link(p, link, name=name)


def move_majorana(p: Program, v, link, fix: bool | float = True, rng=None) -> Program:
    """Move the Majorana at v across a link; it lands two vertices along (see
    majorana_moves). A -1 outcome is fixed by default."""
    lat = p.lattice
    e = resolve(lat, link, Link)
    q = resolve(lat, v, Vertex)
    if not any(resolve(lat, ee, Link) == e for ee, _ in Q.majorana_moves(p, q)):
        raise ValueError(f"the Majorana at {v} cannot move across link {link}")
    return measure_link(p, e, fix=fix, rng=rng)


def route_majorana(p: Program, src, dst, avoid=(), fix: bool | float = True, rng=None) -> Program:
    """Walk the Majorana at src to dst by moves, avoiding given vertices."""
    if not Q.is_majorana(p, src):
        raise ValueError(f"vertex {src} is not a Majorana")
    plan = Q.plan_route(p, src, dst, avoid)
    if plan is None:
        raise ValueError(f"no route from {src} to {dst}")
    for q, e in plan:
        p = move_majorana(p, q, e, fix, rng)
    return p


def move_fermion(p: Program, v, link) -> Program:
    """Move the fermion on v's pair across a link (see fermion_moves)."""
    e = resolve(p.lattice, link, Link)
    if not any(resolve(p.lattice, ee, Link) == e for ee, _ in Q.fermion_moves(p, v)):
        raise ValueError(f"the pair holding {v} has no move across link {link}")
    return apply_link(p, e)


def move_anyon(p: Program, f, link) -> Program:
    """Move the e or m anyon on plaquette f across a link (see anyon_moves)."""
    e = resolve(p.lattice, link, Link)
    if not any(resolve(p.lattice, ee, Link) == e for ee, _ in Q.anyon_moves(p, f)):
        raise ValueError(f"the anyon on plaquette {f} cannot move across link {link}")
    return create_anyons(p, e)


def create_anyons(p: Program, link) -> Program:
    """Apply a link's own Pauli to its first endpoint, flipping exactly the two
    plaquettes either side of it."""
    E = p.lattice.edges[resolve(p.lattice, link, Link)]
    return apply_pauli(p, E.u, E.label.upper())


# --------------------------------------------------------------- replay
def from_circuit(lattice: Lattice, circuit: stim.Circuit) -> Program:
    """Rebuild a Program from a circuit made by Program.circuit."""
    insts = list(circuit)
    if not insts or not insts[0].tag.startswith("prepare:"):
        raise ValueError("circuit does not start with a tagged vacuum preparation")
    p = prepare(lattice, insts[0].tag.split(":")[1])
    rest = insts[len([x for x in insts if x.tag.startswith("prepare:")]):]
    for inst in rest:
        kind, _, arg = inst.tag.partition(":")
        groups = inst.target_groups()
        if kind == "link" and inst.name == "MPP":
            e, _, name = arg.partition(":")
            for _ in groups:
                p = measure_link(p, int(e), name=name or None)
        elif kind == "link":
            for _ in range(len(inst.targets_copy()) // 2):     # one-qubit gate on both ends
                p = apply_link(p, int(arg))
        elif kind == "fix":
            p = _conditional_link(p, int(arg), p.num_records + groups[0][0].value)
        elif kind == "pauli":
            for (t,) in groups:
                p = apply_pauli(p, t.value, inst.name)
        elif kind == "release":
            p = release(p, int(arg))
        elif kind == "readout":
            for g in groups:
                op = stim.PauliString(lattice.n)
                for t in g:
                    op[t.value] = t.pauli_type
                if g[0].is_inverted_result_target:
                    op = -op
                p = readout(p, arg, op)
        else:
            raise ValueError(f"unrecognised instruction {inst}")
    return p


# ---------------------------------------------------------------- helpers
def _conditional_link(p: Program, f: int, k: int) -> Program:
    """Apply link f only if record k reads -1: toggles k in every pair the
    link anticommutes with."""
    lat = p.lattice
    F = lat.edges[f]
    rec = stim.target_rec(k - p.num_records)
    inst = stim.CircuitInstruction(_CTRL[F.label.upper()], [rec, F.u, rec, F.v], tag=f"fix:{f}")
    L = link_pauli(lat, f)
    pairs = tuple(P if P.op.commutes(L) or P.records is None
                  else Pair(P.a, P.b, P.path, P.op, P.kind, P.flip, P.records ^ {k})
                  for P in p.pairs)
    assert all(o.commutes(L) for o in p.plaq_ops)
    return p._evolve(instructions=p.instructions + (inst,), pairs=pairs)


def _flip_by(p: Program, op: stim.PauliString, inst) -> Program:
    """Apply a Pauli: flips every pair and plaquette it anticommutes with."""
    pairs = tuple(P if P.op.commutes(op) or P.records is None
                  else Pair(P.a, P.b, P.path, P.op, P.kind, not P.flip, P.records)
                  for P in p.pairs)
    flips = tuple(f if f is None or o.commutes(op) else not f for o, f in zip(p.plaq_ops, p.plaq_flips))
    return p._evolve(instructions=p.instructions + (inst,), pairs=pairs, plaq_flips=flips)


def _mpp(op: stim.PauliString, tag: str) -> stim.CircuitInstruction:
    return stim.CircuitInstruction("MPP", _mpp_targets(op), tag=tag)


def _decide(fix: bool | float, rng) -> bool:
    """Whether to fix: True/False as given, or a coin with chance `fix`."""
    if isinstance(fix, bool):
        return fix
    if fix >= 1 or fix <= 0:
        return fix >= 1
    return (rng or random).random() < fix
