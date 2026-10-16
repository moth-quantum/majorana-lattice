"""Questions about a Program, answered from its bookkeeping alone.

Inputs are addresses (Vertex, Link, Plaquette) or plain indices; outputs are
addresses when the lattice has them, indices otherwise.
"""
from __future__ import annotations

from collections import deque

from .lattice import Link, Plaquette, Vertex
from .program import Program, present, resolve


def majoranas(p: Program) -> list:
    """Vertices not in a dimer: ends of longer strings and of released pairs."""
    qs = sorted(q for P in p.pairs if not P.is_dimer for q in (P.a, P.b))
    return [present(p.lattice, Vertex, q) for q in qs]


def is_majorana(p: Program, v) -> bool:
    return not p.pair_at(v).is_dimer


def dimers(p: Program) -> list:
    """The links that are dimers right now."""
    return [present(p.lattice, Link, P.path[0]) for P in p.pairs if P.is_dimer]


def is_dimer(p: Program, link) -> bool:
    E = p.lattice.edges[resolve(p.lattice, link, Link)]
    return p.owner[E.u] == p.owner[E.v] and p.pairs[p.owner[E.u]].is_dimer


def partner(p: Program, v):
    """The other end of the pair holding v."""
    q = resolve(p.lattice, v, Vertex)
    return present(p.lattice, Vertex, p.pairs[p.owner[q]].other(q))


def plaquette_values(p: Program) -> list[int]:
    """+1 or -1 for every plaquette, or 0 if a readout left it undetermined.
    These never depend on measurement outcomes."""
    return [0 if f is None else (-1 if f else 1) for f in p.plaq_flips]


def majorana_moves(p: Program, v) -> list[tuple]:
    """Moves for the Majorana at v: [(link to measure, destination)]."""
    lat = p.lattice
    q = resolve(lat, v, Vertex)
    out = []
    for e in lat.incident[q]:
        j = lat.edges[e].other(q)
        Pj = p.pairs[p.owner[j]]
        if Pj.is_dimer and p.owner[j] != p.owner[q]:
            out.append((present(lat, Link, e), present(lat, Vertex, Pj.other(j))))
    return out


def fermion_moves(p: Program, v) -> list[tuple]:
    """Moves for the fermion on the pair holding v: [(link, vertex reached)].
    Every link leaving the pair (from either end) to a different pair."""
    lat = p.lattice
    q = resolve(lat, v, Vertex)
    P = p.pairs[p.owner[q]]
    out = []
    for end in (P.a, P.b):
        for e in lat.incident[end]:
            j = lat.edges[e].other(end)
            if p.owner[j] != p.owner[q]:
                out.append((present(lat, Link, e), present(lat, Vertex, j)))
    return out


def anyon_moves(p: Program, f) -> list[tuple]:
    """Moves for the anyon on plaquette f: [(link to cross, destination face)].
    Only moves that keep its type onto an empty plaquette are listed."""
    lat = p.lattice
    k = resolve(lat, f, Plaquette)
    v = plaquette_values(p)
    if v[k] != -1:
        return []
    col, _ = colouring(p)
    return [(present(lat, Link, e), present(lat, Plaquette, g))
            for g, e in lat.face_neighbours()[k] if col[g] == col[k] and v[g] == 1]


def colouring(p: Program) -> tuple[list[int], bool]:
    """Bicolour faces (0 = e, 1 = m). Faces sharing a link used by an odd
    number of pair paths get the same colour, otherwise different. Returns
    (colours, consistent); inconsistent means an even-link loop winds round
    the torus, so e and m are only defined locally."""
    lat = p.lattice
    count = [0] * len(lat.edges)
    for P in p.pairs:
        for e in P.path:
            count[e] += 1
    nb = lat.face_neighbours()
    col = [-1] * len(lat.faces)
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


def plan_route(p: Program, src, dst, avoid=()) -> list[tuple] | None:
    """Shortest route of Majorana moves from src to dst: [(from, link)], or
    None. Every move rewires the vertex it leaves and the vertex it moves
    across, so a route never reuses a vertex, as a stop or as a vertex moved
    across; dimers away from the route stay put, so its moves stay valid."""
    lat = p.lattice
    s, d = resolve(lat, src, Vertex), resolve(lat, dst, Vertex)
    avoid = {resolve(lat, a, Vertex) for a in avoid}
    best = {s: ((), frozenset((s,)))}
    dq = deque([s])
    while dq:
        q = dq.popleft()
        steps, used = best[q]
        if q == d:
            return [(present(lat, Vertex, a), present(lat, Link, e)) for a, e in steps]
        for e in lat.incident[q]:
            j = lat.edges[e].other(q)
            Pj = p.pairs[p.owner[j]]
            if not Pj.is_dimer or p.owner[j] == p.owner[q]:
                continue
            k = Pj.other(j)
            if j in used or k in used or k in avoid or j in avoid or k in best:
                continue
            best[k] = (steps + ((q, e),), used | {j, k})
            dq.append(k)
    return None
