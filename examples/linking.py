"""Linked Majorana pairs, the experiment that exposed the vacuum convention.

Two pairs are made from vacuum. Linked order: open V, open H, close V,
close H. Unlinked order: open V, close V, open H, close H. Ising anyon theory
says linked pairs can never both come back empty, and unlinked pairs always
do. Moves stay inside a disc so paths never wind round the torus.

    python examples/linking.py
"""
from majorana_lattice import MatchingCode, link_pauli, star

L = star(10, 10)
W, H = L.period
CX, CY, R = W / 2, H / 2, 2.2
OUTSIDE = {q for q, (x, y) in enumerate(L.positions) if (x - CX) ** 2 + (y - CY) ** 2 > 3.6 ** 2}


def nearest(x, y, exclude):
    return min((q for q in range(L.n) if q not in exclude),
               key=lambda q: (L.positions[q][0] - x) ** 2 + (L.positions[q][1] - y) ** 2)


def create_pair(c, x, y):
    """Measure the link between two dimers nearest (x, y); fix a -1."""
    e = min((e for e, E in enumerate(L.edges)
             if c.owner[E.u] != c.owner[E.v]
             and all(c.is_dimer(c.pairs[c.owner[q]]) for q in (E.u, E.v))),
            key=lambda e: (L.positions[L.edges[e].u][0] - x) ** 2 + (L.positions[L.edges[e].u][1] - y) ** 2)
    E = L.edges[e]
    A, B = c.pairs[c.owner[E.u]], c.pairs[c.owner[E.v]]
    ends = [A.other(E.u), B.other(E.v)]
    old = A.path[0]
    if c.measure_link(e) == -1:
        c.apply_link(old)
    return ends


def move(c, p, x, y, avoid):
    blocked = set(avoid) | OUTSIDE
    for extra in range(8):
        t = nearest(x + 0.3 * extra, y, blocked)
        try:
            return c.route_majorana(p, t, avoid=blocked)[-1]
        except ValueError:
            blocked = blocked | {t}
    raise RuntimeError("no route")


def fuse(c, pair, avoid):
    """Bring the two Majoranas together and measure the link between them."""
    p, q = pair
    for f in L.incident[q]:
        t = L.edges[f].other(q)
        if t in avoid or t in OUTSIDE:
            continue
        try:
            if p != t:
                p = c.route_majorana(p, t, avoid=set(avoid) | {q} | OUTSIDE)[-1]
            break
        except ValueError:
            continue
    e = L.edge_between(p, q)
    return c.measure_link(e)


def run(seed, linked):
    c = MatchingCode(L, seed=seed)
    V = create_pair(c, CX, CY)
    V = [move(c, V[0], CX, CY + R, V), move(c, V[1], CX, CY - R, V)]
    if linked:
        Hp = create_pair(c, CX, CY)
        Hp = [move(c, Hp[0], CX + R, CY, V + Hp), move(c, Hp[1], CX - R, CY, V + Hp)]
        v, h = fuse(c, V, Hp), fuse(c, Hp, [])
    else:
        v = fuse(c, V, [])
        Hp = create_pair(c, CX, CY)
        Hp = [move(c, Hp[0], CX + R, CY, Hp), move(c, Hp[1], CX - R, CY, Hp)]
        h = fuse(c, Hp, [])
    return v, h


if __name__ == "__main__":
    for linked in (False, True):
        results = [run(s, linked) for s in range(12)]
        empty = sum(v == 1 and h == 1 for v, h in results)
        both = sum(v == -1 and h == -1 for v, h in results)
        name = "linked  " if linked else "unlinked"
        print(f"{name}: both empty {empty}/12, both fermions {both}/12")
