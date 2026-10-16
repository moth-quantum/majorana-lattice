"""Exchange two Majoranas twice and watch both pair parities flip.

Two dimers are split into Majorana pairs A = (a1, a2) and B = (b1, b2). The
outer Majoranas a1 and b2 are moved well away, then a2 and b1 are exchanged
inside a small disc, so the exchange loops around nothing else. The paper's
exchange moves one Majorana to a temporary spot, moves the other into its
place, then moves the first into the other's old place. One exchange puts
both parities into superposition; two flip them: R^2 |kA,kB> = |-kA,-kB>.

The whole experiment is built as a circuit first and run at the end, once
per stage, with the two original parities read out.

    python examples/braid.py
"""
from majorana_lattice import (Vertex, dimers, majorana_moves, prepare, readout,
                              route_majorana, run, split_dimer, star, to_axial, to_xy)

L = star(10, 10)
C = (5.0, 5.0)                                   # the centre, in axial coordinates
VERTICES = [L.address(Vertex, q) for q in range(L.n)]


def shifted(uv, dx):
    """The point dx to the right of uv (Cartesian units)."""
    x, y = to_xy(*uv)
    return to_axial(x + dx, y)


def go_near(p, v, uv, avoid):
    """Move the Majorana at v to the closest reachable vertex near uv."""
    for t in sorted(VERTICES, key=lambda t: L.distance(t, uv))[:40]:
        if t in avoid:
            continue
        try:
            return route_majorana(p, v, t, avoid=avoid), t
        except ValueError:
            continue
    raise RuntimeError("no reachable vertex near the target")


p = prepare(L)

# two neighbouring dimers near the centre, split into Majorana pairs
A = min(dimers(p), key=lambda l: L.distance(l, shifted(C, -0.6)))
B = min((d for d in dimers(p) if d != A), key=lambda l: L.distance(l, shifted(C, 0.6)))
x = lambda v: to_xy(*L.coords(v))[0]
ends = lambda link: sorted((L.address(Vertex, L.edges[L.index(link)].u),
                            L.address(Vertex, L.edges[L.index(link)].v)), key=x)
a1, a2 = ends(A)                                 # a2 is the inner end
b1, b2 = ends(B)                                 # b1 is the inner end
p = split_dimer(split_dimer(p, A), B)

# move the outer ends far out to the left and right
far = {v for v in VERTICES if L.distance(v, C) > 5.0}
p, a1 = go_near(p, a1, shifted(C, -4.0), {a2, b1, b2} | far)
p, b2 = go_near(p, b2, shifted(C, 4.0), {a1, a2, b1} | far)

# the parity operators of the two pairs, as they stand now
piA, piB = p.pair_at(a1), p.pair_at(b2)

# the exchange stays within a small disc round the two inner Majoranas
mid = to_axial(*[(to_xy(*L.coords(a2))[k] + to_xy(*L.coords(b1))[k]) / 2 for k in (0, 1)])
outside = {v for v in VERTICES if L.distance(v, mid) > 3.0}


def exchange(p, u, w):
    others = {a1, b2} | outside
    t = next(t for _, t in majorana_moves(p, u) if t not in others | {w})
    p = route_majorana(p, u, t, avoid=others | {w})
    p = route_majorana(p, w, u, avoid=others | {t})
    return route_majorana(p, t, w, avoid=others | {u})


def parities(stage, shots=2000):
    res = run(readout(readout(stage, "A", piA), "B", piB), shots=shots, seed=1)
    return [f"{(res[k] == 1).mean():.0%} empty" for k in "AB"]


once = exchange(p, a2, b1)
twice = exchange(once, a2, b1)
for name, stage in [("start", p), ("one exchange", once), ("two exchanges", twice)]:
    a, b = parities(stage)
    print(f"{name + ':':15} pair A {a:11} | pair B {b}")
