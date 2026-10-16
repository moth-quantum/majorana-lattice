"""Linked Majorana pairs, the experiment that exposed the vacuum convention.

Two pairs are made from vacuum. Linked order: open V, open H, close V,
close H. Unlinked order: open V, close V, open H, close H. Ising anyon theory
says linked pairs can never both come back empty, and unlinked pairs always
do. Moves stay inside a disc so paths never wind round the torus.

Each order is built once as a circuit and run for many shots.

    python examples/linking.py
"""
from majorana_lattice import (Link, Vertex, create_majoranas, fuse, partner, prepare, route_majorana,
                              run, star, to_axial, to_xy)

L = star(10, 10)
C, R = (5.0, 5.0), 2.2                              # centre (axial) and radius
VERTICES = [L.address(Vertex, q) for q in range(L.n)]
OUTSIDE = {v for v in VERTICES if L.distance(v, C) > 3.6}


def at(dx, dy):
    """The point (dx, dy) away from the centre, in Cartesian units."""
    x, y = to_xy(*C)
    return to_axial(x + dx, y + dy)


def ends(link):
    E = L.edges[L.index(link)]
    return L.address(Vertex, E.u), L.address(Vertex, E.v)


def create_pair(p, uv):
    """Measure the link between two dimers nearest uv (fixing a -1); returns
    the program and the two new Majoranas."""
    def joins_two_dimers(link):
        u, v = ends(link)
        return partner(p, u) != v and all(p.pair_at(w).is_dimer for w in (u, v))
    link = min((L.address(Link, e) for e in range(len(L.edges))
                if joins_two_dimers(L.address(Link, e))), key=lambda l: L.distance(l, uv))
    u, v = ends(link)
    new = [partner(p, u), partner(p, v)]
    return create_majoranas(p, link), new


def move(p, v, uv, avoid):
    """Move the Majorana at v to a reachable vertex near uv, inside the disc."""
    blocked = set(avoid) | OUTSIDE
    for extra in range(8):
        x, y = to_xy(*uv)
        t = min((w for w in VERTICES if w not in blocked),
                key=lambda w: L.distance(w, to_axial(x + 0.3 * extra, y)))
        try:
            return route_majorana(p, v, t, avoid=blocked), t
        except ValueError:
            blocked |= {t}
    raise RuntimeError("no route")


def close(p, pair, avoid, name):
    """Bring the two Majoranas together and fuse them; the outcome is `name`."""
    a, b = pair
    for w in L.incident[L.index(b)]:
        t = L.address(Vertex, L.edges[w].other(L.index(b)))
        if t in avoid or t in OUTSIDE:
            continue
        try:
            if a != t:
                p = route_majorana(p, a, t, avoid=set(avoid) | {b} | OUTSIDE)
                a = t
            break
        except ValueError:
            continue
    link = L.address(Link, L.edge_between(L.index(a), L.index(b)))
    return fuse(p, link, name)


def experiment(linked):
    p = prepare(L)
    p, V = create_pair(p, C)
    p, v0 = move(p, V[0], at(0, R), V)
    p, v1 = move(p, V[1], at(0, -R), [v0, V[1]])
    V = [v0, v1]
    if linked:
        p, H = create_pair(p, C)
        p, h0 = move(p, H[0], at(R, 0), V + H)
        p, h1 = move(p, H[1], at(-R, 0), V + [h0, H[1]])
        p = close(p, V, [h0, h1], "V")
        p = close(p, [h0, h1], [], "H")
    else:
        p = close(p, V, [], "V")
        p, H = create_pair(p, C)
        p, h0 = move(p, H[0], at(R, 0), H)
        p, h1 = move(p, H[1], at(-R, 0), [h0, H[1]])
        p = close(p, [h0, h1], [], "H")
    return p


if __name__ == "__main__":
    shots = 1000
    for linked in (False, True):
        res = run(experiment(linked), shots=shots, seed=1)
        v, h = res["V"], res["H"]
        name = "linked  " if linked else "unlinked"
        print(f"{name}: both empty {((v == 1) & (h == 1)).sum()}/{shots}, "
              f"both fermions {((v == -1) & (h == -1)).sum()}/{shots}")
