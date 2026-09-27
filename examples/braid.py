"""Exchange two Majoranas twice and watch both pair parities flip.

Two dimers are split into Majorana pairs A = (a1, a2) and B = (b1, b2). The
outer Majoranas a1 and b2 are moved well away, then a2 and b1 are exchanged
inside a small disc, so the exchange loops around nothing else. The paper's
exchange moves one Majorana to a temporary spot, moves the other into its
place, then moves the first into the other's old place. One exchange puts
both parities into superposition; two flip them: R^2 |kA,kB> = |-kA,-kB>.

    python examples/braid.py
"""
from majorana_lattice import MatchingCode, star

L = star(10, 10)
W, H = L.period
CX, CY = W / 2, H / 2


def dist2(q, x, y):
    px, py = L.positions[q]
    return (px - x) ** 2 + (py - y) ** 2


def nearest(x, y, exclude=()):
    return min((q for q in range(L.n) if q not in exclude), key=lambda q: dist2(q, x, y))


def go_near(c, p, x, y, avoid):
    """Move the Majorana at p to the closest reachable vertex near (x, y)."""
    for t in sorted(range(L.n), key=lambda q: dist2(q, x, y))[:40]:
        if t in avoid:
            continue
        try:
            return c.move_majorana(p, t, avoid=avoid)[-1]
        except ValueError:
            continue
    raise RuntimeError("no reachable vertex near the target")


c = MatchingCode(L, seed=3)
peek = c.sim.peek_observable_expectation
label = {1: "no fermion", -1: "fermion", 0: "superposition"}

# two neighbouring dimers near the centre, split into Majorana pairs
A = c.pairs[c.owner[nearest(CX - 0.6, CY)]]
B = c.pairs[c.owner[nearest(CX + 0.6, CY, exclude={A.a, A.b})]]
a1, a2 = sorted((A.a, A.b), key=lambda q: L.positions[q][0])   # a2 is the inner end
b1, b2 = sorted((B.a, B.b), key=lambda q: L.positions[q][0])   # b1 is the inner end
c.release(a1)
c.release(b1)

# move the outer ends far out to the left and right
far = {q for q in range(L.n) if dist2(q, CX, CY) > 5.0 ** 2}
a1 = go_near(c, a1, CX - 4.0, CY, {a2, b1, b2} | far)
b2 = go_near(c, b2, CX + 4.0, CY, {a1, a2, b1} | far)

# parity operators of the two pairs, as they stand now
piA = c.pairs[c.owner[a1]].op.copy()
piB = c.pairs[c.owner[b2]].op.copy()

# the exchange stays within a small disc round the two inner Majoranas
mx, my = [(L.positions[a2][k] + L.positions[b1][k]) / 2 for k in (0, 1)]
outside = {q for q in range(L.n) if dist2(q, mx, my) > 3.0 ** 2}


def exchange(p, q):
    others = {a1, b2} | outside
    t = next(t for _, t in c.hop_options(p) if t not in others | {q})
    t = c.move_majorana(p, t, avoid=others | {q})[-1]
    c.move_majorana(q, p, avoid=others | {t})
    c.move_majorana(t, q, avoid=others | {p})


print("start:          pair A", label[peek(piA)], "| pair B", label[peek(piB)])
exchange(a2, b1)
print("one exchange:   pair A", label[peek(piA)], "| pair B", label[peek(piB)])
exchange(a2, b1)
print("two exchanges:  pair A", label[peek(piA)], "| pair B", label[peek(piB)])
