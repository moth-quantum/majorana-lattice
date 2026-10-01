"""Checks of the Stim engine against the physics in arXiv:1501.07779."""
import random

import stim

from majorana_lattice import PRESETS, MatchingCode, link_pauli, star


def check(name, cond):
    assert cond, name


def test_lattices():
    for name, fn in PRESETS.items():
        L = fn(4, 4)
        comm = all(link_pauli(L, e).commutes(_plaq(L, f))
                   for e in range(len(L.edges)) for f in range(len(L.faces)))
        check(f"{name}: every plaquette commutes with every link", comm)


def _plaq(L, f):
    p = stim.PauliString(L.n)
    for q, P in L.faces[f].pauli.items():
        p[q] = P
    return p


def test_vacuum():
    c = MatchingCode(star(4, 4), seed=1)
    v = c.values()
    check("start: all dimers empty", all(x == 1 for x in v["pairs"].values()))
    check("start: all plaquettes empty", all(x == 1 for x in v["plaquettes"]))
    check("start: 48 blue (z) dimers on 96 qubits", len(c.pairs) == 48)
    check("start: colouring consistent", c.colouring()[1])


def test_link_clicks():
    rng = random.Random(7)
    c = MatchingCode(star(4, 4), seed=2)
    ok_parity, ok_plaq, ok_det = True, True, True
    for _ in range(300):
        e = rng.randrange(len(c.lat.edges))
        c.measure_link(e)
        v = c.values()
        # link measurements never disturb plaquettes
        ok_plaq &= all(x == 1 for x in v["plaquettes"])
        # every tracked pair has a definite value and fermions come in pairs
        ok_det &= all(x != 0 for x in v["pairs"].values())
        ok_parity &= sum(1 for x in v["pairs"].values() if x == -1) % 2 == 0
    check("link clicks never populate plaquettes", ok_plaq)
    check("link clicks leave every pair definite", ok_det)
    check("fermions are created in pairs", ok_parity)


def test_pauli_errors():
    c = MatchingCode(star(4, 4), seed=4)
    c.apply_pauli(0, "X")
    v, s = c.values(), c.summary()
    flipped = sum(1 for x in v["plaquettes"] if x == -1)
    check("one X error flips two plaquettes", flipped == 2)
    # The error's label decides whether the anyons match (e,e or m,m + even
    # fermions) or differ (e,m + odd fermions), as the paper shows.
    ok = (s["e"] + s["m"] == 2) and ((s["e"] == s["m"]) == (s["fermions"] % 2 == 1))
    check("e x m = fermion: mixed anyon pair comes with an odd number of fermions", ok)


def test_floquet_repeat():
    c = MatchingCode(star(4, 4), seed=5)
    first = c.measure_label("x")
    second = c.measure_label("x")
    check("a Floquet round repeated gives the same outcomes", first == second)
    c.measure_label("y")
    check("after x then y rounds, pairs are still all definite",
          all(x != 0 for x in c.values()["pairs"].values()))


def _setup_two_pairs(seed):
    c = MatchingCode(star(6, 6), seed=seed)
    L = c.lat
    # pick two z dimers a few steps apart and release them
    pids = sorted(c.pairs)
    A = c.pairs[pids[0]]
    # find another dimer at graph distance ~4 from A
    from collections import deque
    dist = {A.b: 0}
    dq = deque([A.b])
    while dq:
        q = dq.popleft()
        for e in L.incident[q]:
            w = L.edges[e].other(q)
            if w not in dist:
                dist[w] = dist[q] + 1
                dq.append(w)
    b1 = next(q for q in range(L.n) if dist.get(q) == 5 and c.owner[q] != c.owner[A.a])
    B = c.pairs[c.owner[b1]]
    a1, a2 = A.a, A.b
    b2 = B.other(b1)
    piA = A.op.copy()
    piB = B.op.copy()
    c.release(a1)
    c.release(b1)
    return c, (a1, a2, b1, b2), piA, piB


def _exchange(c, p, q, fixed):
    """Swap the Majoranas at p and q using a temporary spot (the paper's
    three-leg teleportation exchange)."""
    free = set(fixed) - {p, q}
    # choose a temporary vertex reachable from p
    for e, t in c.majorana_moves(p):
        if t not in fixed:
            break
    t = c.route_majorana(p, t, avoid=free | {q})[-1]
    c.route_majorana(q, p, avoid=free | {t})
    c.route_majorana(t, q, avoid=free | {p})


def test_braiding():
    flips, singles = 0, 0
    trials = 20
    for s in range(trials):
        c, (a1, a2, b1, b2), piA, piB = _setup_two_pairs(100 + s)
        peek = c.sim.peek_observable_expectation
        assert peek(piA) == 1 and peek(piB) == 1
        fixed = {a1, a2, b1, b2}
        _exchange(c, a2, b1, fixed)
        singles += (peek(piA) == 0 and peek(piB) == 0 and peek(piA * piB) != 0)
        _exchange(c, a2, b1, fixed)
        flips += (peek(piA) == -1 and peek(piB) == -1)
        assert sorted(c.majoranas()) == sorted(fixed), (c.majoranas(), fixed)
    check(f"single exchange puts both parities in superposition ({singles}/{trials})", singles == trials)
    check(f"double exchange flips both parities, R^2 |kA,kB> = |-kA,-kB> ({flips}/{trials})", flips == trials)


def test_route_independence():
    """A pair made from vacuum, wandered around and fused must come back empty,
    whatever route it took (checked on a patch, away from torus wraparound)."""
    import math
    from collections import deque
    L = star(10, 10)
    W, H = L.period
    out = {q for q, (x, y) in enumerate(L.positions) if (x - W / 2) ** 2 + (y - H / 2) ** 2 > 3.6 ** 2}
    empty = done = 0
    for s in range(15):
        c = MatchingCode(L, seed=200 + s)
        r = random.Random(s)
        e = next(e for e, E in enumerate(L.edges)
                 if E.u not in out and E.v not in out and c.owner[E.u] != c.owner[E.v])
        E = L.edges[e]
        A, B = c.pairs[c.owner[E.u]], c.pairs[c.owner[E.v]]
        p, q = A.other(E.u), B.other(E.v)
        if c.measure_link(e) == -1:
            c.apply_link(A.path[0])
        for _ in range(60):
            mover = r.choice([0, 1])
            x = (p, q)[mover]
            opts = [(f, t) for f, t in c.majorana_moves(x) if t not in out and t not in (p, q)
                    and L.edges[f].other(x) not in (p, q)]
            if opts:
                f, t = r.choice(opts)
                new = c.move_majorana(x, f)
                p, q = (new, q) if mover == 0 else (p, new)
        for f in L.incident[q]:
            target = L.edges[f].other(q)
            if target in out:
                continue
            try:
                if p != target:
                    p = c.route_majorana(p, target, avoid=out | {q})[-1]
                break
            except ValueError:
                continue
        link = L.edge_between(p, q)
        if link is None:
            continue  # walled in by the patch edge; skip this run
        done += 1
        empty += c.measure_link(link) == 1
    check(f"a pair from vacuum fuses back to vacuum whatever its route ({empty}/{done})", empty == done and done >= 10)




def _load_example(name):
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parent.parent / "examples" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_linked_pairs():
    """Two pairs from vacuum. Linked (open V, open H, close V, close H): Ising
    theory says both fuse to a fermion. Unlinked: both fuse to vacuum."""
    ex = _load_example("linking")
    linked = [ex.run(s, True) for s in range(6)]
    unlinked = [ex.run(s, False) for s in range(6)]
    check("linked pairs both fuse to a fermion", all(r == (-1, -1) for r in linked))
    check("unlinked pairs both fuse to vacuum", all(r == (1, 1) for r in unlinked))
