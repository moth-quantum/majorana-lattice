"""Test to add to majorana-lattice/tests (with docs/patches/majorana-lattice-fix-collisions.patch): fix_value on a measurement
between two DIFFERENT Majorana pairs (a collision). The fix uses the string of one of the pairs, so the new dimer and the new
(a, b) pair become deterministic for fix_value +-1 and keep P(+1) = (1 + v) / 2 otherwise, and every pair value the
bookkeeping reports matches a real readout of that pair's operator."""
import random

import numpy as np
import pytest
from majorana_lattice import majorana_moves, majoranas, measure_link, move_majorana, prepare, readout, release, run, star


def two_majorana_pairs(lat, seed):
    rng = random.Random(seed)
    p = prepare(lat)
    d1, d2 = rng.sample([P for P in p.pairs if P.is_dimer], 2)
    p = release(release(p, d1.a), d2.a)
    for _ in range(rng.randrange(0, 5)):          # longer strings: the fixer may be a long path
        v = rng.choice(majoranas(p))
        moves = majorana_moves(p, v)
        if moves:
            p = move_majorana(p, v, rng.choice(moves)[0], 1.0)
    return p


def collisions(p):
    return [e for e, E in enumerate(p.lattice.edges)
            if p.pair_at(E.u) is not p.pair_at(E.v) and not p.pair_at(E.u).is_dimer and not p.pair_at(E.v).is_dimer]


@pytest.mark.parametrize("v", [1.0, -1.0, 0.0, 0.7])
def test_fix_value_on_a_collision(v):
    lat = star(4, 4)
    checked = 0
    for seed in range(60):
        p = two_majorana_pairs(lat, seed)
        cs = collisions(p)
        if not cs:
            continue
        e = random.Random(seed).choice(cs)
        q = measure_link(p, e, v)
        for k, P in enumerate(q.pairs):
            q = readout(q, f"p{k}", P)
        res = run(q, shots=400, seed=1)
        for k, P in enumerate(q.pairs):
            val = res.value(P)
            if val is not None:
                assert np.array_equal(val, res[f"p{k}"]), f"seed {seed}: bookkeeping disagrees with a real readout"
        plus = float((res.value(q.pair_at(lat.edges[e].u)) == 1).mean())
        assert abs(plus - (1 + v) / 2) < 0.08
        checked += 1
    assert checked >= 5, checked
