"""Saving a Program without its vacuum prep, and the prep itself: unchanged
circuits, and fast enough to rebuild on load."""
import hashlib
import time

import numpy as np
import pytest
import stim

from majorana_lattice import (PRESETS, Program, dimers, from_instructions, majorana_moves, majoranas,
                              move_majorana, plaquette_values, prepare, readout, run, split_dimer, star)

# sha256 of the prep circuit text and the signed plaquettes, recorded before
# prepare() was made faster. Any change to the prep circuit shows up here.
GOLDEN = {
    "star 4x4 x": "71a23cf9ce2d33fa6366eadee5c84a4a48b828dcd448af569d68bafa5a191774",
    "star 4x4 y": "30e49783870b5a9370920c7051a665bbc40ec88eaa1017f38bb241a69673f7fc",
    "star 4x4 z": "01e9b280b1f6d42e9108db3c278b3ed046e2c2df49eb020b65cf5e0ef96aa1a3",
    "star 6x6 x": "4a64cb57d67d74122f1ab9ab72b0a6368d981cfd55e3720af0ca31611e483b6d",
    "star 6x6 y": "84dd3d5230238ad94d9a32e173cb618159fa364bc4a89db075e56f50bde58083",
    "star 6x6 z": "369f31dddf120f42e475a5142a0b21613190293e0e1d2c09e513f67fec955490",
    "honeycomb 4x4 x": "9523f2e695d0ec2c19c839eb17a4741422cfbe61076deaa585cdc8bb737b11ac",
    "honeycomb 4x4 y": "67746242231e0d7e5184e07789814e5d9704228edd5d5a528164aabd9e040ca2",
    "honeycomb 4x4 z": "24d4e1f12d936cef14ccdcbf66f50da874a74dfdfec019b85ecae6e3af179ffe",
    "modified_honeycomb 4x4 x": "416eb7282483d7ba0cb7b1731c0aa901075f16aee4b8733dad719fdb96184621",
    "modified_honeycomb 4x4 y": "4a1eab908d3b1aa2e4b36e47ae5ce9edcbdef173e6e98c7740ceeccb28bb9568",
    "modified_honeycomb 4x4 z": "92d87cb7081f41cadd40f2a72e29ba05d27f5ea98e6b2a1724fd3b1cb7453ff5",
}


@pytest.mark.parametrize("key", GOLDEN)
def test_prep_circuit_is_unchanged(key):
    name, size, background = key.split()
    cols, rows = map(int, size.split("x"))
    p = prepare(PRESETS[name](cols, rows), background)
    text = str(p.prep) + "\n#plaq\n" + "\n".join(str(o) for o in p.plaq_ops)
    assert hashlib.sha256(text.encode()).hexdigest() == GOLDEN[key]


def test_prepare_is_fast():
    L = star(12, 12)                       # a fresh lattice: the vacuum is cached on it
    t = time.perf_counter()
    prepare(L)
    assert time.perf_counter() - t < 1.0


def some_program(background="z"):
    """A split, a move with a conditional fix, and a readout."""
    p = prepare(star(4, 4), background)
    p = split_dimer(p, dimers(p)[5])
    v = majoranas(p)[0]
    p = readout(p, "before", p.pair_at(v))
    link, _ = majorana_moves(p, v)[0]
    p = move_majorana(p, v, link, fix_value=1)
    return readout(p, "after", p.pair_at(majoranas(p)[0]))


@pytest.mark.parametrize("background", "xyz")
def test_from_instructions_rebuilds_the_program(background):
    p = some_program(background)
    assert any(inst.name.startswith("C") and inst.tag.startswith("fix:") for inst in p.instructions)
    text = str(p.instructions_circuit())                   # what a store would keep
    q = Program.from_instructions(p.lattice, background, stim.Circuit(text))
    assert str(q.circuit) == str(p.circuit)
    assert q.readouts == p.readouts and q.measurements == p.measurements
    assert majoranas(q) == majoranas(p) and plaquette_values(q) == plaquette_values(p)
    a, b = run(p, shots=200, seed=11), run(q, shots=200, seed=11)
    assert np.array_equal(a.bits, b.bits)
    assert a["before"].tolist() == b["before"].tolist() and a["after"].tolist() == b["after"].tolist()


def test_from_instructions_matches_from_circuit():
    p = some_program()
    full = Program.from_circuit(p.lattice, p.circuit)
    part = from_instructions(p.lattice, "z", p.instructions_circuit())
    assert str(full.circuit) == str(part.circuit)
    assert [(P.a, P.b, P.kind, P.path, str(P.op), P.flip, P.records) for P in full.pairs] == \
        [(P.a, P.b, P.kind, P.path, str(P.op), P.flip, P.records) for P in part.pairs]


def test_instructions_circuit_holds_no_prep():
    p = some_program()
    c = p.instructions_circuit()
    assert not any(inst.tag.startswith("prepare:") for inst in c)
    assert str(p.prep) + "\n" + str(c) == str(p.circuit)


def test_from_instructions_refuses_a_full_circuit():
    p = some_program()
    with pytest.raises(ValueError):
        Program.from_instructions(p.lattice, "z", p.circuit)
