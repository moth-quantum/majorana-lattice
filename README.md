# majorana-lattice

Matching codes, anyons and Majorana modes on trivalent lattices, written as
[Stim](https://github.com/quantumlib/Stim) circuits.

This package implements the matching codes of J. R. Wootton,
[arXiv:1501.07779](https://arxiv.org/abs/1501.07779). Every vertex of a
lattice holds a qubit and a Majorana mode. Vertices are paired into dimers,
and you act on the code by measuring two-qubit link operators. That one
action creates, moves and fuses Majoranas, and the braiding comes out of the
real quantum state rather than from a model of it.

The default lattice is the **star lattice** (every honeycomb vertex replaced
by a triangle). Majoranas can reach every vertex using only two-qubit
measurements, and every vertex is equivalent. The engine works on any
lattice where every vertex has three links labelled x, y and z with no two
matching labels meeting; honeycomb and the paper's modified honeycomb are
included as presets.

## What's here

| Path | What it is |
| --- | --- |
| `src/majorana_lattice/` | The Python package: lattices, programs, operations, queries and `run` |
| `examples/` | Runnable experiments: the braid, and linked Majorana pairs |
| `tests/` | Physics tests, and checks of the bookkeeping against Stim |
| `docs/conventions.md` | How to read the board, and the conventions that matter |

## Install

```bash
pip install git+https://github.com/moth-quantum/majorana-lattice
```

or, from a clone, `pip install -e ".[test]"`.

## Quick start

An experiment is built as a `Program` (a Stim circuit plus the bookkeeping
needed to decide what to do next) and nothing runs until `run`. Every
operation takes a Program and returns a new one, so you can branch from any
point. Elements are addressed by unit cell, in axial coordinates on the
lattice vectors, and a name within the cell (see `docs/conventions.md`).

```python
from majorana_lattice import (Link, Vertex, majorana_moves, majoranas, move_majorana,
                              prepare, route_majorana, run, split_dimer, star)

p = prepare(star(4, 4))                   # vacuum: every dimer and plaquette +1
p = split_dimer(p, Link(0, 1, "z"))       # two Majoranas, no gates applied
a, b = majoranas(p)                       # Vertex(0, 1, "D.z"), Vertex(0, 1, "U.z")

link, dest = majorana_moves(p, a)[0]      # moves on offer: [(link, destination)]
p = move_majorana(p, a, link)             # measures the link and fixes it: fix_value=1
p = route_majorana(p, dest, Vertex(2, 3, "D.x"))

p.circuit                                 # a plain stim.Circuit
res = run(p, shots=1000)                  # the only step that executes anything
res.pairs()                               # every pair's value, per shot
res.summary()                             # fermions, e and m, per shot
```

### Fixing a move

`fix_value`, from -1 to 1, is the expected value of the dimer a move leaves
behind: `move_majorana(p, a, link, fix_value=-0.8)` leaves a fermion in 90%
of shots. 1 (the default) always empties it, -1 always leaves a fermion, and
0 gives 50/50. `create_majoranas` and `route_majorana` take it too.

Every fix is built the same way, with one extra qubit after the lattice's,
the fix qubit. It is rotated by RY(θ) with sin²(θ/2) = (1 - fix_value) / 2 and
measured: reading 0, the move's dimer is corrected to +1; reading 1, to -1.
Stim has no RY(θ), so the circuit writes the fix qubit as a heralded channel
that gives the same record (see `docs/conventions.md`).

`measure_link(p, link)` without a `fix_value` keeps the raw outcome and uses
no fix qubit. That is deliberate, and it is not the same as `fix_value=0`:

- a fusion's outcome is the result (in `examples/linking.py` it is -1 in every
  shot), and two Majoranas have no dimer to correct with anyway;
- a Floquet round repeated with `measure_label` repeats its outcomes;
- a link whose ends are already paired has a determined outcome.

Coupling a 50/50 fix qubit would randomise all of these. Pass a `fix_value`
to `measure_link` only when it rewires a dimer and you want the fix.

### Reading results

A pair's value is never measured to report it. The bookkeeping knows which
measurement records it is the parity of, so `run` computes it per shot. To
measure something else, such as a pair's parity from an earlier stage, add
`readout(p, "name", operator)` and read `res["name"]`. `fuse(p, link, "name")`
names a fusion outcome the same way. `res.summary(at=earlier)` gives values
as they were at an earlier stage of the same program.

`Program.from_circuit(lattice, circuit)` rebuilds a Program from its circuit,
so saving the circuit text saves the whole run. The vacuum prep depends only
on the lattice and the background, so a store can keep just
`p.instructions_circuit()` and rebuild with
`Program.from_instructions(lattice, background, circuit)`.

### Executors

`run` takes an executor: any function from a `stim.Circuit` and a shot count
to a bool array of measurement records, one row per shot. Local Stim is the
default.

```python
res = run(p, shots=1000, executor=my_executor)
```

## Tests

```bash
pytest
```

The physics tests check the paper's results on the star lattice, including:

- link measurements never create e or m anyons;
- fermions are created in pairs;
- a single-qubit error gives e/e or m/m with an even number of fermions, or
  e/m with an odd number;
- exchanging two Majoranas twice flips both pair parities;
- a pair made from vacuum fuses back to vacuum whatever route it took;
- two pairs from vacuum whose worldlines are linked both fuse to a fermion,
  and unlinked pairs both fuse back to vacuum, as Ising anyon theory predicts.

The vacuum is checked against the postselected state on every preset and
background, and every pair value the bookkeeping reports is confirmed by
Stim's own flow checker (`stim.Circuit.has_flow`) and by real readouts.

## Citing

If you use this in research, please cite the papers it is built on; see
`CITATION.cff`.

- J. R. Wootton, *A family of stabilizer codes for D(Z₂) anyons and Majorana
  modes*, [arXiv:1501.07779](https://arxiv.org/abs/1501.07779).
- J. R. Wootton, *Measurements of Floquet code plaquette stabilizers*,
  [arXiv:2210.13154](https://arxiv.org/abs/2210.13154).

## Licence

Apache 2.0. See `LICENSE`.
