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
| `src/majorana_lattice/` | The Python package: lattices, the circuit engine, and the interactive engine |
| `js/engine.js` | The interactive engine in plain JavaScript, for browsers and games |
| `explorer/` | An interactive page for exploring the board, and its build script |
| `examples/` | Runnable experiments: the braid, and linked Majorana pairs |
| `tests/` | Physics tests, and checks that the engines agree with each other and with Stim |
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
p = move_majorana(p, a, link)             # measures the link, fixes a -1 by feedforward
p = route_majorana(p, dest, Vertex(2, 3, "D.x"))

p.circuit                                 # a plain stim.Circuit
res = run(p, shots=1000)                  # the only step that executes anything
res.pairs()                               # every pair's value, per shot
res.summary()                             # fermions, e and m, per shot
```

A pair's value is never measured to report it. The bookkeeping knows which
measurement records it is the parity of, so `run` computes it per shot. To
measure something else, such as a pair's parity from an earlier stage, add
`readout(p, "name", operator)` and read `res["name"]`. `fuse(p, link, "name")`
names a fusion outcome the same way. `res.summary(at=earlier)` gives values
as they were at an earlier stage of the same program.

`Program.from_circuit(lattice, circuit)` rebuilds a Program from its circuit,
so saving the circuit text saves the whole run.

### Running elsewhere

`run` takes any executor: a function from a `stim.Circuit` and a shot count to
measurement records. Local Stim is the default. To run through Moth Quantum's
archaeo library ([MothQ_interconnect](https://github.com/moth-quantum/MothQ_interconnect),
installed from its repository):

```python
from majorana_lattice import archaeo_executor
res = run(p, shots=1000, executor=archaeo_executor())   # target "stim:stabilizer"
```

### The interactive engine

`MatchingCode` simulates as it goes, one shot at a time: every measurement
returns its outcome at once. The explorer and `js/engine.js` work this way.
It is also the reference the circuit engine is tested against.

```python
from majorana_lattice import MatchingCode, star

code = MatchingCode(star(4, 4), seed=1)
outcome = code.measure_link(0)             # +1, or -1 (a pair of fermions)
record = code.record()                     # replay with MatchingCode.replay
```

## Explorer

`explorer/explorer.html` is a self-contained page: click links to measure
them, apply X, Y or Z errors to create e and m anyons, move Majoranas and
watch fermions appear. It runs `js/engine.js` in the browser. Rebuild it after
changing the engine or the template:

```bash
python explorer/build.py
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

The circuit engine is checked step by step against `MatchingCode`, and every
pair value it reports is confirmed by Stim's own flow checker
(`stim.Circuit.has_flow`). If Node.js is installed, `tests/test_crosscheck.py`
also runs random sessions in the JavaScript engine and replays them in Stim,
comparing every value after every step.

## Citing

If you use this in research, please cite the papers it is built on; see
`CITATION.cff`.

- J. R. Wootton, *A family of stabilizer codes for D(Z₂) anyons and Majorana
  modes*, [arXiv:1501.07779](https://arxiv.org/abs/1501.07779).
- J. R. Wootton, *Measurements of Floquet code plaquette stabilizers*,
  [arXiv:2210.13154](https://arxiv.org/abs/2210.13154).

## Licence

Apache 2.0. See `LICENSE`.
