# majorana-lattice

Matching codes, anyons and Majorana modes on trivalent lattices, simulated
exactly with [Stim](https://github.com/quantumlib/Stim).

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
| `src/majorana_lattice/` | The Python package: lattices and the Stim engine |
| `js/engine.js` | The same engine in plain JavaScript, for browsers and games |
| `explorer/` | An interactive page for exploring the board, and its build script |
| `examples/` | Runnable experiments: the braid, and linked Majorana pairs |
| `tests/` | Physics tests, and a check that the JavaScript engine matches Stim |
| `docs/conventions.md` | How to read the board, and the conventions that matter |

## Install

```bash
pip install git+https://github.com/moth-quantum/majorana-lattice
```

or, from a clone, `pip install -e ".[test]"`.

## Quick start

```python
from majorana_lattice import MatchingCode, star

code = MatchingCode(star(4, 4), seed=1)   # z dimers and plaquettes all empty

# Measure a link between two dimers: a new dimer forms, and a pair of
# Majoranas appears at the two far ends.
edge = next(e for e, E in enumerate(code.lat.edges) if E.label == "x")
outcome = code.measure_link(edge)          # +1, or -1 (a pair of fermions)
print(code.majoranas())                    # the two new Majoranas
print(code.summary())                      # counts of fermions, e, m, Majoranas

# Move a Majorana: measures the link and fixes a -1 outcome.
m = code.majoranas()[0]
link, destination = code.majorana_moves(m)[0]
code.move_majorana(m, link)

# Every run is a log you can save and replay exactly.
record = code.record()
same = MatchingCode.replay(code.lat, record)
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

If Node.js is installed, `tests/test_crosscheck.py` also runs random sessions
in the JavaScript engine and replays them in Stim, comparing every value
after every step.

## Citing

If you use this in research, please cite the papers it is built on; see
`CITATION.cff`.

- J. R. Wootton, *A family of stabilizer codes for D(Z₂) anyons and Majorana
  modes*, [arXiv:1501.07779](https://arxiv.org/abs/1501.07779).
- J. R. Wootton, *Measurements of Floquet code plaquette stabilizers*,
  [arXiv:2210.13154](https://arxiv.org/abs/2210.13154).

## Licence

Apache 2.0. See `LICENSE`.
