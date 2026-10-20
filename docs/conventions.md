# Conventions

This page records the choices the engine makes. Several of them were found
the hard way, and getting them wrong gives results that look plausible but
aren't physical.

## Lattice and operators

A lattice has one qubit per vertex. Every vertex has three links, labelled x,
y and z, with no two links of the same label meeting at a vertex. Lattices
live on a torus; each link stores the displacement between its ends, so faces
are traced from the geometry.

- A **link operator** on an x link between qubits j and k is X_j X_k, and
  similarly for y and z.
- A **plaquette operator** is the product of the link operators around a
  face. At each vertex of the face this leaves the Pauli of the one link that
  is not part of the face, so the engine stores it as that tensor product,
  with sign +1.
- A **string operator** is the product of link operators along a path. Two
  strings anticommute only when they share exactly one endpoint.

## Addresses

Preset lattices give every vertex, link and plaquette an address: a unit
cell (i, j) in axial coordinates on the lattice vectors a1 = (√3, 0) and
a2 = (√3/2, 3/2), plus a name within the cell. Each honeycomb cell has a
down site D and an up site U one unit above it, joined by the cell's z link.
On the star lattice, each site becomes a triangle:

| Address | What it is |
| --- | --- |
| `Vertex(i, j, "D.x")` | corner of triangle D that the x link leaves from (`D.x D.y D.z U.x U.y U.z`) |
| `Link(i, j, "z")` | D(i, j) to U(i, j); `"x"` is D(i, j) to U(i+1, j−1), `"y"` is D(i, j) to U(i, j−1) |
| `Link(i, j, "D.z")` | the z edge of triangle D, between its x and y corners |
| `Plaquette(i, j, "D")` | triangle D; `"U"` is triangle U |
| `Plaquette(i, j, "hex")` | the large face centred at (i + 1/3, j + 1/3) |

On the honeycomb the sites are plain `Vertex(i, j, "D")` and `"U"`; on the
modified honeycomb only D becomes a triangle. A dimer is not a fixed element,
so it is addressed by its link.

Cells wrap round the torus. Because a2 leans right, (i, j + rows) is the same
cell as (i + rows/2, j). `lat.index(address)` gives the vertex, edge or face
id, `lat.address(kind, id)` goes back, `lat.coords(address)` gives an
element's centre in axial coordinates, `lat.nearest(kind, (u, v))` finds the
closest element of a kind and `lat.distance(a, b)` measures on the torus.
Every operation and query accepts an address or a plain index, and returns
addresses on preset lattices.

## Programs

An experiment is a `Program`: an immutable value holding the circuit so far
and the bookkeeping needed to decide what can be done next (who is paired
with whom, which plaquettes hold anyons, and what each measurement record
means). Operations take a Program and return a new one; `run` is the only
step that executes anything.

That works because nothing the bookkeeping decides depends on a measurement
outcome: which pairs a link measurement rewires, the moves on offer and the
plaquette values all follow from the operations alone. Outcomes only change
the values of pairs, and those are stored as expressions (below).

## Pairs: dimers and Majoranas

Every vertex is an endpoint of exactly one tracked pair. Each pair has two
endpoints, a path of links between them, and a signed Pauli operator: its
fermion parity.

- A **dimer** is a pair whose path is a single link and which is a
  stabilizer (a "string"). It reads +1 when empty and -1 when it holds a
  fermion.
- Every vertex that is not in a dimer is a **Majorana**. These are the ends
  of longer strings (created by measuring links) and the ends of released
  dimers. A Majorana pair's parity is shared by both Majoranas and is not
  located anywhere; it becomes visible when the two are fused into a dimer.
- The path of a Majorana pair is mostly bookkeeping. It does decide the e/m
  colouring (below): an anyon carried across it switches between e and m.
- A pair made by measuring links is a string, so if moves bring its two ends
  next to each other and its path collapses to that one link, it is a dimer:
  the two Majoranas have fused by themselves.

`release` (and `split_dimer`, which takes the dimer's link) is bookkeeping
only. It relabels a dimer's ends as Majoranas so they can be moved, adds no
gates, and is marked in the circuit by a tagged `I`.

## Pair values are parities of records

The sign of every tracked operator is fixed by how it came about, so that
vacuum always reads +1. Each pair also stores which measurement records its
value is the parity of, plus a sign that flips when a Pauli is applied.

**Measuring a link** between i and j, where i is paired with a and j with b
(operators A and B), gives:

- the new dimer (i, j), with operator +K_ij, whose value is the new record k;
- the new pair (a, b), with operator A B K_ij and path a to i to j to b,
  whose value is (value of A) × (value of B) × record k.

So a -1 outcome puts a fermion on both new pairs; fermions appear in pairs.
If i and j are already paired, the link commutes with everything and the
pair becomes the dimer (i, j), with value record k.

The same rule covers every case: between two dimers it creates a Majorana
pair; next to a Majorana it moves the Majorana two vertices along; between
two Majoranas it fuses them.

**Applying a Pauli** (a link operator, or one qubit's X, Y or Z) flips the
sign of every pair and plaquette it anticommutes with.

Every tracked pair therefore stays determined by the records, and `run`
reports each pair's value in every shot without measuring it. The tests
confirm each one, sign included, with `stim.Circuit.has_flow`.

**Moving a Majorana** (`move_majorana`) is a link measurement next to a
Majorana, followed by a fix. `fix_value`, from -1 to 1, is the expected value
of the dimer the move leaves behind: it reads +1 with probability
(1 + fix_value) / 2 in every shot, so 1 always empties it (the default), -1
always leaves a fermion and 0 is a fair coin. The far pair always matches it,
so fermions still come in pairs.

Every value is built the same way. First the paper's fix: the link operators
along the old dimer, conditioned on the outcome record k (`CZ[fix:7]
rec[-1] ...`), which makes both new pairs +1 and removes k from them. Then a
coin: a record c that is 1 with probability (1 - fix_value) / 2, from a
dedicated qubit after the lattice's (`R[coin:12]` and
`HERALDED_PAULI_CHANNEL_1[coin:12](p, 0, 0, 0)`). The same links,
conditioned on c (`CZ[flip:7]`), flip both pairs back, so their values
become c. The coin is the only randomness, and it is independent of the
measurement. On hardware it is a qubit rotated by RY(θ), with
sin²(θ/2) = p, and measured; Stim has no such rotation, so it is written as
the heralded channel, which gives the same record. `create_majoranas` and
`route_majorana` take `fix_value` too; `measure_link` without one leaves the
raw outcome.

**Routing** (`route_majorana`) plans a route that never reuses a vertex,
either as a stop or as a vertex moved across. Each move rewires the dimers it
passes, so a route that crosses its own wake would land somewhere other than
planned, or even create a stray Majorana pair.

**Readouts** (`readout(p, name, operator)`) measure any operator at that point
in the circuit, and `run` reports it as `res[name]`. Tracked pairs and
plaquettes that anticommute with it become undetermined afterwards. A readout
of a pair from an earlier stage, as in the braid, is the usual use.
`fuse(p, link, name)` and `measure_link(..., name=...)` name a link outcome
the same way.

## Vacuum

Vacuum means **every dimer and every plaquette reads +1**, where each
plaquette operator has sign +1 as defined above.

An earlier version measured each plaquette once and declared whatever came
out to be "empty". That is fine for counting e and m anyons, but a Majorana
carried around a plaquette sees its real value. Roughly half the plaquettes
then carried hidden flux, and a pair made from vacuum could fuse to a fermion
depending on the route it took. With the definition above, fusion results
are route independent; `test_route_independence` checks this.

`prepare` makes the vacuum the way hardware would: reset every qubit (and
rotate for an x or y background, so that background's links start at +1),
measure every plaquette with `MPP`, then undo each -1 by applying that
plaquette's destabilizer, conditioned on its record. The destabilizers come
from a Stim tableau whose stabilizers are the dimers, the plaquettes, and the
products of the background Pauli that commute with every plaquette. Those
last ones must be pure products of the background Pauli, since the start
state already fixes them. If one were a logical times a plaquette, the
correction for that plaquette would flip the logical, and the torus's
logical qubits would come out wrong. A plaquette fixed by the others may end
at -1; its sign is flipped so vacuum still reads +1.

The result is exactly the state you would get by postselecting every
plaquette onto +1, which the tests check on every preset and background. The
plaquette records stay at the start of the measurement record, labelled
`("prepare", f)`.

## e and m

A plaquette at -1 holds an e or an m, depending on a bicolouring of the
faces. Count how many pair paths use each link. Two faces sharing a link used
an odd number of times get the same colour; faces sharing a link used an
even number of times get different colours. Colour 0 is e and colour 1 is m.

On a torus, a loop of even links can wind all the way round. When that
happens the colouring is only defined locally, and `colouring()` reports it
as inconsistent.

## Moves

Majoranas, fermions and anyons all move the same way. Each has a list of the
moves on offer, as `[(link, destination)]`, and a function that makes one
move across a given link. The function raises `ValueError` for a move that is
not on the list.

| particle | list | make one move | what it does |
| --- | --- | --- | --- |
| Majorana | `majorana_moves(p, v)` | `move_majorana(p, v, link)` | measures the link; the Majorana lands two vertices along |
| fermion | `fermion_moves(p, v)` | `move_fermion(p, v, link)` | applies the link, toggling the fermion on both pairs it joins |
| e or m anyon | `anyon_moves(p, f)` | `move_anyon(p, f, link)` | applies the link's Pauli, flipping the two plaquettes either side |

`route_majorana(p, a, b)` chains Majorana moves into a route.

For anyons, applying a link's own Pauli (X on an x link, and so on) to one of
its endpoints flips exactly the two plaquettes either side of that link
(`create_anyons(p, link)` does just that). A move is on offer only if both
faces have the same colour, so e stays e and m stays m, and the destination
holds no anyon, so nothing is annihilated. Crossing a link whose faces differ
in colour turns e into m and leaves a fermion on the dimer at the endpoint;
that is available through `apply_pauli` but is not offered as a move. In
vacuum the same-colour links are the z links.

## Torus effects

The torus has two logical qubits. If a Majorana pair's path winds round the
torus, fusing it measures a logical operator and the outcome can be random.
Experiments that need clean fusion results should stay within a patch; the
tests and `examples/linking.py` do this by keeping moves inside a disc.

## Saving, replaying and running

Instructions are tagged (`MPP[link:12]`, `X[pauli]`, `CX[fix:7]`,
`R[coin:12]` and `HERALDED_PAULI_CHANNEL_1[coin:12]`, `CX[flip:7]`,
`I[release:3]`, `MPP[readout:name]`, and `prepare` for the vacuum), so
`Program.from_circuit(lattice, circuit)` rebuilds the bookkeeping from the
circuit alone. Saving the circuit text saves the whole program.

`run` takes an executor: any function from a `stim.Circuit` and a shot count
to a bool array of measurement records, one row per shot, passed as
`run(p, executor=...)`. Local Stim is the default.
