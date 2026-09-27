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

## Pairs: dimers and Majoranas

Every vertex is an endpoint of exactly one tracked pair. Each pair has two
endpoints, a path of links between them, and a signed Pauli operator: its
fermion parity.

- A **dimer** is a pair whose path is a single link and which is a
  stabilizer. It reads +1 when empty and -1 when it holds a fermion.
- Every vertex that is not in a dimer is a **Majorana**. These are the ends
  of longer strings (created by measuring links) and the ends of split
  dimers. A Majorana pair's parity is shared by both Majoranas and is not
  located anywhere; it becomes visible when the two are fused into a dimer.
- The path of a Majorana pair is mostly bookkeeping. It does decide the e/m
  colouring (below): an anyon carried across it switches between e and m.

`release` (called "Split dimer" in the explorer) is bookkeeping only. It
relabels a dimer's ends as Majoranas so they can be moved, and does nothing
to the qubits.

## Vacuum

Vacuum means **every dimer and every plaquette reads +1**, where each
plaquette operator has sign +1 as defined above. The engine starts from
|0...0> (which is +1 for every z dimer) and then postselects every plaquette
onto +1. On hardware this would be a measurement followed by corrections.

An earlier version measured each plaquette once and declared whatever came
out to be "empty". That is fine for counting e and m anyons, but a Majorana
carried around a plaquette sees its real value. Roughly half the plaquettes
then carried hidden flux, and a pair made from vacuum could fuse to a fermion
depending on the route it took. With the definition above, fusion results
are route independent; `test_route_independence` checks this.

## Signs follow history

The sign of every tracked operator is fixed by how it came about, so that
vacuum always reads +1.

**Measuring a link** between i and j, where i is paired with a and j with b
(operators A and B), gives:

- the new dimer (i, j), with operator +K_ij, whose value is the outcome;
- the new pair (a, b), with operator A B K_ij and path a to i to j to b,
  whose value is (value of A) × (value of B) × outcome.

So a -1 outcome puts a fermion on both new pairs; fermions appear in pairs.
If i and j are already paired, the link commutes with everything, the outcome
is determined, and the pair becomes the dimer (i, j).

The same rule covers every case: between two dimers it creates a Majorana
pair; next to a Majorana it moves the Majorana two vertices along; between
two Majoranas it fuses them.

**Hopping** (`hop`, "Move Majorana" in the explorer) is a link measurement
next to a Majorana, followed by the paper's fix when the outcome is -1:
applying the link operators along the old dimer, which removes both fermions.

**Routing** (`move_majorana`) plans a route that never reuses a vertex,
either as a stop or as a vertex hopped across. Each hop rewires the dimers it
passes, so a route that crosses its own wake would land somewhere other than
planned, or even create a stray Majorana pair.

## e and m

A plaquette at -1 holds an e or an m, depending on a bicolouring of the
faces. Count how many pair paths use each link. Two faces sharing a link used
an odd number of times get the same colour; faces sharing a link used an
even number of times get different colours. Colour 0 is e and colour 1 is m.

On a torus, a loop of even links can wind all the way round. When that
happens the colouring is only defined locally, and `colouring()` reports it
as inconsistent (the explorer shows a warning).

## Torus effects

The torus has two logical qubits. If a Majorana pair's path winds round the
torus, fusing it measures a logical operator and the outcome can be random.
Experiments that need clean fusion results should stay within a patch; the
tests and `examples/linking.py` do this by keeping moves inside a disc.

## Random outcomes and replay

Setting `force = 1` postselects every random outcome to +1 (the explorer's
"Force random outcomes to +1"). Deterministic outcomes are never changed.

Every run is a log of primitive actions (`measure_link` with its outcome,
`apply_link`, `apply_pauli`, `release`), plus the initial plaquette outcomes.
`MatchingCode.replay` rebuilds a run exactly from `record()`, in Python or
JavaScript, which is also how the two engines are checked against each other.
