// Matching-code engine for the browser. Mirrors matching_code.py (the Stim core):
// same pairing rules, same sign conventions, same action log, so any run can
// be replayed exactly in Stim with MatchingCode.replay.

// ---------------------------------------------------------------- PRNG
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------------------------------------------------------- Paulis
// A Pauli is {x, z, r}: bit arrays (x=z=1 means Y) and a sign bit r (1 = minus).
function pauli(n) { return { x: new Uint8Array(n), z: new Uint8Array(n), r: 0 }; }
function pcopy(p) { return { x: p.x.slice(), z: p.z.slice(), r: p.r }; }
function setP(p, q, P) {
  p.x[q] = P === "X" || P === "Y" ? 1 : 0;
  p.z[q] = P === "Z" || P === "Y" ? 1 : 0;
}
function g(x1, z1, x2, z2) {
  if (!x1 && !z1) return 0;
  if (x1 && z1) return z2 - x2;
  if (x1) return z2 * (2 * x2 - 1);
  return x2 * (1 - 2 * z2);
}
function anticommutes(a, b) {
  let s = 0;
  for (let q = 0; q < a.x.length; q++) s ^= (a.x[q] & b.z[q]) ^ (a.z[q] & b.x[q]);
  return s === 1;
}
// a * b for commuting Paulis (real sign)
function pmul(a, b) {
  const n = a.x.length, out = pauli(n);
  let e = 2 * a.r + 2 * b.r;
  for (let q = 0; q < n; q++) {
    e += g(a.x[q], a.z[q], b.x[q], b.z[q]);
    out.x[q] = a.x[q] ^ b.x[q];
    out.z[q] = a.z[q] ^ b.z[q];
  }
  e = ((e % 4) + 4) % 4;
  if (e === 1 || e === 3) throw new Error("pmul of anticommuting Paulis");
  out.r = e === 2 ? 1 : 0;
  return out;
}

// ---------------------------------------------------------------- Tableau
// Aaronson-Gottesman tableau: rows 0..n-1 destabilizers, n..2n-1 stabilizers.
class Tableau {
  constructor(n, rand) {
    this.n = n; this.rand = rand;
    this.rows = [];
    for (let i = 0; i < 2 * n; i++) {
      const p = pauli(n);
      if (i < n) p.x[i] = 1; else p.z[i - n] = 1;   // |0...0>
      this.rows.push(p);
    }
  }
  _rowmul(h, i) { // row h := row i * row h
    const a = this.rows[i], b = this.rows[h], n = this.n;
    let e = 2 * a.r + 2 * b.r;
    for (let q = 0; q < n; q++) {
      e += g(a.x[q], a.z[q], b.x[q], b.z[q]);
      b.x[q] ^= a.x[q]; b.z[q] ^= a.z[q];
    }
    b.r = (((e % 4) + 4) % 4) === 2 ? 1 : 0;
  }
  _deterministic(P) { // returns sign bit of P's value (0: +1, 1: -1)
    const n = this.n;
    let acc = pauli(n);
    for (let i = 0; i < n; i++) if (anticommutes(this.rows[i], P)) {
      const s = this.rows[i + n];
      let e = 2 * s.r + 2 * acc.r;
      for (let q = 0; q < n; q++) {
        e += g(s.x[q], s.z[q], acc.x[q], acc.z[q]);
        acc.x[q] ^= s.x[q]; acc.z[q] ^= s.z[q];
      }
      acc.r = (((e % 4) + 4) % 4) === 2 ? 1 : 0;
    }
    return acc.r ^ P.r;
  }
  peek(P) {
    for (let i = this.n; i < 2 * this.n; i++) if (anticommutes(this.rows[i], P)) return 0;
    return this._deterministic(P) ? -1 : 1;
  }
  // measure P; force = +1/-1 postselects a random outcome
  measure(P, force = null) {
    const n = this.n;
    let p = -1;
    for (let i = n; i < 2 * n; i++) if (anticommutes(this.rows[i], P)) { p = i; break; }
    if (p < 0) return this._deterministic(P) ? -1 : 1;
    for (let i = 0; i < 2 * n; i++)
      if (i !== p && anticommutes(this.rows[i], P)) this._rowmul(i, p);
    this.rows[p - n] = pcopy(this.rows[p]);
    const out = force != null ? force : (this.rand() < 0.5 ? 1 : -1);
    const np = pcopy(P); np.r = P.r ^ (out === -1 ? 1 : 0);
    this.rows[p] = np;
    return out;
  }
  apply(P) { // apply Pauli P as a unitary
    for (const row of this.rows) if (anticommutes(row, P)) row.r ^= 1;
  }
}

// ---------------------------------------------------------------- Code
class MatchingCode {
  constructor(lat, { seed = 1, plaquetteRaw = null } = {}) {
    this.lat = lat;
    this.n = lat.positions.length;
    this.edges = lat.edges.map(([u, v, label, dx, dy]) => ({ u, v, label, d: [dx, dy] }));
    this.incident = Array.from({ length: this.n }, () => []);
    this.edges.forEach((E, e) => { this.incident[E.u].push(e); this.incident[E.v].push(e); });
    this.seed = seed;
    this.tab = new Tableau(this.n, mulberry32(seed));
    this.pairs = new Map(); this.owner = new Array(this.n).fill(-1); this.next = 0;
    this.log = []; this.force = null;
    const seen = new Set();
    for (let q = 0; q < this.n; q++) {
      const e = this.incident[q].find(f => this.edges[f].label === "z");
      if (seen.has(e)) continue;
      seen.add(e);
      const E = this.edges[e];
      this._add({ a: E.u, b: E.v, path: [e], op: this.linkPauli(e), kind: "string" });
    }
    this.plaq = []; this.plaquetteRaw = [];
    lat.faces.forEach((f, k) => {
      const p = pauli(this.n);
      for (const [q, P] of Object.entries(f.pauli)) setP(p, +q, P);
      // vacuum = every plaquette at +1, so empty plaquettes carry no flux
      const raw = this.tab.measure(p, plaquetteRaw ? plaquetteRaw[k] : 1);
      this.plaquetteRaw.push(raw);
      if (raw === -1) p.r ^= 1;
      this.plaq.push(p);
    });
  }
  linkPauli(e) {
    const E = this.edges[e], p = pauli(this.n), L = E.label.toUpperCase();
    setP(p, E.u, L); setP(p, E.v, L);
    return p;
  }
  other(P, q) { return q === P.a ? P.b : P.a; }
  pathFrom(P, q) { return q === P.a ? P.path.slice() : P.path.slice().reverse(); }
  _add(P) { const id = this.next++; this.pairs.set(id, P); this.owner[P.a] = this.owner[P.b] = id; return id; }
  _measure(op) {
    if (this.force != null && this.tab.peek(op) === 0) return this.tab.measure(op, this.force);
    return this.tab.measure(op);
  }

  // ---------------------------------------------------------- actions
  measureLink(e) {
    const E = this.edges[e], i = E.u, j = E.v, K = this.linkPauli(e);
    const pi = this.owner[i], pj = this.owner[j];
    let m;
    if (pi === pj) {
      this.pairs.delete(pi);
      m = this._measure(K);
      this._add({ a: i, b: j, path: [e], op: K, kind: "string" });
    } else {
      const A = this.pairs.get(pi), B = this.pairs.get(pj);
      this.pairs.delete(pi); this.pairs.delete(pj);
      const a = this.other(A, i), b = this.other(B, j);
      const op = pmul(pmul(A.op, B.op), K);
      m = this._measure(K);
      const kind = A.kind === "string" && B.kind === "string" ? "string" : "majorana";
      const path = simplify([...this.pathFrom(A, a), e, ...this.pathFrom(B, j)]);
      this._add({ a: i, b: j, path: [e], op: K, kind: "string" });
      this._add({ a, b, path, op, kind });
    }
    this.log.push({ action: "measure_link", edge: e, outcome: m });
    return m;
  }
  applyLink(e) { this.tab.apply(this.linkPauli(e)); this.log.push({ action: "apply_link", edge: e }); }
  applyPauli(q, P) {
    const p = pauli(this.n); setP(p, q, P);
    this.tab.apply(p); this.log.push({ action: "apply_pauli", qubit: q, pauli: P });
  }
  release(q) {
    const P = this.pairs.get(this.owner[q]);
    P.kind = P.kind === "string" ? "majorana" : "string";
    this.log.push({ action: "release", qubit: q });
  }
  measureLabel(L) {
    const out = [];
    this.edges.forEach((E, e) => { if (E.label === L) out.push(this.measureLink(e)); });
    return out;
  }
  hopOptions(q) {
    const out = [];
    for (const e of this.incident[q]) {
      const j = this.edges[e].u === q ? this.edges[e].v : this.edges[e].u;
      const Pj = this.pairs.get(this.owner[j]);
      if (Pj.kind === "string" && Pj.path.length === 1 && this.owner[j] !== this.owner[q]) out.push({ edge: e, via: j, to: this.other(Pj, j) });
    }
    return out;
  }
  hop(q, e, correct = true) {
    const j = this.edges[e].u === q ? this.edges[e].v : this.edges[e].u;
    const Pj = this.pairs.get(this.owner[j]);
    const k = this.other(Pj, j), old = this.pathFrom(Pj, j);
    const m = this.measureLink(e);
    if (m === -1 && correct) for (const f of old) this.applyLink(f);
    return { to: k, outcome: m };
  }

  // Neighbouring pairs a fermion on q's pair can step to: applying the link between them
  // toggles both, so it only moves the fermion if the target pair definitely has none.
  fermionMoves(q) {
    const id = this.owner[q], P = this.pairs.get(id), out = [];
    if (this.pairValue(id) !== -1) return out;
    for (const end of [P.a, P.b]) for (const e of this.incident[end]) {
      const j = this.edges[e].u === end ? this.edges[e].v : this.edges[e].u;
      if (this.owner[j] !== id && this.pairValue(this.owner[j]) === 1) out.push({ edge: e, to: j, pair: this.owner[j] });
    }
    return out;
  }

  // ---------------------------------------------------------- reading
  isDimer(P) { return P.kind === "string" && P.path.length === 1; }
  majoranas() {
    const out = [];
    for (const P of this.pairs.values()) if (!this.isDimer(P)) out.push(P.a, P.b);
    return out.sort((a, b) => a - b);
  }
  pairValue(id) { return this.tab.peek(this.pairs.get(id).op); }
  plaqValues() { return this.plaq.map(p => this.tab.peek(p)); }
  colouring() {
    const count = new Array(this.edges.length).fill(0);
    for (const P of this.pairs.values()) for (const e of P.path) count[e]++;
    const owners = new Map();
    this.lat.faces.forEach((f, k) => f.edges.forEach(e => {
      if (!owners.has(e)) owners.set(e, []); owners.get(e).push(k);
    }));
    const nb = this.lat.faces.map(() => []);
    for (const [e, fs] of owners) if (fs.length === 2 && fs[0] !== fs[1]) { nb[fs[0]].push([fs[1], e]); nb[fs[1]].push([fs[0], e]); }
    const col = new Array(this.lat.faces.length).fill(-1);
    let ok = true;
    for (let s = 0; s < col.length; s++) {
      if (col[s] >= 0) continue;
      col[s] = 0; const st = [s];
      while (st.length) {
        const f = st.pop();
        for (const [h, e] of nb[f]) {
          const want = count[e] % 2 ? col[f] : 1 - col[f];
          if (col[h] < 0) { col[h] = want; st.push(h); } else if (col[h] !== want) ok = false;
        }
      }
    }
    return { col, ok };
  }
  record() { return { lattice: this.lat.name, seed: this.seed, plaquette_raw: this.plaquetteRaw, log: this.log }; }

  static replay(lat, rec, upTo = rec.log.length) {
    const c = new MatchingCode(lat, { seed: rec.seed, plaquetteRaw: rec.plaquette_raw });
    for (const s of rec.log.slice(0, upTo)) {
      if (s.action === "measure_link") { c.force = s.outcome; c.measureLink(s.edge); c.force = null; }
      else if (s.action === "apply_link") c.applyLink(s.edge);
      else if (s.action === "apply_pauli") c.applyPauli(s.qubit, s.pauli);
      else if (s.action === "release") c.release(s.qubit);
    }
    return c;
  }
}

function simplify(path) {
  const out = [];
  for (const e of path) { if (out.length && out[out.length - 1] === e) out.pop(); else out.push(e); }
  return out;
}

if (typeof module !== "undefined") module.exports = { MatchingCode, Tableau, pauli, setP };
