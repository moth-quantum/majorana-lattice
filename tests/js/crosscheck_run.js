// Random sessions in the browser engine, saved for replay in Stim.
// usage: node crosscheck_run.js lattices.json out.json
const path = require("path");
const { MatchingCode } = require(path.join(__dirname, "..", "..", "js", "engine.js"));
const lats = JSON.parse(require("fs").readFileSync(process.argv[2], "utf8"));
const lat = { ...lats["star 4x4"], name: "star 4x4" };
const runs = [];
for (let s = 1; s <= 12; s++) {
  const c = new MatchingCode(lat, { seed: s });
  let r = s * 7919;
  const rnd = () => (r = (r * 16807) % 2147483647) / 2147483647;
  const snaps = [];
  const snap = () => snaps.push({
    steps: c.log.length,
    plaq: c.plaqValues(),
    pairs: [...c.pairs.entries()].map(([id, P]) => [Math.min(P.a, P.b), Math.max(P.a, P.b), P.kind, c.pairValue(id)]).sort((x, y) => x[0] - y[0]),
    col: c.colouring(),
  });
  if (s % 3 === 0) c.force = 1;
  for (let t = 0; t < 120; t++) {
    const u = rnd();
    if (u < 0.55) c.measureLink(Math.floor(rnd() * c.edges.length));
    else if (u < 0.65) c.applyLink(Math.floor(rnd() * c.edges.length));
    else if (u < 0.72) c.applyPauli(Math.floor(rnd() * c.n), "XYZ"[Math.floor(rnd() * 3)]);
    else if (u < 0.77) c.release(Math.floor(rnd() * c.n));
    else if (u < 0.8) c.measureLabel("xyz"[Math.floor(rnd() * 3)]);
    else {
      const ms = [...c.pairs.values()].filter(P => P.kind === "majorana").flatMap(P => [P.a, P.b]);
      if (ms.length) { const q = ms[Math.floor(rnd() * ms.length)]; const o = c.hopOptions(q); if (o.length) c.hop(q, o[Math.floor(rnd() * o.length)].edge); }
    }
    snap();
  }
  runs.push({ record: c.record(), snaps });
}
require("fs").writeFileSync(process.argv[3], JSON.stringify(runs));
console.log("runs", runs.length, "steps", runs.map(r => r.record.log.length).join(","));
