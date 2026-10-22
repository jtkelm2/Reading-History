// Geometry of the degree-4 encyclopedia: a port of src/geom.py (Num) and src/cells.py replay/lines.
// Directions in radians; a card's type is replayed in its own frame and mapped to the stratum's frame.
const EPS = 1e-7, PI = Math.PI, TP = 2 * Math.PI;
function sgn(x) { return Math.abs(x) < EPS ? 0 : (x > 0 ? 1 : -1); }
function floorDiv(x, p) { let q = Math.floor(x / p); if (p - (x - q * p) < EPS) q += 1; return q; }
function nrm(x) { return x - floorDiv(x, TP) * TP; }
function modpi(x) { return x - floorDiv(x, PI) * PI; }
function eq(x, y) { return sgn(nrm(x - y)) === 0; }
function between(x, a, b) {
  if (eq(a, b)) return !eq(x, a);
  const da = nrm(x - a), db = nrm(b - a);
  return sgn(da) > 0 && sgn(db - da) > 0;
}
function sortDirs(xs) { return xs.slice().sort((x, y) => sgn(x - y)); }
function maps(S) {
  if (!S.length) return [[1, 0]];
  const out = [[1, 0]];
  for (const c of S.slice(1)) { const [s, o] = out[out.length - 1]; out.push([-s, nrm(2 * s * c + o)]); }
  const [s, o] = out[out.length - 1];
  return (-s === 1 && eq(nrm(2 * s * S[0] + o), 0)) ? out : null;
}
function apply(f, x) { return nrm(f[0] * x + f[1]); }
function invApply(f, y) { return nrm(f[0] * (y - f[1])); }
function sectorOf(x, S) {
  if (!S.length) return 0;
  const m = S.length;
  for (let j = 0; j < m; j++) if (eq(x, S[j]) || between(x, S[j], S[(j + 1) % m])) return j;
  throw new Error('no sector');
}
// cells.lines: the scaffold lines of branch u at the vertex a: {step: direction mod pi}
function lines(card, a, u) {
  const sc = card.scaffold, steps = sc.steps, piv = sc.pivots, lam = new Array(steps.length).fill(null);
  for (let i = piv.length - 1; i >= 0; i--) {
    const [j, c, d] = piv[i];
    let x = u[j] * PI;
    c.forEach((cq, q) => { if (cq) x -= cq * a[q + 1]; });
    for (let q = 0; q < steps.length; q++) if (q !== j && d[q]) x -= d[q] * lam[q];
    lam[j] = x / d[j];
  }
  const out = {}; steps.forEach((st, j) => { out[st] = modpi(lam[j]); }); return out;
}
// cells.replay (numeric, without (C)): directions of all rays, or null
function replay(card, a, lams) {
  const rays = card.rays, N = rays.length;
  const d = rays.map(r => r[0] === 'v' ? a[r[1]] : null);
  const status = new Array(N).fill(0);
  for (let i = 0; i < card.steps.length; i++) {
    const T = card.steps[i];
    const S = []; for (let j = 0; j < N; j++) if (status[j] === 1) S.push(j);
    const Sd = sortDirs(S.map(j => d[j]));
    const phis = maps(Sd);
    if (phis === null) return null;
    let anchor = T.find(j => rays[j][0] === 'v' || status[j] === 1);
    if (anchor === undefined) {
      anchor = T[0]; let x = lams[i];
      if (S.length) {
        const below = S.filter(q => q < anchor), above = S.filter(q => q > anchor);
        const p = below.length ? Math.max(...below) : Math.max(...S), q = above.length ? Math.min(...above) : Math.min(...S);
        const c = [x, x + PI].filter(y => between(y, d[p], d[q]));
        if (c.length !== 1) return null;
        x = c[0];
      }
      d[anchor] = x;
    }
    const lam = S.length ? apply(phis[sectorOf(d[anchor], Sd)], d[anchor]) : d[anchor];
    for (const j of T) {
      if (j !== anchor && status[j] === 0 && rays[j][0] === 'h') {
        if (!S.length) { d[j] = nrm(lam + PI); continue; }
        const below = S.filter(q => q < j), above = S.filter(q => q > j);
        const p = below.length ? Math.max(...below) : Math.max(...S), q = above.length ? Math.min(...above) : Math.min(...S);
        const f = phis[sectorOf(d[p], Sd)];
        const c = [invApply(f, lam), invApply(f, nrm(lam + PI))].filter(x => between(x, d[p], d[q]));
        if (c.length !== 1) return null;
        d[j] = c[0];
      }
    }
    const imgs = T.map(j => S.length ? modpi(apply(phis[sectorOf(d[j], Sd)], d[j])) : modpi(d[j]));
    if (imgs.some(x => !eq(modpi(x - imgs[0]), 0) && !eq(modpi(x - imgs[0]), PI))) return null;
    for (const j of T) status[j] = status[j] === 0 ? 1 : 2;
    const now = []; for (let j = 0; j < N; j++) if (status[j] === 1) now.push(d[j]);
    if (maps(sortDirs(now)) === null) return null;
  }
  for (let j = 0; j < N - 1; j++) if (sgn(d[j + 1] - d[j]) <= 0) return null;
  if (sgn(TP - d[N - 1]) <= 0) return null;
  return d;
}
// the dihedral group on (x, y) = (theta_1, theta_2) in units of pi
function act(g, p) {
  const [x, y] = p, s = [x, y, 1 - x, 1 - y], [r, m] = g, md = k => ((k % 4) + 4) % 4;
  return m === 0 ? [s[md(r)], s[md(r + 1)]] : [s[md(r - 1)], s[md(r - 2)]];
}
const GROUP = [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [3, 1]];
function inverse(g) {
  const probe = [[1 / 7, 2 / 9], [3 / 11, 1 / 5]];
  return GROUP.find(h => probe.every(p => { const q = act(h, act(g, p)); return Math.abs(q[0] - p[0]) + Math.abs(q[1] - p[1]) < 1e-12; }));
}
function raysOf(Q) { return [0, Q[0] * PI, (Q[0] + Q[1]) * PI, (1 + Q[1]) * PI]; }
// directions (stratum frame) of every ray of a card's type at the stratum point Q, or null
function cardDirs(card, Q) {
  const g = card.g, q = act(inverse(g), Q), a = raysOf(q);
  let lams = null;
  if (card.scaffold) {
    const L = lines(card, a, card.u), byStep = {};
    for (const k in L) byStep[k] = L[k];
    lams = byStep;
  }
  const d = replay(card, a, lams);
  if (d === null) return null;
  const ar = a[g[0]];
  return d.map(x => g[1] === 0 ? nrm(x - ar) : nrm(ar - x));
}
if (typeof module !== 'undefined') module.exports = { cardDirs, act, inverse, GROUP, replay, lines, raysOf, nrm, PI, TP };
