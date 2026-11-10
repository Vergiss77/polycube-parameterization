// Parser of the G+Smo MultiPatch XML written by the pipeline (07_bspline.xml, see
// src/pcp/io/gismo_xml.cpp). String based so that it runs in browsers and in Node (design D2).
//
// Result:
//   { patches: [{ id, degree: [p, p, p], knots: [Float64Array x3], n: [nu, nv, nw],
//                 coefs: Float64Array(3 * nu * nv * nw) }],   // u fastest, w slowest
//     interfaces: [[p1, s1, p2, s2, m0, m1, m2, o0, o1, o2], ...],
//     boundary: [[patch, side], ...] }                         // sides 1..6
// Throws GismoError with the offending patch and the reason.

export class GismoError extends Error {
  constructor(message) {
    super(message);
    this.name = 'GismoError';
  }
}

function fail(msg) {
  throw new GismoError(msg);
}

// Content of the first element <tag ...>...</tag> at or after `from`; null if absent.
function element(s, tag, from = 0) {
  let a = from;
  for (;;) {
    a = s.indexOf('<' + tag, a);
    if (a < 0) return null;
    const next = s[a + tag.length + 1];
    if (next === ' ' || next === '>' || next === '\n' || next === '\t') break;
    a += tag.length + 1;
  }
  const gt = s.indexOf('>', a);
  const close = '</' + tag + '>';
  const e = s.indexOf(close, gt);
  if (gt < 0 || e < 0) fail(`<${tag}> is not terminated`);
  return { open: s.slice(a, gt), body: s.slice(gt + 1, e), start: a, end: e + close.length };
}

function attribute(open, name) {
  const m = new RegExp(`\\b${name}="([^"]*)"`).exec(open);
  return m ? m[1] : null;
}

function numbers(text, what) {
  const t = text.trim();
  if (t === '') return [];
  return t.split(/\s+/).map((x) => {
    const v = Number(x);
    if (!Number.isFinite(v)) fail(`invalid number '${x}' in ${what}`);
    return v;
  });
}

function integers(text, what) {
  return numbers(text, what).map((v) => {
    if (!Number.isInteger(v)) fail(`expected an integer in ${what}, got ${v}`);
    return v;
  });
}

function parseKnots(body, open, where) {
  const degree = Number(attribute(open, 'degree'));
  if (!Number.isInteger(degree) || degree < 1) fail(`${where}: invalid degree '${attribute(open, 'degree')}'`);
  const k = numbers(body, `${where} knots`);
  for (let i = 1; i < k.length; ++i) if (k[i] < k[i - 1]) fail(`${where}: knots are decreasing`);
  if (k.length < 2 * (degree + 1)) fail(`${where}: ${k.length} knots are too few for degree ${degree}`);
  for (let i = 0; i <= degree; ++i)
    if (k[i] !== k[0] || k[k.length - 1 - i] !== k[k.length - 1])
      fail(`${where}: end knots must be repeated degree + 1 times (open knot vector)`);
  if (!(k[k.length - 1] > k[0])) fail(`${where}: empty parameter range`);
  return { degree, knots: Float64Array.from(k) };
}

export function parseMultiPatch(text) {
  const mpEl = element(text, 'MultiPatch');
  if (!mpEl) fail('no <MultiPatch> element: not a multipatch file written by the pipeline');
  if (attribute(mpEl.open, 'parDim') !== '3') fail(`MultiPatch parDim must be 3, got '${attribute(mpEl.open, 'parDim')}'`);

  // Geometries before the MultiPatch, ids 0..N-1 in order.
  const patches = [];
  let pos = 0;
  for (;;) {
    const g = element(text, 'Geometry', pos);
    if (!g || g.start > mpEl.start) break;
    pos = g.end;
    const id = Number(attribute(g.open, 'id'));
    const where = `patch ${patches.length}`;
    if (id !== patches.length) fail(`${where}: geometry ids must be 0..N-1 in order, got id '${attribute(g.open, 'id')}'`);
    if (attribute(g.open, 'type') !== 'TensorBSpline3')
      fail(`${where}: unsupported geometry type '${attribute(g.open, 'type')}' (only TensorBSpline3)`);
    const degree = [], knots = [], n = [];
    let kp = 0;
    for (let d = 0; d < 3; ++d) {
      const kv = element(g.body, 'KnotVector', kp);
      if (!kv) fail(`${where}: missing knot vector ${d}`);
      kp = kv.end;
      const r = parseKnots(kv.body, kv.open, `${where} direction ${'uvw'[d]}`);
      degree.push(r.degree);
      knots.push(r.knots);
      n.push(r.knots.length - r.degree - 1);
    }
    const ce = element(g.body, 'coefs');
    if (!ce) fail(`${where}: missing <coefs>`);
    if (attribute(ce.open, 'geoDim') !== '3') fail(`${where}: coefs geoDim must be 3`);
    const c = numbers(ce.body, `${where} coefs`);
    const need = 3 * n[0] * n[1] * n[2];
    if (c.length !== need)
      fail(`${where}: ${c.length / 3} control points, but the knot vectors require ${n[0]} x ${n[1]} x ${n[2]} = ${need / 3}`);
    patches.push({ id, degree, knots, n, coefs: Float64Array.from(c) });
  }
  if (patches.length === 0) fail('no <Geometry> elements before <MultiPatch>');

  const mp = mpEl.body;
  const pe = element(mp, 'patches');
  if (!pe || attribute(pe.open, 'type') !== 'id_range') fail('MultiPatch needs <patches type="id_range">');
  const range = integers(pe.body, 'patches');
  if (range.length !== 2 || range[0] !== 0 || range[1] !== patches.length - 1)
    fail(`patches id_range must be 0 ${patches.length - 1}, got '${pe.body.trim()}'`);

  const N = patches.length;
  const seen = new Uint8Array(N * 6);
  const mark = (p, s, what) => {
    if (!(p >= 0 && p < N)) fail(`${what} refers to patch ${p}, but there are ${N} patches`);
    if (!(s >= 1 && s <= 6)) fail(`${what}: side ${s} of patch ${p} is not in 1..6`);
    if (seen[p * 6 + s - 1]) fail(`side ${s} of patch ${p} appears more than once`);
    seen[p * 6 + s - 1] = 1;
  };

  const interfaces = [];
  const ie = element(mp, 'interfaces');
  if (ie) {
    const v = integers(ie.body, 'interfaces');
    if (v.length % 10 !== 0) fail(`interfaces: ${v.length} numbers is not a multiple of 10`);
    for (let i = 0; i < v.length; i += 10) {
      const it = v.slice(i, i + 10);
      mark(it[0], it[1], `interface ${i / 10}`);
      mark(it[2], it[3], `interface ${i / 10}`);
      interfaces.push(it);
    }
  }
  const boundary = [];
  const be = element(mp, 'boundary');
  if (be) {
    const v = integers(be.body, 'boundary');
    if (v.length % 2 !== 0) fail('boundary: odd number of entries');
    for (let i = 0; i < v.length; i += 2) {
      mark(v[i], v[i + 1], `boundary entry ${i / 2}`);
      boundary.push([v[i], v[i + 1]]);
    }
  }
  for (let p = 0; p < N; ++p)
    for (let s = 1; s <= 6; ++s)
      if (!seen[p * 6 + s - 1]) fail(`side ${s} of patch ${p} is neither an interface nor a boundary side`);
  return { patches, interfaces, boundary };
}
