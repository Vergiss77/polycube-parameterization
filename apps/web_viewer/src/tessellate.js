// Discretization of patch sides (boundary sides or interfaces) into triangle meshes (design D4).
// A side (patch b, side s = 1..6) fixes parameter a = floor((s - 1) / 2) at 0 (odd s) or 1 (even
// s); the other two directions d1 = (a + 1) % 3, d2 = (a + 2) % 3 get n_e * subdiv + 1 values,
// element boundaries included, so neighboring patches evaluate shared edges at identical
// parameters (no cracks: their edge control points are shared).
import { Volume, basisTable, breakpoints, scaledJacobian } from './bspline.js';

// Parameter values of one direction: every element split into `subdiv` equal parts.
export function sideParams(knots, subdiv) {
  const b = breakpoints(knots);
  const t = [];
  for (let e = 0; e + 1 < b.length; ++e)
    for (let k = 0; k < subdiv; ++k) t.push(k === 0 ? b[e] : b[e] + (b[e + 1] - b[e]) * (k / subdiv));
  t.push(b[b.length - 1]);
  return { params: Float64Array.from(t), elements: b.length - 1 };
}

// sides: [[patch, side], ...]. Returns typed arrays for one merged mesh:
//   positions (Float64Array, 3 per vertex), normals (Float32Array), patch (Int32Array per vertex),
//   scaledJacobian (Float64Array per vertex), index (Uint32Array, 3 per triangle),
//   elementLines / allLines (Uint32Array, vertex pairs), sides: [{patch, side, firstVertex,
//   vertexCount, firstTriangle, triangleCount}].
export function tessellateSides(mp, sides, subdiv) {
  // Size first, then fill.
  const plans = sides.map(([b, s]) => {
    const patch = mp.patches[b];
    const a = Math.floor((s - 1) / 2), d1 = (a + 1) % 3, d2 = (a + 2) % 3;
    const p1 = sideParams(patch.knots[d1], subdiv), p2 = sideParams(patch.knots[d2], subdiv);
    return { b, s, a, d1, d2, p1, p2, m1: p1.params.length, m2: p2.params.length };
  });
  let nv = 0, nt = 0, nel = 0, nal = 0;
  for (const P of plans) {
    nv += P.m1 * P.m2;
    nt += 2 * (P.m1 - 1) * (P.m2 - 1);
    nal += (P.m1 - 1) * P.m2 + P.m1 * (P.m2 - 1);
    nel += (P.m1 - 1) * (P.p2.elements + 1) + (P.p1.elements + 1) * (P.m2 - 1);
  }
  const out = {
    positions: new Float64Array(3 * nv), normals: new Float32Array(3 * nv), patch: new Int32Array(nv),
    scaledJacobian: new Float64Array(nv), index: new Uint32Array(3 * nt),
    elementLines: new Uint32Array(2 * nel), allLines: new Uint32Array(2 * nal), sides: [],
  };
  let v0 = 0, t0 = 0, ie = 0, ia = 0;
  const x = new Float64Array(3), J = new Float64Array(9), par = [0, 0, 0];
  const volumes = new Map();
  for (const P of plans) {
    const patch = mp.patches[P.b];
    if (!volumes.has(P.b)) volumes.set(P.b, new Volume(patch));
    const V = volumes.get(P.b);
    const high = (P.s - 1) % 2 === 1;
    const lists = [];
    lists[P.a] = Float64Array.of(high ? patch.knots[P.a][patch.knots[P.a].length - 1] : patch.knots[P.a][0]);
    lists[P.d1] = P.p1.params;
    lists[P.d2] = P.p2.params;
    const tab = lists.map((l, d) => basisTable(patch.degree[d], patch.knots[d], l));
    const sign = high ? 1 : -1;
    for (let j = 0; j < P.m2; ++j)
      for (let i = 0; i < P.m1; ++i) {
        par[P.a] = 0; par[P.d1] = i; par[P.d2] = j;
        V.evalTables(tab[0], par[0], tab[1], par[1], tab[2], par[2], x, J);
        const v = v0 + i + P.m1 * j;
        out.positions.set(x, 3 * v);
        // Outward normal: +-(J_d1 x J_d2).
        const ax = J[3 * P.d1], ay = J[3 * P.d1 + 1], az = J[3 * P.d1 + 2];
        const bx = J[3 * P.d2], by = J[3 * P.d2 + 1], bz = J[3 * P.d2 + 2];
        let nx = ay * bz - az * by, ny = az * bx - ax * bz, nz = ax * by - ay * bx;
        const l = Math.hypot(nx, ny, nz);
        if (l > 0) { nx *= sign / l; ny *= sign / l; nz *= sign / l; }
        out.normals[3 * v] = nx; out.normals[3 * v + 1] = ny; out.normals[3 * v + 2] = nz;
        out.patch[v] = P.b;
        out.scaledJacobian[v] = scaledJacobian(J);
      }
    let t = t0;
    for (let j = 0; j + 1 < P.m2; ++j)
      for (let i = 0; i + 1 < P.m1; ++i) {
        const a = v0 + i + P.m1 * j, b = a + 1, c = a + 1 + P.m1, d = a + P.m1;
        // Counter-clockwise seen from outside (e_d1 x e_d2 = e_a).
        if (high) out.index.set([a, b, c, a, c, d], 3 * t);
        else out.index.set([a, c, b, a, d, c], 3 * t);
        t += 2;
      }
    const s = subdiv;
    for (let j = 0; j < P.m2; ++j)
      for (let i = 0; i + 1 < P.m1; ++i) {
        const a = v0 + i + P.m1 * j;
        out.allLines[ia++] = a; out.allLines[ia++] = a + 1;
        if (j % s === 0) { out.elementLines[ie++] = a; out.elementLines[ie++] = a + 1; }
      }
    for (let j = 0; j + 1 < P.m2; ++j)
      for (let i = 0; i < P.m1; ++i) {
        const a = v0 + i + P.m1 * j;
        out.allLines[ia++] = a; out.allLines[ia++] = a + P.m1;
        if (i % s === 0) { out.elementLines[ie++] = a; out.elementLines[ie++] = a + P.m1; }
      }
    out.sides.push({ patch: P.b, side: P.s, firstVertex: v0, vertexCount: P.m1 * P.m2, firstTriangle: t0, triangleCount: t - t0 });
    v0 += P.m1 * P.m2;
    t0 = t;
  }
  return out;
}

// Boundary sides and one side of every interface, as [patch, side] lists.
export function boundarySides(mp) {
  return mp.boundary.map(([b, s]) => [b, s]);
}

export function interfaceSides(mp) {
  return mp.interfaces.map((it) => [it[0], it[1]]);
}
