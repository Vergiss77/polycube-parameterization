// Side discretization and control net of the web viewer.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseMultiPatch } from '../../apps/web_viewer/src/gismo.js';
import { elementCount } from '../../apps/web_viewer/src/bspline.js';
import { boundarySides, interfaceSides, tessellateSides } from '../../apps/web_viewer/src/tessellate.js';
import { controlNet } from '../../apps/web_viewer/src/controlnet.js';
import { readFixture } from './fixture.mjs';

const mp = parseMultiPatch(readFixture('07_bspline.xml'));

function diagonal(pos) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < pos.length; i += 3)
    for (let d = 0; d < 3; ++d) { lo[d] = Math.min(lo[d], pos[i + d]); hi[d] = Math.max(hi[d], pos[i + d]); }
  return Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
}

test('vertex and triangle counts follow the subdivision', () => {
  for (const s of [1, 2, 5]) {
    const T = tessellateSides(mp, boundarySides(mp), s);
    let nv = 0, nt = 0;
    for (const side of T.sides) {
      const a = Math.floor((side.side - 1) / 2), k = mp.patches[side.patch].knots;
      const n1 = elementCount(k[(a + 1) % 3]), n2 = elementCount(k[(a + 2) % 3]);
      assert.equal(side.vertexCount, (n1 * s + 1) * (n2 * s + 1));
      assert.equal(side.triangleCount, 2 * n1 * n2 * s * s);
      nv += side.vertexCount;
      nt += side.triangleCount;
    }
    assert.equal(T.positions.length, 3 * nv);
    assert.equal(T.index.length, 3 * nt);
    assert.ok(T.index.every((i) => i < nv));
    assert.ok(T.elementLines.length > 0 && T.allLines.length >= T.elementLines.length);
  }
});

test('the boundary mesh has no cracks: every side edge vertex is shared with another side', () => {
  const T = tessellateSides(mp, boundarySides(mp), 2);
  const diag = diagonal(T.positions);
  const perim = [];  // [vertex, sideIndex]
  T.sides.forEach((side, si) => {
    const a = Math.floor((side.side - 1) / 2), k = mp.patches[side.patch].knots;
    const m1 = elementCount(k[(a + 1) % 3]) * 2 + 1, m2 = side.vertexCount / m1;
    for (let j = 0; j < m2; ++j)
      for (let i = 0; i < m1; ++i)
        if (i === 0 || j === 0 || i === m1 - 1 || j === m2 - 1) perim.push([side.firstVertex + i + m1 * j, si]);
  });
  let worst = 0;
  for (const [v, si] of perim) {
    let best = Infinity;
    for (const [w, sj] of perim) {
      if (sj === si) continue;
      const d = Math.hypot(T.positions[3 * v] - T.positions[3 * w], T.positions[3 * v + 1] - T.positions[3 * w + 1],
        T.positions[3 * v + 2] - T.positions[3 * w + 2]);
      if (d < best) best = d;
    }
    worst = Math.max(worst, best);
  }
  assert.ok(worst <= 1e-12 * diag, `largest gap ${worst}`);
});

test('both sides of every interface give the same vertices', () => {
  const diag = diagonal(tessellateSides(mp, boundarySides(mp), 1).positions);
  for (const it of mp.interfaces) {
    const A = tessellateSides(mp, [[it[0], it[1]]], 3), B = tessellateSides(mp, [[it[2], it[3]]], 3);
    assert.equal(A.positions.length, B.positions.length);
    let d = 0;
    for (let i = 0; i < A.positions.length; ++i) d = Math.max(d, Math.abs(A.positions[i] - B.positions[i]));
    assert.ok(d <= 1e-12 * diag);
  }
});

test('normals point outwards; interfaces are taken once', () => {
  const T = tessellateSides(mp, boundarySides(mp), 2);
  for (const side of T.sides) {
    const a = Math.floor((side.side - 1) / 2), sign = (side.side - 1) % 2 === 1 ? 1 : -1;
    for (let v = side.firstVertex; v < side.firstVertex + side.vertexCount; ++v) {
      assert.ok(T.normals[3 * v + a] * sign > 0.5);  // near-identity map in the fixture
      assert.ok(Math.abs(Math.hypot(T.normals[3 * v], T.normals[3 * v + 1], T.normals[3 * v + 2]) - 1) < 1e-5);
    }
    // Triangle winding agrees with the normal.
    const t = side.firstTriangle, [i0, i1, i2] = T.index.slice(3 * t, 3 * t + 3);
    const e1 = [0, 1, 2].map((c) => T.positions[3 * i1 + c] - T.positions[3 * i0 + c]);
    const e2 = [0, 1, 2].map((c) => T.positions[3 * i2 + c] - T.positions[3 * i0 + c]);
    const n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
    assert.ok(n[a] * sign > 0);
  }
  assert.equal(interfaceSides(mp).length, mp.interfaces.length);
});

test('control net: points and segments', () => {
  const all = controlNet(mp);
  let np = 0, ns = 0;
  for (const p of mp.patches) {
    const [u, v, w] = p.n;
    np += u * v * w;
    ns += (u - 1) * v * w + u * (v - 1) * w + u * v * (w - 1);
  }
  assert.equal(all.pointCount, np);
  assert.equal(all.segmentCount, ns);
  assert.equal(all.index.length, 2 * ns);
  const one = controlNet(mp, [2]);
  const [u, v, w] = mp.patches[2].n;
  assert.equal(one.pointCount, u * v * w);
  assert.deepEqual(Array.from(one.positions), Array.from(mp.patches[2].coefs));
  assert.ok(one.patch.every((b) => b === 2));
});
