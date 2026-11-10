// Evaluation and Jacobian statistics of the web viewer against the C++ pipeline.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseMultiPatch } from '../../apps/web_viewer/src/gismo.js';
import { basisFunctions, findSpan, scaledJacobian, Volume } from '../../apps/web_viewer/src/bspline.js';
import { forEachEvalPoint, jacobianQuality } from '../../apps/web_viewer/src/quality.js';
import { readFixture, readSamples } from './fixture.mjs';

const mp = parseMultiPatch(readFixture('07_bspline.xml'));
const stats = JSON.parse(readFixture('fit_stats.json'));

function openUniform(p, ne) {
  const k = [];
  for (let i = 0; i < p; ++i) k.push(0);
  for (let e = 0; e <= ne; ++e) k.push(e / ne);
  for (let i = 0; i < p; ++i) k.push(1);
  return Float64Array.from(k);
}

test('partition of unity and derivatives match finite differences', () => {
  for (let p = 1; p <= 4; ++p)
    for (const ne of [1, 3]) {
      const k = openUniform(p, ne);
      const N = new Float64Array(p + 1), dN = new Float64Array(p + 1), Np = new Float64Array(p + 1), Nm = new Float64Array(p + 1);
      for (let i = 0; i <= 50; ++i) {
        const u = i / 50, s = findSpan(p, k, u);
        assert.ok(s >= p && s < k.length - p - 1);
        basisFunctions(p, k, s, u, N, dN);
        assert.ok(Math.abs(N.reduce((a, b) => a + b, 0) - 1) <= 1e-14);
        assert.ok(Math.abs(dN.reduce((a, b) => a + b, 0)) <= 1e-12);
        const h = 1e-6;
        if (u - h <= 0 || u + h >= 1 || findSpan(p, k, u - h) !== s || findSpan(p, k, u + h) !== s) continue;
        basisFunctions(p, k, s, u + h, Np, null);
        basisFunctions(p, k, s, u - h, Nm, null);
        for (let r = 0; r <= p; ++r) assert.ok(Math.abs((Np[r] - Nm[r]) / (2 * h) - dN[r]) <= 1e-7);
      }
    }
});

test('positions and scaled Jacobians match the C++ sample mesh point by point', () => {
  const { points, scaledJacobian: sjRef } = readSamples();
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < points.length; i += 3)
    for (let d = 0; d < 3; ++d) { lo[d] = Math.min(lo[d], points[i + d]); hi[d] = Math.max(hi[d], points[i + d]); }
  const diag = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  let q = 0, maxPos = 0, maxSj = 0;
  forEachEvalPoint(mp, (b, i, j, k, x, J) => {
    maxPos = Math.max(maxPos, Math.hypot(x[0] - points[3 * q], x[1] - points[3 * q + 1], x[2] - points[3 * q + 2]));
    maxSj = Math.max(maxSj, Math.abs(scaledJacobian(J) - sjRef[q]));
    ++q;
  });
  assert.equal(q, points.length / 3);
  assert.ok(maxPos <= 1e-12 * diag, `max position difference ${maxPos}`);
  assert.ok(maxSj <= 1e-9, `max scaled Jacobian difference ${maxSj}`);
});

test('Jacobian statistics match the C++ report', () => {
  const qj = jacobianQuality(mp);
  assert.equal(qj.evalPoints, stats.eval_points);
  assert.equal(qj.nonpositive, stats.nonpositive_jacobian);
  assert.equal(qj.blocksWithNonpositive, stats.blocks_with_nonpositive);
  assert.ok(Math.abs(qj.min - stats.min_scaled_jacobian) <= 1e-9);
  assert.ok(Math.abs(qj.mean - stats.mean_scaled_jacobian) <= 1e-9);
  assert.equal(qj.bad.length, 3 * qj.nonpositive);
});

test('single point evaluation agrees with the tensor table path', () => {
  const V = new Volume(mp.patches[1]);
  const a = V.evalPoint(0.3, 0.7, 0.45);
  const { x, J } = V.evalJacobian(0.3, 0.7, 0.45);
  assert.deepEqual(Array.from(a), Array.from(x));
  // Finite-difference check of the Jacobian.
  const h = 1e-6;
  for (let d = 0; d < 3; ++d) {
    const t = [0.3, 0.7, 0.45], tm = [...t];
    t[d] += h; tm[d] -= h;
    const xp = V.evalPoint(...t), xm = V.evalPoint(...tm);
    for (let c = 0; c < 3; ++c) assert.ok(Math.abs((xp[c] - xm[c]) / (2 * h) - J[3 * d + c]) <= 1e-6);
  }
});
