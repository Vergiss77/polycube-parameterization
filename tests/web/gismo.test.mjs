// XML parser of the web viewer (apps/web_viewer/src/gismo.js).
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseMultiPatch, GismoError } from '../../apps/web_viewer/src/gismo.js';
import { readFixture } from './fixture.mjs';

const xml = readFixture('07_bspline.xml');
const stats = JSON.parse(readFixture('fit_stats.json'));

test('parses the multipatch written by the C++ pipeline', () => {
  const mp = parseMultiPatch(xml);
  assert.equal(mp.patches.length, stats.blocks);
  assert.equal(mp.interfaces.length, stats.interior_faces);
  assert.equal(mp.boundary.length, stats.boundary_faces);
  for (const p of mp.patches) {
    assert.deepEqual(p.degree, [3, 3, 3]);
    for (let d = 0; d < 3; ++d) {
      assert.equal(p.knots[d][0], 0);
      assert.equal(p.knots[d][p.knots[d].length - 1], 1);
    }
    assert.equal(p.coefs.length, 3 * p.n[0] * p.n[1] * p.n[2]);
  }
  // Interface lines as written by gismo_topology: identity map, normal flag 0.
  for (const it of mp.interfaces) {
    assert.deepEqual(it.slice(4, 7), [0, 1, 2]);
    const normal = Math.floor((it[1] - 1) / 2);
    assert.equal(it[7 + normal], 0);
  }
});

test('first coefficient matches the file text exactly', () => {
  const mp = parseMultiPatch(xml);
  const first = /<coefs[^>]*>\s*(\S+)\s+(\S+)\s+(\S+)/.exec(xml);
  assert.deepEqual(Array.from(mp.patches[0].coefs.slice(0, 3)), [Number(first[1]), Number(first[2]), Number(first[3])]);
});

function rejects(text, pattern) {
  assert.throws(() => parseMultiPatch(text), (e) => e instanceof GismoError && pattern.test(e.message));
}

test('rejects files that are not multipatches', () => {
  rejects('<xml><Geometry type="TensorBSpline3" id="0"></Geometry></xml>', /no <MultiPatch>/);
  rejects('hello', /no <MultiPatch>/);
});

test('rejects a patch with too few control points', () => {
  // Drop the last control point line of patch 0.
  const a = xml.indexOf('</coefs>');
  const endOfLast = xml.lastIndexOf('\n', a);
  const startOfLast = xml.lastIndexOf('\n', endOfLast - 1);
  rejects(xml.slice(0, startOfLast + 1) + xml.slice(endOfLast + 1), /patch 0: .* control points, but the knot vectors require/);
});

test('rejects interfaces that refer to a missing patch or side', () => {
  const n = parseMultiPatch(xml).patches.length;
  const withIface = (line) => xml.replace(/<interfaces>\n/, `<interfaces>\n${line}\n`);
  rejects(withIface(`${n + 3} 1 0 2 0 1 2 0 1 1`), new RegExp(`refers to patch ${n + 3}`));
  rejects(withIface('0 7 1 1 0 1 2 0 1 1'), /side 7 of patch 0 is not in 1..6/);
});

test('rejects a side that is neither an interface nor a boundary side', () => {
  const b = /<boundary>\n(\d+ \d+)\n/.exec(xml);
  const [p, s] = b[1].split(' ');
  rejects(xml.replace(`<boundary>\n${b[1]}\n`, '<boundary>\n'), new RegExp(`side ${s} of patch ${p} is neither`));
});

test('rejects decreasing knots and bad geometry ids', () => {
  rejects(xml.replace(/<KnotVector degree="3">0 0 0 0/, '<KnotVector degree="3">0 0 0 0.5'), /knots are decreasing|end knots/);
  rejects(xml.replace('<Geometry type="TensorBSpline3" id="0">', '<Geometry type="TensorBSpline3" id="5">'), /geometry ids must be 0..N-1/);
});
