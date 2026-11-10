// Block coloring and color map of the web viewer.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseMultiPatch } from '../../apps/web_viewer/src/gismo.js';
import { blockColoring, viridis } from '../../apps/web_viewer/src/colors.js';
import { readFixture } from './fixture.mjs';

test('patches sharing an interface get different colors', () => {
  const mp = parseMultiPatch(readFixture('07_bspline.xml'));
  const c = blockColoring(mp);
  assert.equal(c.length, mp.patches.length);
  for (const it of mp.interfaces) assert.notEqual(c[it[0]], c[it[2]]);
});

test('a dense synthetic adjacency still gets a proper coloring', () => {
  // 3 x 3 x 3 grid of patches, face neighbors.
  const id = (i, j, k) => i + 3 * (j + 3 * k);
  const interfaces = [];
  for (let k = 0; k < 3; ++k)
    for (let j = 0; j < 3; ++j)
      for (let i = 0; i < 3; ++i) {
        if (i < 2) interfaces.push([id(i, j, k), 2, id(i + 1, j, k), 1]);
        if (j < 2) interfaces.push([id(i, j, k), 4, id(i, j + 1, k), 3]);
        if (k < 2) interfaces.push([id(i, j, k), 6, id(i, j, k + 1), 5]);
      }
  const c = blockColoring({ patches: new Array(27).fill(null), interfaces });
  for (const it of interfaces) assert.notEqual(c[it[0]], c[it[2]]);
});

test('viridis end points and clamping', () => {
  assert.deepEqual(viridis(0), [0.267, 0.005, 0.329]);
  assert.deepEqual(viridis(1), [0.993, 0.906, 0.144]);
  assert.deepEqual(viridis(-3), viridis(0));
  assert.deepEqual(viridis(7), viridis(1));
});
