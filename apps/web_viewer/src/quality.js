// Scaled Jacobian statistics on the same evaluation grids as the C++ pipeline (design D5):
// every patch, n_elements * 4 + 1 points per direction, points ordered patch by patch with u
// fastest (the order of 07_bspline_samples.vtu).
import { Volume, basisTable, elementCount, scaledJacobian, uniformParams } from './bspline.js';

export function evalGridSize(patch) {
  return patch.knots.map((k) => elementCount(k) * 4 + 1);
}

// Calls fn(patchIndex, i, j, k, x, J) for every evaluation point.
export function forEachEvalPoint(mp, fn) {
  const x = new Float64Array(3), J = new Float64Array(9);
  mp.patches.forEach((patch, b) => {
    const V = new Volume(patch);
    const m = evalGridSize(patch);
    const t = m.map((mi, d) => basisTable(patch.degree[d], patch.knots[d], uniformParams(mi)));
    for (let k = 0; k < m[2]; ++k)
      for (let j = 0; j < m[1]; ++j)
        for (let i = 0; i < m[0]; ++i) {
          V.evalTables(t[0], i, t[1], j, t[2], k, x, J);
          fn(b, i, j, k, x, J);
        }
  });
}

// { evalPoints, nonpositive, blocksWithNonpositive, min, mean, perBlock: [{min, nonpositive}],
//   bad: Float64Array of positions of the non-positive points }
export function jacobianQuality(mp) {
  let count = 0, nonpositive = 0, sum = 0.0, min = Infinity;
  const perBlock = mp.patches.map(() => ({ min: Infinity, nonpositive: 0 }));
  const bad = [];
  forEachEvalPoint(mp, (b, i, j, k, x, J) => {
    const sj = scaledJacobian(J);
    sum += sj;
    ++count;
    if (sj < min) min = sj;
    if (sj < perBlock[b].min) perBlock[b].min = sj;
    if (sj <= 0.0) {
      ++nonpositive;
      ++perBlock[b].nonpositive;
      bad.push(x[0], x[1], x[2]);
    }
  });
  return {
    evalPoints: count,
    nonpositive,
    blocksWithNonpositive: perBlock.filter((q) => q.nonpositive > 0).length,
    min: count > 0 ? min : 0,
    mean: count > 0 ? sum / count : 0,
    perBlock,
    bad: Float64Array.from(bad),
  };
}
