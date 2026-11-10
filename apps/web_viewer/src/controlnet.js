// Control points and control net lines (neighbors along u, v, w within a patch), design D7.
// patchIds: the patches to include (all when omitted).
export function controlNet(mp, patchIds) {
  const ids = patchIds ?? mp.patches.map((_, i) => i);
  let np = 0, ns = 0;
  for (const b of ids) {
    const [nu, nv, nw] = mp.patches[b].n;
    np += nu * nv * nw;
    ns += (nu - 1) * nv * nw + nu * (nv - 1) * nw + nu * nv * (nw - 1);
  }
  const positions = new Float64Array(3 * np), index = new Uint32Array(2 * ns), patch = new Int32Array(np);
  let p0 = 0, s = 0;
  for (const b of ids) {
    const P = mp.patches[b];
    const [nu, nv, nw] = P.n;
    positions.set(P.coefs, 3 * p0);
    patch.fill(b, p0, p0 + nu * nv * nw);
    const id = (i, j, k) => p0 + i + nu * (j + nv * k);
    for (let k = 0; k < nw; ++k)
      for (let j = 0; j < nv; ++j)
        for (let i = 0; i < nu; ++i) {
          if (i + 1 < nu) { index[s++] = id(i, j, k); index[s++] = id(i + 1, j, k); }
          if (j + 1 < nv) { index[s++] = id(i, j, k); index[s++] = id(i, j + 1, k); }
          if (k + 1 < nw) { index[s++] = id(i, j, k); index[s++] = id(i, j, k + 1); }
        }
    p0 += nu * nv * nw;
  }
  return { positions, index, patch, pointCount: np, segmentCount: ns };
}
