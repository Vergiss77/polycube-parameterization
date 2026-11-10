// B-spline evaluation, ported from src/pcp/blocks/bspline.cpp with the same algorithms and loop
// order (Piegl & Tiller A2.2 / A2.3), so that results agree with the C++ pipeline up to round-off
// (design D3). Jacobians are stored column-major: J[3 * d + i] = d x_i / d (u, v, w)_d.

// Knot span s (degree <= s < nBasis) with knots[s] <= u < knots[s + 1]; u at the end of the
// range belongs to the last non-empty span.
export function findSpan(degree, knots, u) {
  const n = knots.length - degree - 1;
  if (u >= knots[n]) return n - 1;
  if (u <= knots[degree]) return degree;
  let lo = degree, hi = n;  // knots[lo] <= u < knots[hi]
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (knots[mid] <= u) lo = mid;
    else hi = mid;
  }
  return lo;
}

// Values N[0..p] and first derivatives dN[0..p] (dN may be null) of the basis functions
// s - p .. s at u.
export function basisFunctions(degree, knots, s, u, N, dN) {
  const p = degree;
  const ndu = [];
  for (let j = 0; j <= p; ++j) ndu.push(new Float64Array(p + 1));
  const left = new Float64Array(p + 1), right = new Float64Array(p + 1);
  ndu[0][0] = 1.0;
  for (let j = 1; j <= p; ++j) {
    left[j] = u - knots[s + 1 - j];
    right[j] = knots[s + j] - u;
    let saved = 0.0;
    for (let r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }
  for (let j = 0; j <= p; ++j) N[j] = ndu[j][p];
  if (!dN) return;
  for (let r = 0; r <= p; ++r) {
    let d = 0.0;
    if (r >= 1) d += ndu[r - 1][p - 1] / ndu[p][r - 1];
    if (r <= p - 1) d -= ndu[r][p - 1] / ndu[p][r];
    dN[r] = p * d;
  }
}

// Distinct knot values (element boundaries) of a knot vector.
export function breakpoints(knots) {
  const b = [knots[0]];
  for (let i = 1; i < knots.length; ++i) if (knots[i] > b[b.length - 1]) b.push(knots[i]);
  return b;
}

export function elementCount(knots) {
  return breakpoints(knots).length - 1;
}

// Basis values of one direction at a list of parameters (reused across a tensor grid).
export function basisTable(degree, knots, params) {
  const m = params.length, p1 = degree + 1;
  const t = { degree, span: new Int32Array(m), N: new Float64Array(m * p1), dN: new Float64Array(m * p1) };
  const N = new Float64Array(p1), dN = new Float64Array(p1);
  for (let i = 0; i < m; ++i) {
    const s = findSpan(degree, knots, params[i]);
    basisFunctions(degree, knots, s, params[i], N, dN);
    t.span[i] = s;
    t.N.set(N, i * p1);
    t.dN.set(dN, i * p1);
  }
  return t;
}

// Evaluation of one patch ({degree, knots, n, coefs} as returned by parseMultiPatch).
export class Volume {
  constructor(patch) {
    this.patch = patch;
    this.degree = patch.degree;
    this.knots = patch.knots;
    this.n = patch.n;
    this.coefs = patch.coefs;
  }

  // Position (x[0..2]) and, if J is given, Jacobian from basis table entries iu, iv, iw.
  evalTables(tu, iu, tv, iv, tw, iw, x, J) {
    const [p0, p1, p2] = this.degree;
    const [nu, nv] = this.n;
    const c = this.coefs;
    const su = tu.span[iu] - p0, sv = tv.span[iv] - p1, sw = tw.span[iw] - p2;
    const ou = iu * (p0 + 1), ov = iv * (p1 + 1), ow = iw * (p2 + 1);
    x[0] = x[1] = x[2] = 0.0;
    if (J) J.fill(0.0);
    for (let k = 0; k <= p2; ++k) {
      const Nw = tw.N[ow + k], dNw = tw.dN[ow + k];
      for (let j = 0; j <= p1; ++j) {
        const Nv = tv.N[ov + j], dNv = tv.dN[ov + j];
        for (let i = 0; i <= p0; ++i) {
          const Nu = tu.N[ou + i];
          const q = 3 * ((su + i) + nu * ((sv + j) + nv * (sw + k)));
          const w = Nu * Nv * Nw;
          x[0] += w * c[q]; x[1] += w * c[q + 1]; x[2] += w * c[q + 2];
          if (J) {
            const wu = tu.dN[ou + i] * Nv * Nw, wv = Nu * dNv * Nw, ww = Nu * Nv * dNw;
            J[0] += wu * c[q]; J[1] += wu * c[q + 1]; J[2] += wu * c[q + 2];
            J[3] += wv * c[q]; J[4] += wv * c[q + 1]; J[5] += wv * c[q + 2];
            J[6] += ww * c[q]; J[7] += ww * c[q + 1]; J[8] += ww * c[q + 2];
          }
        }
      }
    }
  }

  evalPoint(u, v, w, x = new Float64Array(3)) {
    this.evalTables(...this._tables(u, v, w), x, null);
    return x;
  }

  evalJacobian(u, v, w, x = new Float64Array(3), J = new Float64Array(9)) {
    this.evalTables(...this._tables(u, v, w), x, J);
    return { x, J };
  }

  _tables(u, v, w) {
    const t = [u, v, w].map((s, d) => basisTable(this.degree[d], this.knots[d], [s]));
    return [t[0], 0, t[1], 0, t[2], 0];
  }
}

// det J / (|J_u| |J_v| |J_w|), 0 for a degenerate column (as in the C++ code).
export function scaledJacobian(J) {
  const l0 = Math.hypot(J[0], J[1], J[2]), l1 = Math.hypot(J[3], J[4], J[5]), l2 = Math.hypot(J[6], J[7], J[8]);
  const l = l0 * l1 * l2;
  if (!(l > 0)) return 0.0;
  const det = J[0] * (J[4] * J[8] - J[7] * J[5]) - J[3] * (J[1] * J[8] - J[7] * J[2]) + J[6] * (J[1] * J[5] - J[4] * J[2]);
  return det / l;
}

// Uniform parameter values i / (m - 1), i = 0..m-1 (as the C++ evaluation grids).
export function uniformParams(m) {
  const t = new Float64Array(m);
  for (let i = 0; i < m; ++i) t[i] = i / (m - 1);
  return t;
}
