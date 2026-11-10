// Colors (design D6): block palette with a greedy graph coloring over the interface adjacency,
// and the viridis color map for the scaled Jacobian.

export const PALETTE = [
  [0.36, 0.60, 0.90], [0.95, 0.60, 0.25], [0.40, 0.75, 0.45], [0.85, 0.40, 0.45],
  [0.65, 0.50, 0.85], [0.60, 0.45, 0.35], [0.90, 0.55, 0.75], [0.60, 0.60, 0.60],
  [0.75, 0.75, 0.30], [0.30, 0.75, 0.80], [0.70, 0.85, 0.95], [0.95, 0.80, 0.55],
];

export const BAD = [1.0, 0.18, 0.33];  // non-positive scaled Jacobian
export const PLAIN = [0.78, 0.80, 0.83];

// Color index of every patch such that patches sharing an interface differ (greedy, in patch
// order; uses more than PALETTE.length colors only if needed, then wraps).
export function blockColoring(mp) {
  const N = mp.patches.length;
  const adj = Array.from({ length: N }, () => []);
  for (const it of mp.interfaces) {
    adj[it[0]].push(it[2]);
    adj[it[2]].push(it[0]);
  }
  const color = new Int32Array(N).fill(-1);
  for (let b = 0; b < N; ++b) {
    const used = new Set(adj[b].map((o) => color[o]));
    let c = b % PALETTE.length;  // spread colors, then search for a free one
    for (let k = 0; k < PALETTE.length && used.has(c); ++k) c = (c + 1) % PALETTE.length;
    if (used.has(c)) for (c = PALETTE.length; used.has(c); ++c);
    color[b] = c;
  }
  return color;
}

export function paletteColor(i) {
  return PALETTE[i % PALETTE.length];
}

const VIRIDIS = [
  [0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553],
  [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441],
  [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144],
];

// t in [0, 1] -> [r, g, b].
export function viridis(t) {
  const x = Math.min(1, Math.max(0, t)) * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i;
  return [0, 1, 2].map((c) => VIRIDIS[i][c] * (1 - f) + VIRIDIS[i + 1][c] * f);
}

// CSS linear-gradient of the color map (legend).
export function viridisGradient() {
  return 'linear-gradient(to right, ' + VIRIDIS.map((c, i) =>
    `rgb(${c.map((v) => Math.round(v * 255)).join(',')}) ${Math.round((100 * i) / (VIRIDIS.length - 1))}%`).join(', ') + ')';
}
