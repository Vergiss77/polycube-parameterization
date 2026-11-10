// Shared access to the C++-written fixture (see "web viewer fixture" in tests/unit/test_blocks.cpp).
import { readFileSync, existsSync } from 'node:fs';
import { join } from 'node:path';

export const fixtureDir = process.env.PCP_WEB_FIXTURE_DIR;

export function readFixture(name) {
  if (!fixtureDir) throw new Error('PCP_WEB_FIXTURE_DIR is not set (run through ctest)');
  const p = join(fixtureDir, name);
  if (!existsSync(p)) throw new Error(`fixture file ${p} is missing (run the web_fixture test first)`);
  return readFileSync(p, 'utf8');
}

// Points and the scaled_jacobian point field of 07_bspline_samples.vtu.
export function readSamples() {
  const s = readFixture('07_bspline_samples.vtu');
  const pts = /<Points>\s*<DataArray[^>]*>([\s\S]*?)<\/DataArray>/.exec(s)[1].trim().split(/\s+/).map(Number);
  const sj = /Name="scaled_jacobian"[^>]*>([\s\S]*?)<\/DataArray>/.exec(s)[1].trim().split(/\s+/).map(Number);
  return { points: Float64Array.from(pts), scaledJacobian: Float64Array.from(sj) };
}
