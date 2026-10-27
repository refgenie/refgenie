#!/usr/bin/env node
/*
 * stamp-build: writes ../refgenie/server/webui/build-info.json after `vite
 * build` has produced the bundle. Run as the last step of `npm run build`.
 *
 * uv_build ships whatever is sitting in refgenie/server/webui/ at package-build
 * time, .gitignore notwithstanding -- so a maintainer's months-old local build
 * can ride into a release wheel undetected. This stamp is what lets
 * tests/scripts/check-wheel-assets.py (--require-fresh) and a deployed
 * server's /service-info (refgenie.web_ui) tell a fresh build from a stale
 * one, by commit rather than by trust.
 */

import { execSync } from 'node:child_process';
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const FRONTEND_ROOT = resolve(HERE, '..');
const REPO_ROOT = resolve(FRONTEND_ROOT, '..');
const WEBUI_DIR = join(REPO_ROOT, 'refgenie', 'server', 'webui');
const OUT_PATH = join(WEBUI_DIR, 'build-info.json');

function readVersion() {
  const pyproject = readFileSync(join(REPO_ROOT, 'pyproject.toml'), 'utf8');
  const match = pyproject.match(/^version\s*=\s*"([^"]+)"/m);
  return match ? match[1] : 'unknown';
}

function run(cmd) {
  try {
    return execSync(cmd, { cwd: REPO_ROOT, stdio: ['ignore', 'pipe', 'ignore'] })
      .toString()
      .trim();
  } catch {
    return null;
  }
}

function readCommit() {
  // CI sets GITHUB_SHA; a local build falls back to git, and finally to
  // 'unknown' rather than failing the build over provenance metadata.
  return process.env.GITHUB_SHA ?? run('git rev-parse HEAD') ?? 'unknown';
}

function readDirty() {
  // No git available (e.g. a source tarball with no .git) reads as dirty:
  // provenance cannot be verified, so --require-fresh must reject it, not
  // silently pass it.
  const status = run('git status --porcelain frontend');
  return status === null || status.length > 0;
}

function main() {
  if (!existsSync(WEBUI_DIR)) {
    console.error(`stamp-build: ${WEBUI_DIR} does not exist -- run \`vite build\` first`);
    process.exit(1);
  }

  const buildInfo = {
    version: readVersion(),
    commit: readCommit(),
    dirty: readDirty(),
    built_at: new Date().toISOString(),
  };

  writeFileSync(OUT_PATH, JSON.stringify(buildInfo, null, 2) + '\n');
  console.log(`stamp-build: wrote ${OUT_PATH}`);
  console.log(`  ${JSON.stringify(buildInfo)}`);
}

main();
