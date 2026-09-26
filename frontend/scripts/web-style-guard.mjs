#!/usr/bin/env node
/*
 * web-style-guard: mechanical enforcement of the house Web UI Standard.
 *
 * ONE guard for the whole gen2 fleet. It is owned by the web-design-style
 * skill and copied into a project by web-style-sync.mjs -- never edit the
 * project's copy; edit the skill's and re-run sync.
 *
 * Six checks:
 *   a. frozen-file integrity  -- src/styles/{utilities,modal}.css are verbatim
 *      copies of the skill's canonical files (SHA-256, no recorded hashes).
 *   b. source bans            -- no inline style={{ (except a style object whose
 *      keys are ALL custom properties), no Bootstrap/Tailwind import or asset
 *      reference, no !important in src/**\/*.{ts,tsx}.
 *   c. CSS raw-value bans     -- no hex / px / rem literal outside tokens.css.
 *   d. Tailwind-isms          -- theme-color classes, alpha modifiers, state
 *      variants, arbitrary [] values, bg-white/black, config residue.
 *   e. class resolution       -- every className token resolves to a selector
 *      defined in some stylesheet the app actually loads (src/styles/*.css,
 *      anything a module imports, and their @import graph).
 *   f. stack conformance      -- package.json deps satisfy stack.json; banned
 *      deps fail, discouraged deps warn; .node-version >= nodeMin.
 *
 * Usage:
 *   node scripts/web-style-guard.mjs [--path DIR] [--warn-only] [--update]
 *
 *   --path DIR    project root to check (default: the guard's parent dir, or cwd)
 *   --warn-only   report everything but always exit 0 (rollout mode)
 *   --update      re-record the cached canonical copies from src/styles/ after
 *                 a deliberate, reviewed refresh
 *
 * Per-project escapes live in package.json under "webStyleGuard":
 *   { "ignoreClasses": ["ai"], "exclude": ["src/vendor/"] }
 */

import { createHash } from 'node:crypto';
import {
  copyFileSync, existsSync, mkdirSync, readdirSync, readFileSync, statSync,
} from 'node:fs';
import { basename, dirname, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));

/* ------------------------------------------------------------------ args */

function parseArgs(argv) {
  const out = { warnOnly: false, update: false, path: null };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--warn-only') out.warnOnly = true;
    else if (a === '--update') out.update = true;
    else if (a === '--path') out.path = argv[++i];
    else if (a.startsWith('--path=')) out.path = a.slice('--path='.length);
  }
  return out;
}

const ARGS = parseArgs(process.argv.slice(2));

/*
 * Canonical files live either next to the guard (skill layout) or in
 * <guardDir>/.web-style/ (project layout, written by web-style-sync). Looking
 * both up keeps the guard self-contained in CI without a skill checkout.
 */
const CANON_DIRS = [HERE, join(HERE, '.web-style')];

function canonical(name) {
  for (const d of CANON_DIRS) {
    const p = join(d, name);
    if (existsSync(p)) return p;
  }
  return null;
}

const STACK_PATH = canonical('stack.json');
if (!STACK_PATH) {
  console.error('web-style-guard: cannot find stack.json next to the guard or in .web-style/.');
  console.error('  Run web-style-sync.mjs against this project to install it.');
  process.exit(2);
}
const STACK = JSON.parse(readFileSync(STACK_PATH, 'utf8'));

/*
 * Project root. Explicit --path wins; otherwise the guard's parent (it is
 * installed at <root>/scripts/), falling back to cwd when run from the skill.
 */
const ROOT = resolve(
  ARGS.path
    ? ARGS.path
    : basename(HERE) === 'scripts'
      ? resolve(HERE, '..')
      : process.cwd(),
);
const SRC = join(ROOT, 'src');
const STYLES = join(SRC, 'styles');
const TOKENS_FILE = 'src/styles/tokens.css';

const PKG_PATH = join(ROOT, 'package.json');
const PKG = existsSync(PKG_PATH) ? JSON.parse(readFileSync(PKG_PATH, 'utf8')) : {};
const OPTS = PKG.webStyleGuard || {};
const IGNORE_CLASSES = new Set(OPTS.ignoreClasses || []);
const EXCLUDE = OPTS.exclude || [];

const FROZEN = STACK.frozen || [];
const FROZEN_SET = new Set(FROZEN);

/* --------------------------------------------------------------- helpers */

const findings = new Map(); // check label -> string[]

function fail(check, msg) {
  if (!findings.has(check)) findings.set(check, []);
  findings.get(check).push(msg);
}

const warnings = new Map();

function warn(check, msg) {
  if (!warnings.has(check)) warnings.set(check, []);
  warnings.get(check).push(msg);
}

function sha256(path) {
  return createHash('sha256').update(readFileSync(path)).digest('hex');
}

function rel(p) {
  return relative(ROOT, p).split(sep).join('/');
}

function excluded(relPath) {
  return EXCLUDE.some((e) => relPath.startsWith(e) || relPath.includes(e));
}

function walk(dir, exts, out = []) {
  if (!existsSync(dir)) return out;
  for (const entry of readdirSync(dir)) {
    if (entry === 'node_modules' || entry === 'dist' || entry === 'generated') continue;
    const full = join(dir, entry);
    if (statSync(full).isDirectory()) walk(full, exts, out);
    else if (exts.some((e) => entry.endsWith(e))) out.push(full);
  }
  return out;
}

/* ------------------------------------------- (a) frozen-file integrity */

function checkFrozen() {
  const CHECK = 'frozen files';
  for (const relPath of FROZEN) {
    const projFile = join(ROOT, relPath);
    const canonFile = canonical(basename(relPath));
    if (!canonFile) {
      fail(CHECK, `no canonical copy of ${basename(relPath)} available to compare against`);
      continue;
    }
    if (!existsSync(projFile)) {
      fail(CHECK, `${relPath}: missing. Run web-style-sync.mjs.`);
      continue;
    }
    if (sha256(projFile) !== sha256(canonFile)) {
      fail(
        CHECK,
        `${relPath}: differs from the canonical copy in the web-design-style skill.\n` +
          '    Do not edit this file -- project rules belong in src/styles/components.css.\n' +
          '    To restore it:  node <skill>/web-style-sync.mjs --path .',
      );
    }
  }
}

/* --------------------------------------------------- (b)+(c) file bans */

const SOURCE_BANS = [
  // Match the frameworks where they would actually be pulled in (an import, a
  // stylesheet reference, a CDN URL) rather than anywhere the words appear:
  // "bootstrap config" is a legitimate phrase for the startup handshake.
  {
    pattern:
      /(?:from\s*['"]|require\(\s*['"]|@import\s*['"]|url\(\s*['"]?|href=['"]|src=['"])[^'")]*\b(bootstrap|tailwind)/i,
    message: 'Bootstrap/Tailwind import (this project uses utilities.css + BEM only)',
  },
  {
    pattern: /(@tailwind\b|tailwind\.config|bootstrap-icons|bootstrap(\.min)?\.(css|js)\b)/i,
    message: 'Bootstrap/Tailwind asset reference (this project uses utilities.css + BEM only)',
  },
  { pattern: /!important/, message: '!important (raise specificity with a BEM block instead)' },
];

const CSS_BANS = [
  { pattern: /#[0-9a-fA-F]{3,8}\b/, message: 'raw hex color (define an RGB-triplet token in tokens.css)' },
  {
    pattern: /(?<![\w-])\d+(\.\d+)?px\b/,
    message: 'raw px length (use a --space-* / --radius-* token)',
    skipInMedia: true,
  },
  {
    pattern: /(?<![\w-])\d+(\.\d+)?rem\b/,
    message: 'raw rem length (use a --space-* / --font-size-* token)',
    skipInMedia: true,
  },
];

/**
 * Blank out /* ... *\/ comment bodies, keeping newlines so line numbers still
 * line up. Prose about a 1px border is not a 1px border.
 */
function stripCssComments(css) {
  return css.replace(/\/\*[\s\S]*?\*\//g, (m) => m.replace(/[^\n]/g, ' '));
}

/**
 * Report `style={{ ... }}` objects that set real CSS properties.
 *
 * A style object whose keys are ALL custom properties is allowed: it hands a
 * value to CSS rather than bypassing it, and it is the only honest way to
 * express geometry that is computed at runtime (a popover anchored to a
 * selection rect, a bar sized from a ratio). Everything else belongs in a
 * utility class or a BEM block.
 */
function checkInlineStyles(file, relPath) {
  const src = readFileSync(file, 'utf8');
  const re = /style=\{\{/g;
  let m;
  while ((m = re.exec(src)) !== null) {
    let depth = 0;
    let k = m.index + 'style='.length;
    for (; k < src.length; k++) {
      if (src[k] === '{') depth++;
      else if (src[k] === '}' && --depth === 0) break;
    }
    const inner = src.slice(m.index + 'style={'.length, k)
      .replace(/\{[^{}]*\}/g, ' ')  // drop nested objects/expressions
      .replace(/\([^()]*\)/g, ' ');
    const keys = [...inner.matchAll(/(?:^|[,{])\s*(?:(['"`])(.*?)\1|([A-Za-z_$][\w$]*))\s*:/g)]
      .map((km) => km[2] ?? km[3]);
    const spread = /\.\.\./.test(inner);
    const bad = keys.filter((key) => !key.startsWith('--'));
    if (!keys.length || spread || bad.length) {
      const line = src.slice(0, m.index).split('\n').length;
      const what = keys.length
        ? `sets ${bad.length ? bad.join(', ') : 'spread properties'}`
        : 'sets styles';
      fail(
        'source bans',
        `${relPath}:${line}: inline style attribute (${what})\n` +
          '    Use a utility class or a BEM block. For runtime geometry, pass the\n' +
          "    value as a custom property -- style={{ '--x': `${x}px` }} -- and read\n" +
          '    it from CSS with var(--x).',
      );
    }
  }
}

function checkFileBans() {
  for (const file of walk(SRC, ['.ts', '.tsx', '.css'])) {
    const r = rel(file);
    if (FROZEN_SET.has(r) || excluded(r)) continue;
    const isCss = r.endsWith('.css');
    const bans = isCss && r !== TOKENS_FILE ? [...SOURCE_BANS, ...CSS_BANS] : SOURCE_BANS;
    const check = isCss ? 'CSS raw values' : 'source bans';
    if (!isCss) checkInlineStyles(file, r);
    const text = isCss ? stripCssComments(readFileSync(file, 'utf8')) : readFileSync(file, 'utf8');
    text
      .split('\n')
      .forEach((line, i) => {
        // Media-query breakpoints are the one legitimate raw length in CSS:
        // custom properties are not allowed in a media prelude.
        const isMediaPrelude = /^\s*@media\b/.test(line);
        for (const { pattern, message, skipInMedia } of bans) {
          if (isMediaPrelude && skipInMedia) continue;
          if (pattern.test(line)) fail(check, `${r}:${i + 1}: ${message}\n    ${line.trim()}`);
        }
      });
  }
}

/* ------------------------------------------ className extraction (d)+(e) */

/**
 * Pull every className= value out of a source file. Literals yield their raw
 * string; expressions yield each string/template literal found inside, which
 * covers the `cond ? 'a' : 'b'` and `` `x ${y}` `` shapes we actually write.
 */
function classNameValues(src) {
  const out = [];
  const re = /className\s*=\s*/g;
  let m;
  while ((m = re.exec(src)) !== null) {
    let j = m.index + m[0].length;
    if (j >= src.length) continue;
    const c = src[j];
    const line = src.slice(0, j).split('\n').length;
    if (c === '"' || c === "'") {
      const k = src.indexOf(c, j + 1);
      if (k !== -1) out.push({ kind: 'literal', text: src.slice(j + 1, k), line });
    } else if (c === '{') {
      let depth = 0;
      let k = j;
      for (; k < src.length; k++) {
        if (src[k] === '{') depth++;
        else if (src[k] === '}' && --depth === 0) break;
      }
      out.push({ kind: 'expr', text: src.slice(j + 1, k), line });
    }
  }
  return out;
}

const STR_RE = /"([^"]*)"|'([^']*)'|`([^`]*)`/g;
const TOK_RE = /^[a-z0-9][\w:/.\-[\]]*$/;
/** Stands in for a `${...}` hole so an interpolated token stays one token. */
const HOLE = '\u0001';

/**
 * Blank out /* ... *\/ bodies, newlines preserved. A usage example in a JSDoc
 * block is documentation, not markup, and its classes need not resolve.
 */
function stripBlockComments(src) {
  return src.replace(/\/\*[\s\S]*?\*\//g, (m) => m.replace(/[^\n]/g, ' '));
}

/**
 * Blank out string literals that sit on either side of a comparison inside a
 * className expression. `cn('x', align === 'right' && 'x--numeric')` mentions
 * 'right' as data, not as a class, and the fleet writes this constantly.
 */
function stripComparands(expr) {
  const blank = (m) => ' '.repeat(m.length);
  return expr
    .replace(/(?:===|!==|==|!=)\s*(['"`])(?:[^'"`\\]|\\.)*?\1/g, blank)
    .replace(/(['"`])(?:[^'"`\\]|\\.)*?\1\s*(?:===|!==|==|!=)/g, blank);
}

/** Replace every `${...}` hole -- nested braces included -- with HOLE. */
function maskHoles(text) {
  let out = '';
  for (let i = 0; i < text.length; i++) {
    if (text[i] === '$' && text[i + 1] === '{') {
      let depth = 0;
      let j = i + 1;
      for (; j < text.length; j++) {
        if (text[j] === '{') depth++;
        else if (text[j] === '}' && --depth === 0) break;
      }
      out += HOLE;
      i = j;
    } else {
      out += text[i];
    }
  }
  return out;
}

/** token -> [{file, line}] across the whole project. */
function collectClassTokens() {
  const used = new Map();
  for (const file of walk(SRC, ['.tsx'])) {
    const r = rel(file);
    if (r.startsWith('src/styles/') || excluded(r)) continue;
    const src = stripBlockComments(readFileSync(file, 'utf8'));
    const lines = src.split('\n');
    for (const { kind, text, line } of classNameValues(src)) {
      // Commented-out JSX is not shipped markup.
      if (/^\s*\/\//.test(lines[line - 1] ?? '')) continue;
      const segs = kind === 'literal'
        ? [text]
        : [...stripComparands(text).matchAll(STR_RE)].map((sm) => sm[1] ?? sm[2] ?? sm[3] ?? '');
      for (const seg of segs) {
        for (const tok of maskHoles(seg).split(/\s+/)) {
          // A token built by interpolation (`rg-badge--${variant}`) cannot be
          // resolved statically; only its fully-literal siblings are checked.
          if (!tok || tok.includes(HOLE) || !TOK_RE.test(tok)) continue;
          if (!used.has(tok)) used.set(tok, []);
          used.get(tok).push({ file: r, line });
        }
      }
    }
  }
  return used;
}

/* ------------------------------------------------------- (d) Tailwind-isms */

// Sanctioned semantic helpers in components.css end in -fg / -strong / -soft;
// the negative lookahead keeps them out of the theme-color rule.
const TW_RULES = [
  {
    re: /^(?:bg|text|border|ring|from|via|to|fill|stroke)-(?:brand|danger|success|warning|info|ai|gray|primary|secondary|slate|zinc|neutral|stone|red|orange|amber|yellow|lime|green|emerald|teal|cyan|sky|blue|indigo|violet|purple|fuchsia|pink|rose)(?:-\d+)?$(?<!-(?:fg|strong|soft))/,
    message: 'Tailwind theme color class (use a components.css helper or BEM)',
  },
  {
    re: /^(?:bg|text|border|ring|from|via|to|fill|stroke)-[\w-]+\/\d/,
    message: 'Tailwind alpha-modifier color class',
  },
  {
    re: /^(?:hover|focus-visible|focus|active|disabled|group-hover):/,
    message: 'Tailwind state variant prefix (express state in components.css/BEM)',
  },
  {
    re: /^(?:ring-[\w/]|divide-|space-[xy]-|transition-colors$|antialiased$|tracking-)/,
    message: 'Tailwind-only utility (move into a BEM rule)',
  },
  { re: /[a-z0-9]-\[[^\]]+\]/, message: 'Tailwind arbitrary-value bracket' },
  {
    re: /^(?:bg|text)-(?:white|black)$/,
    message: 'raw bg/text-white|black (use .bg-surface / .text-on-brand)',
  },
];

const TW_CONFIGS = [
  'tailwind.config.js', 'tailwind.config.ts', 'tailwind.config.cjs', 'tailwind.config.mjs',
  'postcss.config.js', 'postcss.config.cjs', 'postcss.config.mjs', 'postcss.config.ts',
];

function checkTailwindisms(used) {
  const CHECK = 'Tailwind-isms';
  for (const [tok, locs] of used) {
    for (const { re, message } of TW_RULES) {
      if (re.test(tok)) {
        fail(CHECK, `${tok}  ${formatLocs(locs)}\n    ${message}`);
        break;
      }
    }
  }
  // @tailwind directives anywhere under src, including plain .css.
  for (const file of walk(SRC, ['.css', '.ts', '.tsx'])) {
    const r = rel(file);
    if (excluded(r)) continue;
    readFileSync(file, 'utf8').split('\n').forEach((line, i) => {
      if (/@tailwind\b/.test(line)) fail(CHECK, `${r}:${i + 1}: @tailwind directive`);
    });
  }
  for (const f of TW_CONFIGS) {
    if (existsSync(join(ROOT, f))) fail(CHECK, `Tailwind/PostCSS config still present: ${f}`);
  }
}

function formatLocs(locs) {
  const shown = locs.slice(0, 6).map((l) => `${l.file}:${l.line}`).join(', ');
  return shown + (locs.length > 6 ? ` (+${locs.length - 6} more)` : '');
}

/* ---------------------------------------------------- (e) class resolution */

const SEL_RE = /\.((?:[A-Za-z0-9_-]|\\[:/.[\]])+)/g;

/**
 * Every stylesheet the app actually loads.
 *
 * Seeds from src/styles/*.css plus any .css a module imports directly, then
 * follows @import transitively. A self-contained package that ships its own
 * stylesheet next to its components (snapgate/snapgate.css) is thereby seen,
 * without blessing a .css file nothing imports.
 */
function loadedStylesheets() {
  const seen = new Set();
  const queue = [];

  const enqueue = (path) => {
    if (!path.endsWith('.css') || seen.has(path) || !existsSync(path)) return;
    seen.add(path);
    queue.push(path);
  };

  if (existsSync(STYLES)) {
    for (const f of readdirSync(STYLES)) enqueue(join(STYLES, f));
  }
  // `import './x.css'` from a module -- the other way a stylesheet gets loaded.
  for (const file of walk(SRC, ['.ts', '.tsx'])) {
    const src = readFileSync(file, 'utf8');
    for (const m of src.matchAll(/import\s+['"]([^'"]+\.css)['"]/g)) {
      if (/^[a-z]+:|^~|^[^./]/.test(m[1])) continue; // bare/remote specifier
      enqueue(resolve(dirname(file), m[1]));
    }
  }

  while (queue.length) {
    const path = queue.shift();
    const css = readFileSync(path, 'utf8');
    for (const m of css.matchAll(/@import\s+(?:url\(\s*)?['"]([^'"]+)['"]/g)) {
      if (/^[a-z]+:|^\/\//.test(m[1])) continue; // remote
      enqueue(resolve(dirname(path), m[1]));
    }
  }
  return seen;
}

function definedClasses() {
  const defined = new Set();
  for (const path of loadedStylesheets()) {
    const css = readFileSync(path, 'utf8');
    for (const m of css.matchAll(SEL_RE)) {
      defined.add(
        m[1]
          .replace(/\\:/g, ':').replace(/\\\//g, '/').replace(/\\\./g, '.')
          .replace(/\\\[/g, '[').replace(/\\\]/g, ']'),
      );
    }
  }
  return defined;
}

function checkClassResolution(used) {
  const CHECK = 'class resolution';
  const defined = definedClasses();
  if (!defined.size) {
    warn(CHECK, 'no loaded stylesheet defines any selector; skipping class-resolution check');
    return;
  }
  for (const [tok, locs] of used) {
    if (defined.has(tok) || IGNORE_CLASSES.has(tok)) continue;
    fail(
      CHECK,
      `${tok}  ${formatLocs(locs)}\n` +
        '    not defined in src/styles/*.css (use an existing utility, a components.css\n' +
        '    helper, or add a BEM rule; list intentional non-classes in\n' +
        '    package.json "webStyleGuard".ignoreClasses)',
    );
  }
}

/* ------------------------------------------------------- (f) stack check */

/** Lower bound of a semver range, as [major, minor, patch]. */
function rangeFloor(range) {
  const m = /(\d+)(?:\.(\d+))?(?:\.(\d+))?/.exec(range);
  if (!m) return null;
  return [Number(m[1]), Number(m[2] ?? 0), Number(m[3] ?? 0)];
}

/** Exclusive upper bound of a range, as [major, minor, patch]; null = unbounded. */
function rangeCeil(range) {
  const r = range.trim();
  const floor = rangeFloor(r);
  if (!floor) return null;
  const [maj, min] = floor;
  if (r.startsWith('^')) return maj > 0 ? [maj + 1, 0, 0] : [0, min + 1, 0];
  if (r.startsWith('~')) return [maj, min + 1, 0];
  if (/^\d/.test(r)) {
    // A bare version pins exactly; a bare major/minor is treated as that band.
    const parts = r.split('.').length;
    if (parts >= 3) return [floor[0], floor[1], floor[2] + 1];
    if (parts === 2) return [maj, min + 1, 0];
    return [maj + 1, 0, 0];
  }
  return null; // >=, *, ranges we do not model: unbounded above
}

function cmp(a, b) {
  for (let i = 0; i < 3; i++) if (a[i] !== b[i]) return a[i] - b[i];
  return 0;
}

/** Is `declared` entirely contained in `policy`? */
function satisfies(declared, policy) {
  const dLo = rangeFloor(declared);
  const pLo = rangeFloor(policy);
  if (!dLo || !pLo) return null; // unparseable -> caller warns
  if (cmp(dLo, pLo) < 0) return false;
  const pHi = rangeCeil(policy);
  const dHi = rangeCeil(declared);
  if (!pHi) return true;
  if (!dHi) return false; // declared is open-ended, policy is not
  return cmp(dHi, pHi) <= 0;
}

function allDeps() {
  return { ...(PKG.dependencies || {}), ...(PKG.devDependencies || {}) };
}

function checkStack() {
  const CHECK = 'stack';
  if (!existsSync(PKG_PATH)) {
    fail(CHECK, 'no package.json found at the project root');
    return;
  }
  const deps = allDeps();

  for (const [name, policy] of Object.entries(STACK.require || {})) {
    const declared = deps[name];
    if (!declared) continue; // only enforced when the dep is actually used
    const ok = satisfies(declared, policy);
    if (ok === null) warn(CHECK, `${name}: cannot parse declared range "${declared}" (policy ${policy})`);
    else if (!ok) fail(CHECK, `${name}: declared "${declared}" does not satisfy house policy "${policy}"`);
  }

  for (const [name, why] of Object.entries(STACK.banned || {})) {
    if (deps[name]) fail(CHECK, `${name}: banned dependency (${why})`);
  }

  for (const [name, why] of Object.entries(STACK.discouraged || {})) {
    if (deps[name]) warn(CHECK, `${name}: discouraged dependency (${why})`);
  }

  if (STACK.nodeMin) {
    const min = rangeFloor(STACK.nodeMin);
    const nvPath = join(ROOT, '.node-version');
    if (!existsSync(nvPath)) {
      fail(CHECK, `.node-version missing (house minimum is ${STACK.nodeMin})`);
    } else {
      const declared = rangeFloor(readFileSync(nvPath, 'utf8').trim());
      if (!declared) fail(CHECK, '.node-version is not a parseable version');
      else if (cmp(declared, min) < 0) {
        fail(CHECK, `.node-version ${readFileSync(nvPath, 'utf8').trim()} is below the house minimum ${STACK.nodeMin}`);
      }
    }
  }
}

/* -------------------------------------------------------------- --update */

function runUpdate() {
  const dest = join(HERE, '.web-style');
  if (HERE === dest || !basename(HERE).startsWith('scripts')) {
    console.error('web-style-guard: --update only applies to a project copy under scripts/.');
    process.exit(2);
  }
  mkdirSync(dest, { recursive: true });
  const done = [];
  for (const relPath of FROZEN) {
    const from = join(ROOT, relPath);
    if (!existsSync(from)) continue;
    copyFileSync(from, join(dest, basename(relPath)));
    done.push(relPath);
  }
  console.log('web-style-guard: re-recorded canonical copies of ' + done.join(', '));
  console.log('  NOTE: the skill remains the source of truth -- run web-style-sync.mjs');
  console.log('  to push a reviewed change back into the fleet.');
}

/* ------------------------------------------------------------------ main */

function main() {
  if (ARGS.update) return runUpdate();

  checkFrozen();
  checkFileBans();
  const used = collectClassTokens();
  checkTailwindisms(used);
  checkClassResolution(used);
  checkStack();

  let total = 0;
  for (const items of findings.values()) total += items.length;

  for (const [check, items] of warnings) {
    console.warn(`\n! ${check} (${items.length} warning${items.length === 1 ? '' : 's'})`);
    for (const i of items) console.warn('    ' + i);
  }

  if (!total) {
    console.log(`web-style-guard: OK (${ROOT})`);
    return;
  }

  console.error(`\nweb-style-guard: ${total} finding${total === 1 ? '' : 's'} in ${ROOT}\n`);
  for (const [check, items] of findings) {
    console.error(`x ${check} (${items.length})`);
    for (const i of items) console.error('    ' + i);
    console.error('');
  }
  if (ARGS.warnOnly) {
    console.error('web-style-guard: --warn-only, not failing the build.');
    return;
  }
  process.exit(1);
}

main();
