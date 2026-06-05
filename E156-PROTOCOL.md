# E156-PROTOCOL — TransfusionSynth (Restrictive vs Liberal Transfusion Synthesis)

- **Project:** blood (GitHub repo `blood`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `Blood.html` dump, "TransfusionSynth V1.2: Validated Engine")
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>`; the app now
  loads no external resource (system fonts fall back). `grep -nE 'https?://'`
  on the HTML returns nothing.
- Extracted the deterministic statistical core (`metaAnalysisDL`,
  `runDiagnostics`) **verbatim** into a pure `engine.js` (single source of
  truth). The page loads it via `<script src="engine.js">`; the worker, a
  separate realm, pulls in the same file via `importScripts` of a page-injected
  absolute URL. The inline duplicates were deleted.
- Added `tests.js` (42 assertions, all passing), each expected value hand-derived
  independently of the engine.
- Added helper functions (`normCDF`, `tQuantile`, `predictionInterval`); verified
  `normCDF(0)=0.5` exactly and a `t_{k−1}` prediction interval.
- Documented (not "fixed") that I²/Egger are computed against the pooled
  random-effects mean — a non-standard but faithful-to-original convention, now
  pinned by a test.
- Added Pages scaffold (`.nojekyll`, `.gitignore`, README); renamed
  `Blood.html` → `index.html`. No statistical method was changed.

## Body (E156 draft — CURRENT BODY)

Does a restrictive red-cell transfusion threshold change mortality relative to a
liberal one, and does that effect depend on whether the patient has acute
coronary syndrome? This dashboard pools published mortality risk ratios from an
expanded randomised-trial set, including MINT, REALITY, TITRe2 and FOCUS, on the
log scale. It fits a DerSimonian–Laird random-effects model and derives
the acute-versus-stable subgroup interaction directly from the trial
differences. Diagnostics report τ², I², an Egger small-study check and, in
Bayesian mode, a Gelman–Rubin convergence statistic. Across the bundled set the
pooled risk ratio sits near unity with wide intervals and high heterogeneity,
while the acute-coronary subgroup leans toward harm from restriction. A revival
audit made the tool fully offline, extracted the deterministic core into one
tested module, and locked it behind a hand-verified forty-two-assertion suite. The honest read is that transfusion strategy is broadly
neutral but plausibly context-dependent, so this is a transparent synthesis aid,
not a bedside rule.

SUBMITTED: [ ]
