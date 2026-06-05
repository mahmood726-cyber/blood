# TransfusionSynth — Restrictive vs Liberal Transfusion Synthesis

A single-file, **fully offline** dashboard that pools mortality risk ratios for
**restrictive vs liberal red-cell transfusion** strategies across an RCT set
(MINT 2024, REALITY, TITRe2, FOCUS, TRICC and others), with a DerSimonian–Laird
random-effects model, a dynamic ACS-vs-stable subgroup interaction, forest /
funnel / interaction plots, an Egger small-study check, and an optional Bayesian
(Gibbs/Metropolis) path with a Gelman–Rubin R̂ convergence read-out.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN — system fonts fall back.

## Layout

```
index.html   single-file UI (loads engine.js; the worker shares it too)
engine.js    pure deterministic statistical core — runs in Node and the browser
tests.js     Node test harness, 42 assertions
LICENSE      Apache-2.0
```

## Statistical core (`engine.js`)

| Function | What it does |
|---|---|
| `metaAnalysisDL(trials)` | DerSimonian–Laird random-effects pooling of logRR: τ² via method-of-moments `(Q−(k−1))/C`, inverse-variance 95% CI on the log scale |
| `runDiagnostics(trials, results)` | I² heterogeneity (computed against the pooled RE mean) and Egger's small-study radial-regression slope/intercept |
| `normCDF(x)` | standard normal CDF, `Φ(x)=½(1+erf(x/√2))` (for two-sided p-values) |
| `tQuantile(df)` | tabulated two-sided 0.975 Student-t critical values (R `qt()` reference, Cornish–Fisher for df>30) |
| `predictionInterval(pooled)` | 95% prediction interval on the RR scale, `t_{k−1}·√(seRE²+τ²)`; undefined (returns `null`) for k<3 |

The `metaAnalysisDL` and `runDiagnostics` bodies are extracted **verbatim** from
the dashboard's inline worker so that the page, the worker, and the Node tests
all share one source of truth. The MCMC / Gibbs sampler and the normal random
generator are intentionally **not** in `engine.js` — they are non-deterministic
and remain inline in the worker.

## Fixes applied during revival (2026-06-05)

- **Made fully offline:** removed the Google Fonts CDN `<link>`; the app now
  loads no external resource. `grep -nE 'https?://' index.html` returns nothing.
- **Single source of truth:** extracted the deterministic core into `engine.js`;
  the page loads it via `<script src="engine.js">`, and the worker (a separate
  realm) pulls in the same file via `importScripts` using an absolute URL the
  page injects before building the worker Blob. The inline duplicates of
  `metaAnalysisDL` and `runDiagnostics` were deleted.
- **Added `tests.js`** (42 assertions, all passing), with every expected value
  hand-derived independently of the engine.
- **Added helpers** (`normCDF`, `tQuantile`, `predictionInterval`) so the pooled
  logRR can yield a correct two-sided p-value and a `t_{k−1}` prediction
  interval; `normCDF(0)` is exactly 0.5.
- **Renamed** `Blood.html` → `index.html` and added the Pages scaffold
  (`.nojekyll`, `.gitignore`, this README).

No methodology was changed. One behaviour worth noting (left **as shipped**,
documented not "fixed"): `runDiagnostics` computes the I² and Egger statistics
against the pooled **random-effects** mean rather than the fixed-effect mean.
This is a defensible-but-non-standard convention; it was preserved for
continuity and is pinned by an explicit test so it cannot drift silently.

## Tests

```
node tests.js
# 42 passed, 0 failed
```

Checks include normal-CDF reference points (Φ(0)=0.5 exactly, Φ(1.96)=0.975),
t-quantile reference values, an empty guard, a single-trial passthrough (k=1,
τ²=0), a two-identical-trial case (τ²=0, I²=0, seRE=√(1/20)), a hand-computed
homogeneous pair (Q<df ⇒ τ²=0, pooled RR≈0.838), a hand-computed heterogeneous
pair (τ²≈0.395, I²≈94.3%), a three-study prediction interval, and a hand-derived
Egger intercept/slope.

## Caveats

DerSimonian–Laird under-estimates τ² for small *k* (REML / Paule–Mandel are
preferred for k<10); this dashboard preserves the original method for continuity
and reports τ² and I² alongside every estimate. The bundled trial vector encodes
published logRR/variance approximations and the patient-facing recommendation
panel is illustrative — treat pooled estimates as hypothesis-generating, not a
clinical decision rule. Apache-2.0 licensed.
