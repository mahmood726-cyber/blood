/*
 * TransfusionSynth engine — pure meta-analysis core for the restrictive-vs-
 * liberal transfusion evidence-synthesis dashboard.
 *
 * Extracted verbatim from the dashboard's inline worker script so the
 * deterministic statistical core is a single source of truth, importable under
 * Node for testing. Browser: functions are globals (plain declarations).
 * Node: module.exports.
 *
 * Method: DerSimonian-Laird random-effects pooling of logRR (log-scale), with
 * the method-of-moments tau^2 = max(0, (Q-(k-1))/C), I^2, an inverse-variance
 * 95% confidence interval, and Egger's small-study radial-regression intercept.
 * The MCMC / stochastic sampling paths in the dashboard are intentionally NOT
 * extracted here (they are non-deterministic and untestable as pure functions).
 *
 * Helpers (normCDF / tQuantile / predictionInterval) are added so the pooled
 * logRR can be turned into a two-sided p-value and a t_{k-1} prediction
 * interval; they are not used by the legacy inline UI but complete the core.
 */

// --- DerSimonian-Laird random-effects pooling on logRR (VERBATIM from worker) ---
function metaAnalysisDL(groupTrials) {
    if(groupTrials.length === 0) return { mean:0, se:0, lo:0, hi:0 };

    let num=0, den=0;
    // FE Weights for Q
    groupTrials.forEach(t => {
      const w = 1/t.varLogRR;
      num += t.logRR * w;
      den += w;
    });
    const muFE = num/den;

    // Q & Tau2
    const Q = groupTrials.reduce((acc, t) => acc + (1/t.varLogRR) * Math.pow(t.logRR - muFE, 2), 0);
    const df = Math.max(1, groupTrials.length - 1);
    const C = den - groupTrials.reduce((acc, t) => acc + Math.pow(1/t.varLogRR, 2), 0) / den;
    const tau2 = Math.max(0, (Q - df) / C);

    // RE Weights
    num=0; den=0;
    groupTrials.forEach(t => {
      const w = 1 / (t.varLogRR + tau2);
      num += t.logRR * w;
      den += w;
    });
    const muRE = num/den;
    const seRE = Math.sqrt(1/den);

    return {
      mean: muRE,
      se: seRE,
      lo: muRE - 1.96*seRE,
      hi: muRE + 1.96*seRE,
      tau2: tau2,
      Q: Q,
      k: groupTrials.length
    };
}

// --- Heterogeneity (I^2) and Egger's small-study test (VERBATIM from worker) ---
function runDiagnostics(trials, results) {
    // I2 Calculation (approx)
    let Q = 0;
    const mu = results.logRR.mean;
    trials.forEach(t => {
       Q += (1/t.varLogRR) * Math.pow(t.logRR - mu, 2);
    });
    const df = Math.max(1, trials.length - 1);
    const I2 = Math.max(0, (Q - df) / Q) * 100;

    // Egger's Test (Calculated)
    const x = [], y = [];
    trials.forEach(t => {
       const se = Math.sqrt(t.varLogRR);
       x.push(1/se);
       y.push(t.logRR/se);
    });
    const n = trials.length;
    const sumX = x.reduce((a,b)=>a+b,0), sumY = y.reduce((a,b)=>a+b,0);
    const sumXY = x.reduce((a,b,i)=>a+b*y[i],0), sumXX = x.reduce((a,b)=>a+b*b,0);
    const slope = (n*sumXY - sumX*sumY) / (n*sumXX - sumX*sumX);
    const intercept = (sumY - slope*sumX) / n;
    const pVal = Math.abs(intercept) > 1.5 ? 0.04 : 0.25; // Approximate p-val mapping for UI

    return {
      I2: I2.toFixed(1) + "%",
      I2num: I2,
      egger: pVal < 0.05 ? `p=${pVal} (Bias)` : `p=${pVal} (No Bias)`,
      eggerP: pVal,
      eggerIntercept: intercept,
      eggerSlope: slope
    };
}

// --- Helpers (added during 2026-06 revival) ---

// Standard normal CDF via the Abramowitz & Stegun 7.1.26 erf approximation.
// Phi(x) = 0.5*(1 + erf(x/sqrt(2))). Phi(0) = 0.5 exactly. |abs error| < 1.5e-7.
function erf(x) {
    const sign = x < 0 ? -1 : 1;
    const ax = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * ax);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax);
    return sign * y;
}
function normCDF(x) {
    return 0.5 * (1 + erf(x / Math.sqrt(2)));
}

// Inverse Student-t at the two-sided 0.975 level. Exact table for df 1..30
// (matches R qt() to 3 dp); Cornish-Fisher expansion for df > 30.
function tQuantile(df) {
    const T = {
        1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
        8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
        15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
        22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056, 27: 2.052, 28: 2.048,
        29: 2.045, 30: 2.042
    };
    if (df <= 0) return Infinity;
    const d = Math.round(df);
    if (d <= 30) return T[d];
    const z = 1.959963985;             // qnorm(0.975)
    const z2 = z * z, z3 = z2 * z, z5 = z3 * z2;
    return z
        + (z3 + z) / (4 * df)
        + (5 * z5 + 16 * z3 + 3 * z) / (96 * df * df);
}

// 95% prediction interval on the RR scale from a DL pooled result.
// Uses t_{k-1} * sqrt(seRE^2 + tau2). Undefined (returns null) for k<3.
function predictionInterval(pooled) {
    if (!pooled || pooled.k < 3) return null;
    const piSE = Math.sqrt(pooled.se * pooled.se + pooled.tau2);
    const tCrit = tQuantile(pooled.k - 1);
    return {
        lo: Math.exp(pooled.mean - tCrit * piSE),
        hi: Math.exp(pooled.mean + tCrit * piSE)
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { metaAnalysisDL, runDiagnostics, normCDF, erf, tQuantile, predictionInterval };
}
