/*
 * Node tests for the TransfusionSynth engine. Run: node tests.js
 * Every expected value is hand-computed independently (see the derivation
 * comments) — NOT a re-run of the engine against itself.
 */
const { metaAnalysisDL, runDiagnostics, normCDF, tQuantile, predictionInterval } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

const L = Math.log;

// --- normCDF reference points (hand: Phi(0)=0.5 exactly; Phi(1.96)=0.975) ---
ok('normCDF(0) == 0.5 (0 would be a BUG)', close(normCDF(0), 0.5, 1e-9), 'got ' + normCDF(0));
ok('normCDF(1.96) ~ 0.975', close(normCDF(1.96), 0.975, 2e-3), 'got ' + normCDF(1.96));
ok('normCDF(-1.96) ~ 0.025', close(normCDF(-1.96), 0.025, 2e-3), 'got ' + normCDF(-1.96));
ok('normCDF symmetric: Phi(x)+Phi(-x)=1', close(normCDF(1.3) + normCDF(-1.3), 1, 1e-9));

// --- tQuantile (two-sided 0.975 t critical values, R qt() reference) ---
ok('tQuantile df=2 -> 4.303', tQuantile(2) === 4.303);
ok('tQuantile df=1 -> 12.706', tQuantile(1) === 12.706);
ok('tQuantile df=10 -> 2.228', tQuantile(10) === 2.228);
ok('tQuantile large df -> ~1.96', close(tQuantile(100000), 1.96, 1e-3), 'got ' + tQuantile(100000));

// --- empty guard (hand: defined return, mean 0) ---
const empty = metaAnalysisDL([]);
ok('empty -> mean 0', empty.mean === 0);
ok('empty -> se 0', empty.se === 0);

// --- single trial (k=1): df guard, tau2 clamps to 0, mean == the study ---
// y=ln(0.8)=-0.223144, var=0.1, seRE=sqrt(0.1)=0.316228
const single = metaAnalysisDL([{ logRR: L(0.8), varLogRR: 0.1 }]);
ok('single: k=1', single.k === 1);
ok('single: mean == ln(0.8)', close(single.mean, L(0.8), 1e-9), 'got ' + single.mean);
ok('single: tau2 == 0', single.tau2 === 0);
ok('single: seRE == sqrt(0.1)', close(single.se, Math.sqrt(0.1), 1e-9), 'got ' + single.se);

// --- two identical trials: pooled == study, tau2=0, I2=0 ---
// y=ln(0.8), var=0.1 each. Q=0, so tau2=0; RE weights 10+10=20; seRE=sqrt(1/20)
const ident = metaAnalysisDL([{ logRR: L(0.8), varLogRR: 0.1 }, { logRR: L(0.8), varLogRR: 0.1 }]);
ok('identical: mean == ln(0.8)', close(ident.mean, L(0.8), 1e-9), 'got ' + ident.mean);
ok('identical: tau2 == 0', close(ident.tau2, 0, 1e-12), 'got ' + ident.tau2);
ok('identical: Q == 0', close(ident.Q, 0, 1e-12), 'got ' + ident.Q);
ok('identical: seRE == sqrt(1/20)', close(ident.se, Math.sqrt(1 / 20), 1e-9), 'got ' + ident.se);
const identDiag = runDiagnostics(
    [{ logRR: L(0.8), varLogRR: 0.1 }, { logRR: L(0.8), varLogRR: 0.1 }],
    { logRR: { mean: ident.mean } });
ok('identical: I2 == 0', close(identDiag.I2num, 0, 1e-9), 'got ' + identDiag.I2num);

// --- two-study hand-worked DL pooling (homogeneous: Q<df -> tau2=0) ---
// t1: logRR=ln(0.7), var=0.1 (w=10); t2: logRR=ln(1.2), var=0.2 (w=5)
// muFE=-0.177009; Q=0.968391 (< df=1) -> tau2=0; muRE=muFE; seRE=0.258199;
// RR=exp(muRE)=0.837772; I2=0
const homo = metaAnalysisDL([{ logRR: L(0.7), varLogRR: 0.1 }, { logRR: L(1.2), varLogRR: 0.2 }]);
ok('homo: Q ~ 0.968391', close(homo.Q, 0.968391, 1e-5), 'got ' + homo.Q);
ok('homo: tau2 == 0 (Q<df)', homo.tau2 === 0);
ok('homo: muRE ~ -0.177009', close(homo.mean, -0.177009, 1e-5), 'got ' + homo.mean);
ok('homo: seRE ~ 0.258199', close(homo.se, 0.258199, 1e-5), 'got ' + homo.se);
ok('homo: RR ~ 0.837772', close(Math.exp(homo.mean), 0.837772, 1e-5), 'got ' + Math.exp(homo.mean));
ok('homo: CI lo ~ -0.683079', close(homo.lo, -0.683079, 1e-5), 'got ' + homo.lo);
ok('homo: CI hi ~ 0.329060', close(homo.hi, 0.329060, 1e-5), 'got ' + homo.hi);
const homoDiag = runDiagnostics(
    [{ logRR: L(0.7), varLogRR: 0.1 }, { logRR: L(1.2), varLogRR: 0.2 }],
    { logRR: { mean: homo.mean } });
ok('homo: I2 == 0 (Q<df)', close(homoDiag.I2num, 0, 1e-9), 'got ' + homoDiag.I2num);

// --- two-study HETEROGENEOUS DL pooling (tau2>0, I2>0) ---
// t1: logRR=ln(0.6), var=0.02 (w=50); t2: logRR=ln(1.5), var=0.03 (w=33.33)
// muFE=-0.144309; Q=16.791774; C=40 -> tau2=0.394794; muRE=-0.058137;
// seRE=0.458113; RR=0.943521; I2=94.0447%
const het = metaAnalysisDL([{ logRR: L(0.6), varLogRR: 0.02 }, { logRR: L(1.5), varLogRR: 0.03 }]);
ok('het: Q ~ 16.791774', close(het.Q, 16.791774, 1e-4), 'got ' + het.Q);
ok('het: tau2 ~ 0.394794', close(het.tau2, 0.394794, 1e-5), 'got ' + het.tau2);
ok('het: muRE ~ -0.058137', close(het.mean, -0.058137, 1e-5), 'got ' + het.mean);
ok('het: seRE ~ 0.458113', close(het.se, 0.458113, 1e-5), 'got ' + het.se);
ok('het: RR ~ 0.943521', close(Math.exp(het.mean), 0.943521, 1e-5), 'got ' + Math.exp(het.mean));
const hetDiag = runDiagnostics(
    [{ logRR: L(0.6), varLogRR: 0.02 }, { logRR: L(1.5), varLogRR: 0.03 }],
    { logRR: { mean: het.mean } });
// NOTE: runDiagnostics computes Q against the pooled RE mean (results.logRR.mean),
// not the FE mean, so I2 here is derived from muRE=-0.058137:
// Q_RE = (1/0.02)(ln0.6-muRE)^2 + (1/0.03)(ln1.5-muRE)^2 = 17.41058 -> I2=94.25637%
ok('het: I2 ~ 94.25637% (Q vs RE mean)', close(hetDiag.I2num, 94.25637, 1e-3), 'got ' + hetDiag.I2num);

// --- three-study case + prediction interval (k=3 -> t_{k-1}=t_2=4.303) ---
// studies: ln(0.6)/0.02, ln(1.5)/0.03, ln(0.9)/0.05
// muRE=-0.075339; seRE=0.29575; tau2=0.229652; piSE=0.563134
// PI(exp): [0.082208, 10.462792]
const tri = [
    { logRR: L(0.6), varLogRR: 0.02 },
    { logRR: L(1.5), varLogRR: 0.03 },
    { logRR: L(0.9), varLogRR: 0.05 }
];
const triPool = metaAnalysisDL(tri);
ok('tri: muRE ~ -0.075339', close(triPool.mean, -0.075339, 1e-5), 'got ' + triPool.mean);
ok('tri: seRE ~ 0.29575', close(triPool.se, 0.29575, 1e-4), 'got ' + triPool.se);
ok('tri: tau2 ~ 0.229652', close(triPool.tau2, 0.229652, 1e-5), 'got ' + triPool.tau2);
const triPI = predictionInterval(triPool);
ok('tri: PI lo ~ 0.082208', close(triPI.lo, 0.082208, 1e-4), 'got ' + triPI.lo);
ok('tri: PI hi ~ 10.462792', close(triPI.hi, 10.462792, 2e-3), 'got ' + triPI.hi);

// --- prediction interval undefined for k<3 ---
ok('PI undefined for k=2', predictionInterval(het) === null);
ok('PI undefined for k=1', predictionInterval(single) === null);

// --- Egger intercept (hand-derived on the 3-study set) ---
// 1/se and logRR/se points -> slope=-1.206888, intercept=6.385671
const eg = runDiagnostics(tri, { logRR: { mean: triPool.mean } });
ok('Egger intercept ~ 6.385671', close(eg.eggerIntercept, 6.385671, 1e-4), 'got ' + eg.eggerIntercept);
ok('Egger slope ~ -1.206888', close(eg.eggerSlope, -1.206888, 1e-4), 'got ' + eg.eggerSlope);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
