/* util.js, shared helpers. General pyCaCTS tool dashboard: the DepMap organizational levels. */
const U = (() => {
  const DIVS = [["lineage", "Lineage"], ["disease", "Primary disease"],
                ["subtype", "Subtype"], ["modeltype", "Model type"], ["line", "Cell line"]];
  const el = id => document.getElementById(id);
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" }[c]));
  const tags = () => "";                         // no lineage-specific badges in the general tool view

  const median = a => { const s = [...a].sort((x, y) => x - y), n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : 0; };
  // standard-normal CDF via Abramowitz-Stegun erf (|err| < 1.5e-7): no scipy in the browser
  const erf = x => { const s = x < 0 ? -1 : 1; x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return s * y; };
  const normcdf = z => 0.5 * (1 + erf(z / Math.SQRT2));

  // Recompute the CaCTS MTF class (specific / non-specific) for a set of {tf,cacts,expr} records,
  // exactly as the pipeline does: top-5% by score ∩ top-5% by expression. Returns rows ranked by score.
  const classify = recs => {
    const N = recs.length, k = Math.round(0.05 * N);
    const byScore = [...recs].sort((a, b) => a.cacts - b.cacts);
    const byExpr = [...recs].sort((a, b) => b.expr - a.expr);
    const er = {}; byExpr.forEach((r, i) => er[r.tf] = i + 1);
    const tS = new Set(byScore.slice(0, k).map(r => r.tf));
    const tE = new Set(byExpr.slice(0, k).map(r => r.tf));
    return byScore.map((r, i) => ({
      rank: i + 1, tf: r.tf, cacts: r.cacts, expr: r.expr, exprrank: er[r.tf],
      cat: (tS.has(r.tf) && tE.has(r.tf)) ? "specific" : (tE.has(r.tf) ? "non_specific" : ""),
    }));
  };

  // Empirical-null FDR for the CaCTS scores of one group (same method as the pipeline's Step-1 gate):
  // model the null from the non-specific (high-JSD) side, left-tail p per TF, Benjamini-Hochberg.
  // Works in NATURAL-LOG space so highly specific TFs stay distinct. The A-S normcdf suffers catastrophic
  // cancellation (1+(−1+ε)→0) for z ≲ −8, so below z=−6 we switch to the asymptotic ln Φ(z) (accurate
  // there; continuous with log(normcdf) to ~2%). Returns {tf -> ln(FDR)}, ≤ 0 (ln 0.05 ≈ -3.0).
  const lnLeftTail = z => z > -6
    ? Math.log(Math.max(normcdf(z), Number.MIN_VALUE))            // reliable central/tail region
    : -0.5 * z * z - Math.log(-z) - 0.5 * Math.log(2 * Math.PI);  // asymptotic ln Φ(z) for z ≪ 0
  const empiricalFDR = recs => {
    const x = recs.map(r => r.cacts), mu = median(x);
    const sig = 1.4826 * median(x.filter(v => v >= mu).map(v => Math.abs(v - mu))) || 1e-9;
    const lp = recs.map(r => lnLeftTail((r.cacts - mu) / sig));   // ln(left-tail p) per TF
    const order = lp.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
    const n = lp.length, out = {}; let prev = 0;                  // prev = ln of running BH min (ln 1 = 0)
    for (let k = n - 1; k >= 0; k--) { const [lv, i] = order[k];  // ln(fdr) = ln p + ln n − ln rank
      prev = Math.min(prev, lv + Math.log(n) - Math.log(k + 1)); out[recs[i].tf] = Math.min(0, prev); }
    return out;
  };

  // Populate a <select> dropdown. opts: [{val, hint}]; option value = val (what selectVal matches on).
  // Native <select> gives an obvious dropdown affordance and type-to-jump even for the ~1,450 cell lines.
  const fillSelect = (selEl, opts) =>
    selEl.innerHTML = opts.map(o => `<option value="${esc(o.val)}">${esc(o.val)}</option>`).join("");

  return { DIVS, el, esc, tags, classify, empiricalFDR, fillSelect };
})();
