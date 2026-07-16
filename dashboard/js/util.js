/* util.js, shared helpers. The two datasets (DepMap cell lines, TCGA tumors) and their grouping levels. */
const U = (() => {
  const DIVS = [["lineage", "Lineage"], ["disease", "Primary disease"],
                ["subtype", "Subtype"], ["modeltype", "Model type"], ["line", "Cell line"]];
  // each dataset: where its data lives (prefix), the unit it groups, and its resolution levels
  const DATASETS = {
    depmap: { label: "DepMap", prefix: "data/", unit: "cell line",
              divs: [["lineage", "Lineage"], ["disease", "Primary disease"], ["subtype", "Subtype"],
                     ["modeltype", "Model type"], ["line", "Cell line"]] },
    tcga:   { label: "TCGA", prefix: "data/tcga/", unit: "sample", breadthFile: "data/tcga/breadth.json",
              divs: [["type", "Tumor type"], ["subtype", "Molecular subtype"], ["sampletype", "Sample type"]] },
  };
  const DSLIST = Object.entries(DATASETS).map(([k, v]) => [k, v.label]);
  const el = id => document.getElementById(id);
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" }[c]));
  const tags = () => "";                         // no lineage-specific badges in the general tool view

  const median = a => { const s = [...a].sort((x, y) => x - y), n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : 0; };
  // standard-normal CDF via Abramowitz-Stegun erf (|err| < 1.5e-7): no scipy in the browser
  const erf = x => { const s = x < 0 ? -1 : 1; x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return s * y; };
  const normcdf = z => 0.5 * (1 + erf(z / Math.SQRT2));

  const FDR_MAX = 0.10, EXPR_FLOOR = 1.0;                       // specific = empirical-null FDR <= 0.10 ...
  const LN_FDR_MAX = Math.log(FDR_MAX);                         // ... AND mean expr >= 1 TPM (log2(TPM+1) units)
  // Recompute the CaCTS MTF class for a set of {tf,cacts,expr} records: specific = empirical-null FDR <= 0.10
  // AND expression >= 1 TPM; non-specific (candidate ubiquitous) = top-5% expression, not specific. Ranked by score.
  const classify = recs => {
    const byScore = [...recs].sort((a, b) => a.cacts - b.cacts);
    const byExpr = [...recs].sort((a, b) => b.expr - a.expr);
    const er = {}; byExpr.forEach((r, i) => er[r.tf] = i + 1);
    const kE = Math.round(0.05 * recs.length);
    const tE = new Set(byExpr.slice(0, kE).map(r => r.tf));     // top-5% expr, for the ubiquitous (non-specific) call
    const lnFdr = empiricalFDR(recs);                          // {tf -> ln(FDR)}
    return byScore.map((r, i) => {
      const sig = lnFdr[r.tf] <= LN_FDR_MAX, floored = r.expr >= EXPR_FLOOR;
      return {
        rank: i + 1, tf: r.tf, cacts: r.cacts, expr: r.expr, exprrank: er[r.tf], fdr: lnFdr[r.tf],
        sig, floored, top5e: tE.has(r.tf),                     // the two gates behind the class + ubiquitous marker
        cat: (sig && floored) ? "specific" : (tE.has(r.tf) ? "non_specific" : ""),
      };
    });
  };

  // Download row-objects as a TSV file. cols = [{label, key} | {label, get(row)}]; rows = [{...}].
  const downloadTSV = (filename, cols, rows) => {
    const cell = v => v == null ? "" : String(v).replace(/[\t\r\n]/g, " ");
    const lines = [cols.map(c => c.label).join("\t")].concat(
      rows.map(r => cols.map(c => cell(c.get ? c.get(r) : r[c.key])).join("\t")));
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/tab-separated-values" }));
    a.download = filename; document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
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

  return { DIVS, DATASETS, DSLIST, el, esc, tags, classify, empiricalFDR, downloadTSV };
})();
