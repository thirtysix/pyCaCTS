/* overview.js — hero stats, benchmark, and a guide to the tabs. */
const Overview = (() => {
  async function init() {
    const [meta, bench] = await Promise.all([
      DataLoader.loadJSON("data/meta.json"),
      DataLoader.loadJSON("data/benchmark.json"),
    ]);
    // sidebar snapshot
    U.el("snap").innerHTML = [
      ["≡ R", "identical to original"],
      [meta.speedup, "faster than R"],
      [meta.n_tfs, "TFs scored"],
      [meta.n_lines, "CCLE cell lines"],
    ].map(([b, l]) => `<div class="s"><b>${b}</b> ${l}</div>`).join("");

    U.el("ov-stats").innerHTML = [
      ["≡ R", "validated", "identical to the original R<br>(max |Δ| &lt; 1e-15)"],
      [meta.speedup, "good", "faster than the R implementation"],
      [meta.n_tfs, "", "transcription factors (CaCTS catalogue)"],
      [meta.n_lines, "", "DepMap / CCLE cell lines"],
    ].map(([k, cls, l]) => `<div class="stat"><div class="k ${cls}">${k}</div><div class="l">${l}</div></div>`).join("");

    // Log scale: a linear axis can't render a >3,000× gap (the fast bar would be an invisible sliver).
    const lg = v => Math.log10(Math.max(v, 1));
    const maxLg = Math.max(...bench.rows.flatMap(r => [lg(r.py_ms), lg(r.r_ms)])) || 1;
    const w = v => Math.max(2, (lg(v) / maxLg) * 100);
    U.el("ov-bench").innerHTML = `<div class="bench">` + bench.rows.map(r =>
      `<div class="brow">
        <div class="bhead"><span class="blbl">${r.input}</span><span class="bspd" title="measured wall-clock speedup">${r.speedup} <span>faster</span></span></div>
        <div class="btrack" title="pyCaCTS: ${r.py_ms} ms"><div class="bfill py" style="width:${w(r.py_ms)}%"></div><span class="bval">pyCaCTS · ${r.py_ms} ms</span></div>
        <div class="btrack" title="original R: ${r.r_ms.toLocaleString()} ms"><div class="bfill r" style="width:${w(r.r_ms)}%"></div><span class="bval">R · ${r.r_ms.toLocaleString()} ms</span></div>
      </div>`).join("") +
      `<div class="bnote">Bars are on a <b>log₁₀(time)</b> scale — on a linear axis the pyCaCTS bar would be an invisible sliver at this ratio. The multiplier at right is the true measured speedup.</div></div>`;

    U.el("ov-guide").innerHTML = [
      ["▦", "atlas", "MTF atlas", "specific &amp; non-specific MTFs for every group, at four resolutions"],
      ["≣", "scores", "TF scores", "every TF's specificity score for a group, sortable"],
      ["⌕", "finder", "TF finder", "every group where a given TF is a specific master regulator"],
      ["ⓘ", "about", "About &amp; methods", "what CaCTS computes, how an MTF is called, validation"],
    ].map(([i, tab, t, s]) => `<a href="#${tab}"><span class="gi">${i}</span><span class="gt"><b>${t}</b><small>${s}</small></span></a>`).join("");
  }
  return { init };
})();
