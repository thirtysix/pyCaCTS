/* scores.js, every TF's CaCTS specificity score + expression rank for a chosen group, sortable.
   Class (specific / non-specific) recomputed as in the pipeline. Works at all five levels incl. single cell line. */
const Scores = (() => {
  let manifest, mtDesc = null, linesIdx = null, curDiv = "subtype", curOpts = [];
  const isLine = () => curDiv === "line";
  const state = { rows: [], sort: { col: "rank", dir: "asc" }, filter: "", fdrMax: 1 };

  async function ensure(div) { if (div === "line") linesIdx ||= await DataLoader.loadJSON("data/lines_index.json"); }

  function optionsFor() {
    if (isLine()) return linesIdx.map(l => ({ val: `${l.n} · ${l.a}`, hint: l.s, key: l.a, name: l.n, sub: l.s }));
    return Object.entries(manifest.divisions[curDiv].groups).sort((a, b) => b[1] - a[1])
      .map(([g, n]) => ({ val: g, hint: `n=${n}`, key: g, n }));
  }

  async function recsFor(o) {
    if (isLine()) {
      const [d, tfNames] = await Promise.all([
        DataLoader.loadJSON(`data/lines/${o.key}.json`), DataLoader.loadJSON("data/tf_names.json")]);
      return tfNames.map((tf, i) => ({ tf, cacts: d.c[i], expr: d.e[i] }));
    }
    const [S, E] = await Promise.all([DataLoader.loadTSV(`data/scores_${curDiv}.tsv`), DataLoader.loadTSV(`data/expr_${curDiv}.tsv`)]);
    const tfCol = S.columns[0], ex = {};
    for (const r of E.rows) ex[r[E.columns[0]]] = r[o.key];
    return S.rows.map(r => ({ tf: r[tfCol], cacts: r[o.key], expr: ex[r[tfCol]] }));
  }

  function setDesc(o) {
    let html;
    if (isLine()) html = `<b>${U.esc(o.name)}</b> &middot; ${U.esc(o.sub || "–")} &middot; <span class="mono">${o.key}</span>`;
    else {
      html = `${o.n.toLocaleString()} cell line${o.n === 1 ? "" : "s"}`;
      if (curDiv === "modeltype" && mtDesc && mtDesc[o.key]) html = `<b>${U.esc(o.key)}</b> &rarr; ${U.esc(mtDesc[o.key])} &middot; ${html}`;
    }
    U.el("scores-desc").innerHTML = html;
  }

  async function load(o) {
    const recs = await recsFor(o);
    const rows = U.classify(recs), fdr = U.empiricalFDR(recs);
    rows.forEach(r => r.fdr = fdr[r.tf]);
    state.rows = rows; setDesc(o); render();
  }

  const CLS = {
    specific: '<span class="cls spec" title="top-5% by CaCTS score AND top-5% by expression">specific</span>',
    non_specific: '<span class="cls nons" title="top-5% expression but not top-5% specificity, a candidate ubiquitous MTF">non-specific</span>',
    "": '<span class="cls none" title="not a called MTF in this group">–</span>',
  };
  const LSIG = Math.log(0.10);                             // ln-FDR significance line (FDR < 0.10)
  const fmtFDR = lf => {                                    // lf = ln(FDR)
    if (lf == null) return "–";
    const l10 = lf / Math.LN10;                            // log10(FDR)
    if (l10 >= -3) return Math.pow(10, l10).toFixed(3);    // ≥ 1e-3
    const e = Math.floor(l10), m = Math.pow(10, l10 - e);  // mantissa ∈ [1,10); representable below float64 range
    return `${m.toFixed(1)}e${e}`;                         // e.g. 2.3e-340
  };

  function render() {
    const { col, dir } = state.sort, f = state.filter.toUpperCase();
    let rows = state.rows;
    if (f) rows = rows.filter(r => r.tf.toUpperCase().includes(f));
    if (state.fdrMax < 1) rows = rows.filter(r => r.fdr != null && r.fdr <= Math.log(state.fdrMax));
    rows = [...rows].sort((a, b) => {
      const va = a[col], vb = b[col], cmp = typeof va === "string" ? va.localeCompare(vb) : (va - vb) || 0;
      return dir === "asc" ? cmp : -cmp;
    });
    U.el("scores-body").innerHTML = rows.map(r =>
      `<tr><td class="num mono">${r.rank}</td>
        <td><span class="tf-sym">${r.tf}</span></td>
        <td class="num mono">${r.cacts.toFixed(4)}</td>
        <td class="num mono">${r.expr.toFixed(1)}</td>
        <td class="num mono">${r.exprrank}</td>
        <td class="num mono${r.fdr != null && r.fdr < LSIG ? " sig" : ""}" title="empirical-null FDR = ${r.fdr == null ? "n/a" : fmtFDR(r.fdr)}">${fmtFDR(r.fdr)}</td>
        <td>${CLS[r.cat]}</td></tr>`).join("");
    document.querySelectorAll("#scores-table th[data-sort]").forEach(th => {
      th.classList.remove("sorted-asc", "sorted-desc");
      if (th.dataset.sort === col) th.classList.add(dir === "asc" ? "sorted-asc" : "sorted-desc");
    });
    const nSpec = state.rows.filter(r => r.cat === "specific").length;
    const nSig = state.rows.filter(r => r.fdr != null && r.fdr < LSIG).length;
    const flt = (f || state.fdrMax < 1) ? " (filtered)" : "";
    U.el("scores-foot").innerHTML = `${rows.length.toLocaleString()} TFs shown${flt} · ${nSpec} specific · ${nSig} significant at empirical-null FDR&lt;0.10 · lower CaCTS score / FDR = more specific.`;
  }

  function fillPicker() {
    curOpts = optionsFor();
    U.fillSelect(U.el("scores-group"), curOpts);
    const def = curOpts[0];                     // largest group (group levels) / first line (line level)
    U.el("scores-group").value = def ? def.val : "";
    if (def) load(def);
  }

  function selectVal(val) {
    const o = curOpts.find(x => x.val === val) || curOpts.find(x => x.val.toUpperCase() === val.toUpperCase());
    if (o) load(o);
  }
  async function pick(div) { curDiv = div; await ensure(div); fillPicker(); }

  async function init() {
    [manifest, mtDesc] = await Promise.all([
      DataLoader.loadJSON("data/manifest.json"), DataLoader.loadJSON("data/modeltype_desc.json")]);
    const seg = U.el("scores-div");
    seg.innerHTML = U.DIVS.map(([k, lab]) =>
      `<button role="tab" data-div="${k}" aria-selected="${k === curDiv}" title="score TFs specifically within each ${lab.toLowerCase()}">${lab}</button>`).join("");
    seg.addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...seg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await pick(b.dataset.div);
    });
    U.el("scores-group").addEventListener("change", e => selectVal(e.target.value));
    U.el("scores-filter").addEventListener("input", e => { state.filter = e.target.value.trim(); render(); });
    const fseg = U.el("scores-fdr");
    fseg.addEventListener("click", e => {
      const b = e.target.closest("button"); if (!b) return;
      [...fseg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      state.fdrMax = parseFloat(b.dataset.fdr); render();
    });
    document.querySelectorAll("#scores-table th[data-sort]").forEach(th => th.addEventListener("click", () => {
      const c = th.dataset.sort;
      if (state.sort.col === c) state.sort.dir = state.sort.dir === "asc" ? "desc" : "asc";
      else state.sort = { col: c, dir: th.dataset.num ? (["rank", "exprrank", "fdr"].includes(c) ? "asc" : "desc") : "asc" };
      render();
    }));
    await ensure(curDiv); fillPicker();
  }
  return { init };
})();
