/* scores.js: every TF's CaCTS specificity score for a chosen group, sortable, for either dataset.
   One table: class (specific / non-specific), the two gates behind it (top-5% score, top-5% expression),
   empirical-null FDR, gene name + family, per-group CRISPR essentiality (DepMap only), cross-group breadth,
   and out-links. DepMap: five levels incl. single cell line. TCGA: tumor type / molecular subtype / sample
   type. The current view is downloadable as TSV. */
const Scores = (() => {
  let dataset = "depmap", manifest, mtDesc = null, linesIdx = null, info = {}, tcgaBreadth = null,
      curDiv = "subtype", curOpts = [], combo;
  const manifests = {};
  const DS = () => U.DATASETS[dataset];
  const dp = f => DS().prefix + f;
  const isLine = () => curDiv === "line";
  const state = { rows: [], sort: { col: "rank", dir: "asc" }, filter: "", fdrMax: 1 };

  async function loadManifest() {
    manifest = manifests[dataset] ||= await DataLoader.loadJSON(dp("manifest.json"));
    if (dataset === "depmap" && !mtDesc) mtDesc = await DataLoader.loadJSON("data/modeltype_desc.json");
    if (dataset === "tcga" && !tcgaBreadth) tcgaBreadth = await DataLoader.loadJSON("data/tcga/breadth.json");
  }
  async function ensure(div) { if (div === "line") linesIdx ||= await DataLoader.loadJSON("data/lines_index.json"); }

  function optionsFor() {                                    // {key, label, search} for the combobox
    if (isLine()) return linesIdx.map(l => ({
      key: l.a, label: `${l.n} · ${l.s || "–"}`, search: `${l.n} ${l.a} ${l.s || ""}`.toLowerCase(),
      name: l.n, sub: l.s }));
    const groups = Object.entries(manifest.divisions[curDiv].groups).sort((a, b) => b[1] - a[1]);
    if (curDiv === "modeltype") return groups.map(([g, n]) => ({   // show + search the code's description
      key: g, label: mtDesc && mtDesc[g] ? `${g} · ${mtDesc[g]}` : g,
      search: `${g} ${mtDesc && mtDesc[g] ? mtDesc[g] : ""}`.toLowerCase(), n }));
    return groups.map(([g, n]) => ({ key: g, label: g, search: g.toLowerCase(), n }));
  }

  async function recsFor(o) {
    if (isLine()) {                                            // per-line essentiality not staged -> null
      const [d, tfNames] = await Promise.all([
        DataLoader.loadJSON(`data/lines/${o.key}.json`), DataLoader.loadJSON("data/tf_names.json")]);
      return tfNames.map((tf, i) => ({ tf, cacts: d.c[i], expr: d.e[i], ess: null }));
    }
    const [S, E, ES] = await Promise.all([
      DataLoader.loadTSV(dp(`scores_${curDiv}.tsv`)), DataLoader.loadTSV(dp(`expr_${curDiv}.tsv`)),
      dataset === "depmap" ? DataLoader.loadTSV(dp(`ess_${curDiv}.tsv`)).catch(() => null) : Promise.resolve(null)]);
    const tfCol = S.columns[0], ex = {}, es = {};
    for (const r of E.rows) ex[r[E.columns[0]]] = r[o.key];
    const hasEss = ES && ES.columns.includes(o.key);
    if (hasEss) for (const r of ES.rows) es[r[ES.columns[0]]] = r[o.key];
    return S.rows.map(r => ({ tf: r[tfCol], cacts: r[o.key], expr: ex[r[tfCol]], ess: hasEss ? es[r[tfCol]] : null }));
  }

  function setDesc(o) {
    let html;
    if (isLine()) html = `<b>${U.esc(o.name)}</b> &middot; ${U.esc(o.sub || "–")} &middot; <span class="mono">${o.key}</span>`;
    else {
      const unit = DS().unit;
      html = `${o.n.toLocaleString()} ${unit}${o.n === 1 ? "" : "s"}`;
      if (curDiv === "modeltype" && mtDesc && mtDesc[o.key]) html = `<b>${U.esc(o.key)}</b> &rarr; ${U.esc(mtDesc[o.key])} &middot; ${html}`;
    }
    U.el("scores-desc").innerHTML = html;
  }

  async function load(o) {
    const recs = await recsFor(o);
    const rows = U.classify(recs), fdr = U.empiricalFDR(recs);
    const ess = {}; recs.forEach(r => ess[r.tf] = r.ess);
    rows.forEach(r => {
      r.fdr = fdr[r.tf]; r.ess = ess[r.tf];
      const gi = info[r.tf] || {};
      r.name = gi.name || ""; r.family = gi.family || "";
      r.breadth = dataset === "tcga" ? (tcgaBreadth[r.tf] || 0) : (gi.breadth || 0);
      r.entrez = gi.entrez || ""; r.ensembl = gi.ensembl || "";
    });
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
  const chk = b => b ? '<span class="chk" title="yes">✓</span>' : '<span class="chk-no" title="no">·</span>';
  const essFmt = e => (e == null || isNaN(e)) ? "–" : (+e).toFixed(2);
  const links = r => {
    const a = (href, t, lab) => `<a href="${href}" target="_blank" rel="noopener" title="${t}">${lab}</a>`;
    const out = [];
    if (r.entrez) out.push(a(`https://www.ncbi.nlm.nih.gov/gene/${r.entrez}`, "NCBI Gene", "N"));
    out.push(a(`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(r.tf)}`, "GeneCards", "G"));
    out.push(a(`https://depmap.org/portal/gene/${encodeURIComponent(r.tf)}?tab=overview`, "DepMap", "D"));
    return `<span class="glinks">${out.join("")}</span>`;
  };

  function currentRows() {                                  // filtered + sorted view (shared by render + download)
    const { col, dir } = state.sort, f = state.filter.toUpperCase();
    let rows = state.rows;
    if (f) rows = rows.filter(r => r.tf.toUpperCase().includes(f));
    if (state.fdrMax < 1) rows = rows.filter(r => r.fdr != null && r.fdr <= Math.log(state.fdrMax));
    return [...rows].sort((a, b) => {
      let va = a[col], vb = b[col];
      if (va == null) va = col === "ess" ? Infinity : -Infinity;   // nulls sort last
      if (vb == null) vb = col === "ess" ? Infinity : -Infinity;
      const cmp = typeof va === "string" ? va.localeCompare(vb) : (va - vb) || 0;
      return dir === "asc" ? cmp : -cmp;
    });
  }

  function render() {
    const { col, dir } = state.sort, f = state.filter;
    const rows = currentRows();
    U.el("scores-body").innerHTML = rows.map(r =>
      `<tr><td class="num mono">${r.rank}</td>
        <td><span class="tf-sym" title="${U.esc(r.name || r.tf)}">${r.tf}</span></td>
        <td class="fam" title="DNA-binding-domain family">${U.esc(r.family || "–")}</td>
        <td class="num mono">${r.cacts.toFixed(4)}</td>
        <td class="ctr">${chk(r.top5s)}</td>
        <td class="num mono">${r.expr.toFixed(1)}</td>
        <td class="ctr">${chk(r.top5e)}</td>
        <td class="num mono">${r.exprrank}</td>
        <td class="num mono${r.fdr != null && r.fdr < LSIG ? " sig" : ""}" title="empirical-null FDR = ${r.fdr == null ? "n/a" : fmtFDR(r.fdr)}">${fmtFDR(r.fdr)}</td>
        <td>${CLS[r.cat]}</td>
        <td class="num mono${(r.ess != null && !isNaN(r.ess) && r.ess < -0.5) ? " ess-dep" : ""}" title="mean CRISPR (Chronos) in this group; lower = stronger dependency (DepMap only)">${essFmt(r.ess)}</td>
        <td class="num mono" title="number of groups in this dataset where this TF is a specific MTF">${r.breadth || 0}</td>
        <td class="ctr">${links(r)}</td></tr>`).join("");
    document.querySelectorAll("#scores-table th[data-sort]").forEach(th => {
      th.classList.remove("sorted-asc", "sorted-desc");
      if (th.dataset.sort === col) th.classList.add(dir === "asc" ? "sorted-asc" : "sorted-desc");
    });
    const nSpec = state.rows.filter(r => r.cat === "specific").length;
    const nSig = state.rows.filter(r => r.fdr != null && r.fdr < LSIG).length;
    const flt = (f || state.fdrMax < 1) ? " (filtered)" : "";
    U.el("scores-foot").innerHTML = `${rows.length.toLocaleString()} TFs shown${flt} · ${nSpec} specific · ${nSig} significant at empirical-null FDR&lt;0.10 · lower CaCTS score / FDR = more specific.`;
  }

  function download() {
    const o = curOpts.find(x => x.key === combo.getKey()) || {};
    const g = o.name || o.key || "group";
    const cols = [
      { label: "rank", key: "rank" }, { label: "tf", key: "tf" }, { label: "gene_name", key: "name" },
      { label: "tf_family", key: "family" }, { label: "cacts_score", get: r => r.cacts.toFixed(6) },
      { label: "top5pct_cacts", get: r => r.top5s ? 1 : 0 }, { label: "expr_log2", get: r => r.expr.toFixed(4) },
      { label: "top5pct_expr", get: r => r.top5e ? 1 : 0 }, { label: "expr_rank", key: "exprrank" },
      { label: "empirical_fdr", get: r => r.fdr == null ? "" : Math.exp(r.fdr).toExponential(3) },
      { label: "class", get: r => r.cat || "" },
      { label: "mean_chronos", get: r => (r.ess == null || isNaN(r.ess)) ? "" : (+r.ess).toFixed(3) },
      { label: "specific_in_n_groups", key: "breadth" }, { label: "entrez", key: "entrez" },
      { label: "ensembl", key: "ensembl" }];
    const safe = g.replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "");
    U.downloadTSV(`pycacts_${dataset}_${curDiv}_${safe}.tsv`, cols, currentRows());
  }

  function fillPicker() {
    curOpts = optionsFor();
    combo.setOptions(curOpts);
    const def = curOpts[0];                     // largest group (group levels) / first line (line level)
    combo.setValue(def ? def.key : null);
    if (def) load(def);
  }
  function selectVal(key) { const o = curOpts.find(x => x.key === key); if (o) load(o); }
  async function pick(div) { curDiv = div; await ensure(div); fillPicker(); }

  function buildLevelSeg() {
    U.el("scores-div").innerHTML = DS().divs.map(([k, lab]) =>
      `<button role="tab" data-div="${k}" aria-selected="${k === curDiv}" title="score TFs specifically within each ${lab.toLowerCase()}">${lab}</button>`).join("");
  }

  async function setDataset(ds) {
    if (ds === dataset) return;
    dataset = ds; curDiv = dataset === "tcga" ? "type" : "subtype";
    await loadManifest();
    buildLevelSeg();
    await ensure(curDiv); fillPicker();
  }

  async function init() {
    dataset = "depmap"; curDiv = "subtype";
    info = await DataLoader.loadJSON("data/gene_info.json");
    await loadManifest();
    const dseg = U.el("scores-ds");
    dseg.innerHTML = U.DSLIST.map(([k, lab]) =>
      `<button role="tab" data-ds="${k}" aria-selected="${k === dataset}" title="score TFs within ${lab} groups">${lab}</button>`).join("");
    dseg.addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...dseg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await setDataset(b.dataset.ds);
    });
    buildLevelSeg();
    U.el("scores-div").addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...e.currentTarget.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await pick(b.dataset.div);
    });
    combo = Combo.make(U.el("scores-group"), key => selectVal(key));
    U.el("scores-filter").addEventListener("input", e => { state.filter = e.target.value.trim(); render(); });
    U.el("scores-dl").addEventListener("click", download);
    const fseg = U.el("scores-fdr");
    fseg.addEventListener("click", e => {
      const b = e.target.closest("button"); if (!b) return;
      [...fseg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      state.fdrMax = parseFloat(b.dataset.fdr); render();
    });
    document.querySelectorAll("#scores-table th[data-sort]").forEach(th => th.addEventListener("click", () => {
      const c = th.dataset.sort;
      if (state.sort.col === c) state.sort.dir = state.sort.dir === "asc" ? "desc" : "asc";
      else state.sort = { col: c, dir: th.dataset.num ? (["rank", "exprrank", "fdr", "ess"].includes(c) ? "asc" : "desc") : "asc" };
      render();
    }));
    await ensure(curDiv); fillPicker();
  }
  return { init };
})();
