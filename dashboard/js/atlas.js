/* atlas.js: browse specific / non-specific MTFs for every group, for either dataset (DepMap or TCGA). */
const Atlas = (() => {
  let dataset = "depmap", manifest, mtDesc = null, typeDesc = null, linesIdx = null, curDiv, curOpts = [], curSpec = [], curNons = [], combo;
  const manifests = {}, byGroup = {};          // caches; byGroup keyed by "dataset:div"
  const DS = () => U.DATASETS[dataset];
  const dp = f => DS().prefix + f;             // data-path for the current dataset
  const isLine = () => curDiv === "line";
  const bg = () => byGroup[`${dataset}:${curDiv}`];
  // code -> full-name map for "coded" levels (DepMap model type, TCGA tumor type), else null
  const codeDesc = () => (dataset === "depmap" && curDiv === "modeltype") ? mtDesc
    : (dataset === "tcga" && curDiv === "type") ? typeDesc : null;

  async function loadManifest() {
    manifest = manifests[dataset] ||= await DataLoader.loadJSON(dp("manifest.json"));
    if (dataset === "depmap" && !mtDesc) mtDesc = await DataLoader.loadJSON("data/modeltype_desc.json");
    if (dataset === "tcga" && !typeDesc) typeDesc = await DataLoader.loadJSON("data/tcga/type_desc.json");
  }

  async function ensure(div) {
    if (div === "line") { linesIdx ||= await DataLoader.loadJSON("data/lines_index.json"); return; }
    const ck = `${dataset}:${div}`;
    if (byGroup[ck]) return;
    const { rows } = await DataLoader.loadTSV(dp(`mtfs_${div}.tsv`));
    const m = {};
    for (const r of rows) (m[r.group] ||= { spec: [], nons: [] })[r.category === "specific" ? "spec" : "nons"]
      .push([r.tf, r.jsd_rank, r.group_expr_log2tpm]);
    for (const g in m) { m[g].spec.sort((a, b) => a[1] - b[1]); m[g].nons.sort((a, b) => a[1] - b[1]); }
    byGroup[ck] = m;
  }

  function optionsFor() {                                    // {key, label, search} for the combobox
    if (isLine()) return linesIdx.map(l => ({
      key: l.a, label: `${l.n} · ${l.s || "–"}`, search: `${l.n} ${l.a} ${l.s || ""}`.toLowerCase(),
      name: l.n, sub: l.s }));
    const groups = Object.entries(manifest.divisions[curDiv].groups).sort((a, b) => b[1] - a[1]);
    const cd = codeDesc();
    if (cd) return groups.map(([g, n]) => ({   // show + search the code's expanded name (e.g. SKCM · Skin Cutaneous Melanoma)
      key: g, label: cd[g] ? `${g} · ${cd[g]}` : g, search: `${g} ${cd[g] || ""}`.toLowerCase(), n }));
    return groups.map(([g, n]) => ({ key: g, label: g, search: g.toLowerCase(), n }));
  }

  function fillPicker() {
    curOpts = optionsFor();
    combo.setOptions(curOpts);
    // land on the most illustrative group (the one with the most specific MTFs); line level: first line
    let def = curOpts[0];
    if (!isLine() && bg()) def = curOpts.reduce((best, o) =>
      (bg()[o.key]?.spec.length || 0) > (bg()[best.key]?.spec.length || 0) ? o : best, curOpts[0]);
    combo.setValue(def ? def.key : null);
    if (def) selectVal(def.key);
  }

  const setDesc = html => U.el("atlas-desc").innerHTML = html;
  const row = ([tf, rank, expr]) =>
    `<div class="mtf" title="${tf}: rank #${rank} by CaCTS specificity · mean expression ${expr.toFixed(1)} log₂">
      <div class="g"><span class="tf-sym">${tf}</span></div><div class="meta">rank <b>#${rank}</b> · ${expr.toFixed(1)}</div></div>`;

  function renderLists(spec, nons) {
    curSpec = spec; curNons = nons;            // keep for TSV download
    U.el("atlas-hint").textContent = `${spec.length} specific · ${nons.length} non-specific`;
    U.el("atlas-spec").innerHTML = spec.length ? spec.map(row).join("")
      : `<div class="empty">No factor clears both the FDR &lt; 0.10 and ≥ 1 TPM cutoffs here.</div>`;
    U.el("atlas-nons").innerHTML = nons.length ? nons.map(row).join("") : `<div class="empty">None.</div>`;
  }

  async function selectVal(key) {
    const o = curOpts.find(x => x.key === key);
    if (!o) return;
    if (isLine()) {
      const [d, tfNames] = await Promise.all([
        DataLoader.loadJSON(`data/lines/${o.key}.json`), DataLoader.loadJSON("data/tf_names.json")]);
      const cl = U.classify(tfNames.map((tf, i) => ({ tf, cacts: d.c[i], expr: d.e[i] })));
      renderLists(cl.filter(r => r.cat === "specific").map(r => [r.tf, r.rank, r.expr]),
                  cl.filter(r => r.cat === "non_specific").map(r => [r.tf, r.rank, r.expr]));
      setDesc(`<b>${U.esc(o.name)}</b> &middot; ${U.esc(o.sub || "–")} &middot; <span class="mono">${o.key}</span>`);
    } else {
      const d = bg()[o.key], unit = DS().unit;
      renderLists(d ? d.spec : [], d ? d.nons : []);
      let ds = `${o.n.toLocaleString()} ${unit}${o.n === 1 ? "" : "s"} in this group`;
      const cd = codeDesc();
      if (cd && cd[o.key]) ds = `<b>${U.esc(o.key)}</b> &rarr; ${U.esc(cd[o.key])} &middot; ${ds}`;
      setDesc(ds);
    }
  }

  function download() {
    const o = curOpts.find(x => x.key === combo.getKey()) || {};
    const g = o.name || o.key || "group";
    const rows = curSpec.map(x => ({ tf: x[0], rank: x[1], expr: x[2], cat: "specific" }))
      .concat(curNons.map(x => ({ tf: x[0], rank: x[1], expr: x[2], cat: "non_specific" })));
    const cols = [{ label: "tf", key: "tf" }, { label: "cacts_rank", key: "rank" },
      { label: "expr_log2", get: r => (+r.expr).toFixed(4) }, { label: "class", key: "cat" }];
    const safe = g.replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "");
    U.downloadTSV(`pycacts_MTFs_${dataset}_${curDiv}_${safe}.tsv`, cols, rows);
  }

  async function pick(div) { curDiv = div; await ensure(div); fillPicker(); }

  function buildLevelSeg() {
    U.el("atlas-div").innerHTML = DS().divs.map(([k, lab]) =>
      `<button role="tab" data-div="${k}" aria-selected="${k === curDiv}" title="group by ${lab.toLowerCase()}">${lab}</button>`).join("");
  }

  async function setDataset(ds) {
    if (ds === dataset) return;
    dataset = ds; curDiv = DS().divs[0][0];
    await loadManifest();
    buildLevelSeg();
    await ensure(curDiv); fillPicker();
  }

  async function init() {
    dataset = "depmap"; curDiv = DS().divs[0][0];
    await loadManifest();
    const dseg = U.el("atlas-ds");
    dseg.innerHTML = U.DSLIST.map(([k, lab]) =>
      `<button role="tab" data-ds="${k}" aria-selected="${k === dataset}" title="show ${lab} groups">${lab}</button>`).join("");
    dseg.addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...dseg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await setDataset(b.dataset.ds);
    });
    buildLevelSeg();
    U.el("atlas-div").addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...e.currentTarget.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await pick(b.dataset.div);
    });
    combo = Combo.make(U.el("atlas-group"), key => selectVal(key));
    U.el("atlas-dl").addEventListener("click", download);
    await ensure(curDiv); fillPicker();
  }
  return { init };
})();
