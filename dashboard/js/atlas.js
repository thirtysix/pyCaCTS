/* atlas.js, browse specific / non-specific MTFs for every group, at five resolutions incl. single cell line. */
const Atlas = (() => {
  let manifest, mtDesc = null, linesIdx = null, curDiv = "lineage", curOpts = [], curSpec = [], curNons = [], combo;
  const byGroup = {};                          // group-level: div -> {group -> {spec, nons}}
  const isLine = () => curDiv === "line";

  async function ensure(div) {
    if (div === "line") { linesIdx ||= await DataLoader.loadJSON("data/lines_index.json"); return; }
    if (byGroup[div]) return;
    const { rows } = await DataLoader.loadTSV(`data/mtfs_${div}.tsv`);
    const m = {};
    for (const r of rows) (m[r.group] ||= { spec: [], nons: [] })[r.category === "specific" ? "spec" : "nons"]
      .push([r.tf, r.jsd_rank, r.group_expr_log2tpm]);
    for (const g in m) { m[g].spec.sort((a, b) => a[1] - b[1]); m[g].nons.sort((a, b) => a[1] - b[1]); }
    byGroup[div] = m;
  }

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

  function fillPicker() {
    curOpts = optionsFor();
    combo.setOptions(curOpts);
    // land on the most illustrative group, the one with the most specific MTFs (line level: first line)
    let def = curOpts[0];
    if (!isLine() && byGroup[curDiv]) def = curOpts.reduce((best, o) =>
      (byGroup[curDiv][o.key]?.spec.length || 0) > (byGroup[curDiv][best.key]?.spec.length || 0) ? o : best, curOpts[0]);
    combo.setValue(def ? def.key : null);
    if (def) selectVal(def.key);
  }

  const setDesc = html => U.el("atlas-desc").innerHTML = html;
  const row = ([tf, rank, expr]) =>
    `<div class="mtf" title="${tf}: rank #${rank} of 1,651 by CaCTS specificity · mean expression ${expr.toFixed(1)} log₂TPM">
      <div class="g"><span class="tf-sym">${tf}</span></div><div class="meta">rank <b>#${rank}</b> · ${expr.toFixed(1)}</div></div>`;

  function renderLists(spec, nons) {
    curSpec = spec; curNons = nons;            // keep for TSV download
    U.el("atlas-hint").textContent = `${spec.length} specific · ${nons.length} non-specific`;
    U.el("atlas-spec").innerHTML = spec.length ? spec.map(row).join("")
      : `<div class="empty">No factor clears both the top-5% specificity and top-5% expression cutoffs here.</div>`;
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
      const d = byGroup[curDiv][o.key];
      renderLists(d ? d.spec : [], d ? d.nons : []);
      let ds = `${o.n.toLocaleString()} cell line${o.n === 1 ? "" : "s"} in this group`;
      if (curDiv === "modeltype" && mtDesc && mtDesc[o.key]) ds = `<b>${U.esc(o.key)}</b> &rarr; ${U.esc(mtDesc[o.key])} &middot; ${ds}`;
      setDesc(ds);
    }
  }

  function download() {
    const o = curOpts.find(x => x.key === combo.getKey()) || {};
    const g = o.name || o.key || "group";
    const rows = curSpec.map(x => ({ tf: x[0], rank: x[1], expr: x[2], cat: "specific" }))
      .concat(curNons.map(x => ({ tf: x[0], rank: x[1], expr: x[2], cat: "non_specific" })));
    const cols = [{ label: "tf", key: "tf" }, { label: "cacts_rank", key: "rank" },
      { label: "expr_log2tpm", get: r => (+r.expr).toFixed(4) }, { label: "class", key: "cat" }];
    const safe = g.replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "");
    U.downloadTSV(`pycacts_MTFs_${curDiv}_${safe}.tsv`, cols, rows);
  }

  async function pick(div) { curDiv = div; await ensure(div); fillPicker(); }

  async function init() {
    [manifest, mtDesc] = await Promise.all([
      DataLoader.loadJSON("data/manifest.json"), DataLoader.loadJSON("data/modeltype_desc.json")]);
    const seg = U.el("atlas-div");
    seg.innerHTML = U.DIVS.map(([k, lab]) =>
      `<button role="tab" data-div="${k}" aria-selected="${k === curDiv}" title="group cell lines by ${lab.toLowerCase()}">${lab}</button>`).join("");
    seg.addEventListener("click", async e => {
      const b = e.target.closest("button"); if (!b) return;
      [...seg.children].forEach(x => x.setAttribute("aria-selected", x === b));
      await pick(b.dataset.div);
    });
    combo = Combo.make(U.el("atlas-group"), key => selectVal(key));
    U.el("atlas-dl").addEventListener("click", download);
    await ensure(curDiv); fillPicker();
  }
  return { init };
})();
