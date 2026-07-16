/* compare.js: Tumor vs model. A cancer's specific master TFs in the TCGA tumor vs its DepMap cell-line
   models, split into shared / tumor-only / model-only. Pairing is the curated data/crosswalk.json. */
const Compare = (() => {
  let crosswalk = null, combo, cur = null, curLists = null;

  function specific(mtfs, group) {              // [tf, rank, expr] specific MTFs for a group, jsd-ordered
    const out = [];
    for (const r of mtfs.rows) if (r.group === group && r.category === "specific")
      out.push([r.tf, r.jsd_rank, r.group_expr_log2tpm]);
    out.sort((a, b) => a[1] - b[1]);
    return out;
  }

  const tfRow = ([tf, rank]) =>
    `<div class="mtf" title="${tf}: CaCTS specificity rank #${rank}"><div class="g"><span class="tf-sym">${tf}</span></div><div class="meta">#${rank}</div></div>`;
  const fillCard = (id, list, empty) =>
    U.el(id).innerHTML = list.length ? list.map(tfRow).join("") : `<div class="empty">${empty}</div>`;

  async function selectVal(key) {
    const c = crosswalk.find(x => x.tcga === key); if (!c) return;
    cur = c;
    const [tM, dM] = await Promise.all([
      DataLoader.loadTSV("data/tcga/mtfs_type.tsv"), DataLoader.loadTSV(`data/mtfs_${c.div}.tsv`)]);
    const tumor = specific(tM, c.tcga), model = specific(dM, c.group);
    const tset = new Set(tumor.map(x => x[0])), mset = new Set(model.map(x => x[0]));
    const shared = tumor.filter(x => mset.has(x[0])), tumorOnly = tumor.filter(x => !mset.has(x[0])),
          modelOnly = model.filter(x => !tset.has(x[0]));
    curLists = { shared, tumorOnly, modelOnly };
    U.el("cmp-desc").innerHTML =
      `<b>TCGA ${U.esc(c.tcga)}</b> (${c.n_tcga.toLocaleString()} tumors) &nbsp;vs&nbsp; <b>DepMap ${U.esc(c.group)}</b> (${c.n_dep} cell line${c.n_dep === 1 ? "" : "s"})`;
    U.el("cmp-summary").innerHTML =
      [[shared.length, "shared"], [tumorOnly.length, "tumor only"], [modelOnly.length, "model only"]]
        .map(([n, l]) => `<div class="ps-item"><span class="ps-n">${n}</span><span class="ps-l">${l}</span></div>`).join("");
    fillCard("cmp-shared", shared, "No master TF is specific in both.");
    fillCard("cmp-tumor", tumorOnly, "None.");
    fillCard("cmp-model", modelOnly, "None.");
  }

  function download() {
    if (!cur || !curLists) return;
    const rows = [];
    const add = (list, cat) => list.forEach(([tf, rank]) => rows.push({ tf, rank, cat }));
    add(curLists.shared, "shared"); add(curLists.tumorOnly, "tumor_only"); add(curLists.modelOnly, "model_only");
    const cols = [{ label: "tf", key: "tf" }, { label: "cacts_rank", key: "rank" }, { label: "class", key: "cat" }];
    U.downloadTSV(`pycacts_compare_${cur.tcga}_vs_${cur.group.replace(/[^A-Za-z0-9]+/g, "_")}.tsv`, cols, rows);
  }

  async function init() {
    crosswalk = await DataLoader.loadJSON("data/crosswalk.json");
    combo = Combo.make(U.el("cmp-group"), key => selectVal(key));
    combo.setOptions(crosswalk.map(c => ({
      key: c.tcga, label: c.tcga_label, search: `${c.tcga} ${c.tcga_label} ${c.group}`.toLowerCase() })));
    const def = crosswalk.find(c => c.tcga === "SKCM") || crosswalk[0];
    combo.setValue(def.tcga);
    U.el("cmp-dl").addEventListener("click", download);
    selectVal(def.tcga);
  }
  return { init };
})();
