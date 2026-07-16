/* within.js: Within-cancer master TFs. For one cancer, CaCTS is scored across ITS OWN groups (reference
   set = that cancer's samples only), on two axes: molecular subtype, or tumor vs adjacent-normal tissue.
   Each group's most-specific expressed TFs are shown side by side, one card per group. */
const Within = (() => {
  let manifest = null, typeDesc = {}, combo, axis = "subtype", cur = null, curRows = null;
  const AXES = [["subtype", "Molecular subtype"], ["tumornormal", "Tumor vs normal"]];
  const file = a => `data/tcga/within_${a}_mtfs.tsv`;
  const DEFAULT = "BRCA";                                // best showcase (luminal ESR1 / FOXA1 / GATA3)

  const cancerLabel = c => typeDesc[c] ? `${c}: ${typeDesc[c]}` : c;
  function shortGroup(g) {                               // "BRCA.LumA" -> "LumA", "GBM_LGG.Classic-like" -> "Classic-like"
    const rest = g.includes(".") ? g.slice(g.indexOf(".") + 1) : g;
    return /^\d+$/.test(rest) ? `Subtype ${rest}` : (rest || g);
  }

  const tfRow = r => {
    const spec = r.category === "specific";
    const badge = spec ? ` <span class="pill spec tiny" title="specific: top-5% CaCTS score AND top-5% expression in this group">S</span>` : "";
    return `<div class="mtf" title="${U.esc(r.tf)}: rank #${r.rank} by CaCTS specificity here; ${r.group_expr_log2tpm} log₂TPM expression">
      <div class="g"><span class="tf-rank">${r.rank}</span><span class="tf-sym">${U.esc(r.tf)}</span>${badge}</div>
      <div class="meta">${Number(r.group_expr_log2tpm).toFixed(1)}</div></div>`;
  };

  async function selectVal(code) {
    cur = code;
    const mtfs = await DataLoader.loadTSV(file(axis));
    const rows = mtfs.rows.filter(r => r.cancer === code);
    curRows = rows;
    const groups = {};
    for (const r of rows) (groups[r.group] ||= []).push(r);
    // biggest group first; tumor/normal axis forced Tumor-then-Normal
    const order = axis === "tumornormal" ? ["Tumor", "Normal"] : null;
    const gnames = Object.keys(groups).sort((a, b) =>
      order ? order.indexOf(a) - order.indexOf(b) : (groups[b].length - groups[a].length) || a.localeCompare(b));

    const axisWord = axis === "subtype" ? "molecular subtypes" : "sample states (tumor vs adjacent normal)";
    U.el("within-desc").innerHTML =
      `<b>${U.esc(cancerLabel(code))}</b> scored across its own ${gnames.length} ${axisWord} (reference set = this cancer only). ` +
      `Each column lists that group's most-specific expressed master TFs (CaCTS rank within the group); <span class="pill spec tiny">S</span> marks a strict specific MTF.`;

    U.el("within-cards").innerHTML = gnames.map(g => {
      const gr = groups[g].slice().sort((a, b) => a.rank - b.rank);
      const n = gr[0] ? gr[0].group_size : 0;
      return `<div class="card"><div class="card-h"><h3 title="${U.esc(g)}">${U.esc(shortGroup(g))}</h3>` +
        `<span class="muted-s" title="samples in this group">${n}</span></div>` +
        `<div class="card-b listbox">${gr.map(tfRow).join("")}</div></div>`;
    }).join("");
  }

  function setAxis(a) {
    axis = a;
    U.el("within-axis").querySelectorAll("button").forEach(b =>
      b.setAttribute("aria-selected", String(b.dataset.axis === a)));
    const cancers = Object.keys(manifest[axis]).sort();
    combo.setOptions(cancers.map(c => ({ key: c, label: cancerLabel(c), search: `${c} ${typeDesc[c] || ""}`.toLowerCase() })));
    const pick = cancers.includes(cur) ? cur : (cancers.includes(DEFAULT) ? DEFAULT : cancers[0]);
    combo.setValue(pick);
    selectVal(pick);
  }

  function download() {
    if (!cur || !curRows) return;
    const cols = ["group", "group_size", "tf", "rank", "category", "cacts_score", "group_expr_log2tpm"]
      .map(k => ({ label: k, key: k }));
    U.downloadTSV(`pycacts_within_${axis}_${cur}.tsv`, cols, curRows);
  }

  async function init() {
    [manifest, typeDesc] = await Promise.all([
      DataLoader.loadJSON("data/tcga/within_manifest.json"),
      DataLoader.loadJSON("data/tcga/type_desc.json").catch(() => ({})),
    ]);
    U.el("within-axis").innerHTML = AXES.map(([k, lab]) =>
      `<button data-axis="${k}" aria-selected="${k === axis}" title="score each cancer across its ${k === "subtype" ? "molecular subtypes" : "tumor vs adjacent-normal samples"}">${lab}</button>`).join("");
    U.el("within-axis").querySelectorAll("button").forEach(b =>
      b.addEventListener("click", () => setAxis(b.dataset.axis)));
    combo = Combo.make(U.el("within-group"), key => selectVal(key));
    U.el("within-dl").addEventListener("click", download);
    setAxis(axis);
  }
  return { init };
})();
