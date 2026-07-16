/* within.js: Within-cancer master TFs. For one cancer, CaCTS is scored across ITS OWN groups (reference
   set = that cancer's samples only), on two axes: molecular subtype, or tumor vs adjacent-normal tissue.
   Each group's specific MTFs (empirical-null FDR < 0.10 AND mean >= 1 TPM) are shown side by side. */
const Within = (() => {
  let manifest = null, typeDesc = {}, combo, axis = "subtype", cur = null, curRows = null;
  const AXES = [["subtype", "Molecular subtype"], ["tumornormal", "Tumor vs normal"]];
  const file = a => `data/tcga/within_${a}_mtfs.tsv`;
  const DEFAULT = "BRCA";
  const CAP = 40;                                        // most-specific TFs shown per card (rest in the TSV)

  const cancerLabel = c => typeDesc[c] ? `${c}: ${typeDesc[c]}` : c;
  function shortGroup(g) {                               // "BRCA.LumA" -> "LumA", "GBM_LGG.Classic-like" -> "Classic-like"
    const rest = g.includes(".") ? g.slice(g.indexOf(".") + 1) : g;
    return /^\d+$/.test(rest) ? `Subtype ${rest}` : (rest || g);
  }
  const fmtFDR = l10 => {                                // l10 = log10(FDR)
    if (l10 == null || isNaN(l10)) return "–";
    if (l10 >= -3) return (10 ** l10).toFixed(3);
    const e = Math.floor(l10), m = 10 ** (l10 - e);
    return `${m.toFixed(1)}e${e}`;
  };

  const tfRow = (r, i) =>
    `<div class="mtf" title="${U.esc(r.tf)}: CaCTS rank #${r.jsd_rank} within this cancer; empirical-null FDR ${fmtFDR(+r.fdr_log10)}; ${r.group_expr_log2tpm} log₂(TPM+1)">
      <div class="g"><span class="tf-rank">${i}</span><span class="tf-sym">${U.esc(r.tf)}</span></div>
      <div class="meta" title="empirical-null FDR">${fmtFDR(+r.fdr_log10)}</div></div>`;

  async function selectVal(code) {
    cur = code;
    const mtfs = await DataLoader.loadTSV(file(axis));
    const rows = mtfs.rows.filter(r => r.cancer === code);
    curRows = rows;
    const groups = {};
    for (const r of rows) if (r.category === "specific") (groups[r.group] ||= []).push(r);
    // every group that exists (even with 0 specific), biggest first; tumor/normal forced Tumor-then-Normal
    const allGroups = [...new Set(rows.map(r => r.group))];
    const order = axis === "tumornormal" ? ["Tumor", "Normal"] : null;
    const sizeOf = g => { const any = rows.find(r => r.group === g); return any ? +any.group_size : 0; };
    const gnames = allGroups.sort((a, b) =>
      order ? order.indexOf(a) - order.indexOf(b) : sizeOf(b) - sizeOf(a) || a.localeCompare(b));

    const axisWord = axis === "subtype" ? "molecular subtypes" : "sample states (tumor vs adjacent normal)";
    U.el("within-desc").innerHTML =
      `<b>${U.esc(cancerLabel(code))}</b> scored across its own ${gnames.length} ${axisWord} (reference set = this cancer only). ` +
      `Each column lists that subgroup's <b>specific</b> master TFs (empirical-null FDR &lt; 0.10 AND mean ≥ 1 TPM), ordered by CaCTS specificity; the value shown is the FDR.`;

    U.el("within-cards").innerHTML = gnames.map(g => {
      const gr = (groups[g] || []).slice().sort((a, b) => (+a.cacts_score) - (+b.cacts_score));
      const n = sizeOf(g), more = gr.length > CAP ? `<div class="empty">…and ${gr.length - CAP} more (in the TSV)</div>` : "";
      const body = gr.length
        ? gr.slice(0, CAP).map((r, i) => tfRow(r, i + 1)).join("") + more
        : `<div class="empty">No FDR-significant specific TFs.</div>`;
      return `<div class="card"><div class="card-h"><h3 title="${U.esc(g)}">${U.esc(shortGroup(g))}</h3>` +
        `<span class="muted-s" title="${gr.length} specific TFs · ${n} samples">${gr.length} · n=${n}</span></div>` +
        `<div class="card-b listbox">${body}</div></div>`;
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
    const cols = ["group", "group_size", "tf", "category", "jsd_rank", "cacts_score", "group_expr_log2tpm", "fdr_log10"]
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
