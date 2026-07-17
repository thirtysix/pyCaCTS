/* within.js: Within-cancer master TFs. For one cancer, CaCTS is scored across ITS OWN groups (reference
   set = that cancer's samples only), on two axes (molecular subtype, or tumor vs adjacent-normal) and with
   a Call toggle: "Significant" = empirical-null FDR < 0.10 AND mean >= 1 TPM (the tool-wide call);
   "Abundant" = top-5% expressed, ranked by specificity (the canonical highly-expressed masters, e.g.
   breast LumA -> ESR1 / FOXA1 / GATA3). Within one lineage the two views diverge; both are shown. */
const Within = (() => {
  let manifest = null, typeDesc = {}, combo, axis = "subtype", call = "fdr", cur = null, curRows = null;
  const AXES = [["subtype", "Molecular subtype"], ["tumornormal", "Tumor vs normal"]];
  const CALLS = [["fdr", "Significant"], ["abundant", "Abundant"]];
  const file = a => `data/tcga/within_${a}_mtfs.tsv`;
  const DEFAULT = "BRCA";
  const CAP = 40;                                        // TFs shown per card (rest in the TSV)

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

  const tfRow = (r, i, abundant) => {
    const dot = abundant && r.category === "specific"
      ? ` <span class="pill spec tiny" title="also empirical-null FDR &lt; 0.10">S</span>` : "";
    const meta = abundant ? (+r.group_expr_log2tpm).toFixed(1) : fmtFDR(+r.fdr_log10);
    return `<div class="mtf" title="${U.esc(r.tf)}: CaCTS rank #${r.jsd_rank} within this cancer; empirical-null FDR ${fmtFDR(+r.fdr_log10)}; ${r.group_expr_log2tpm} log₂(TPM+1)">
      <div class="g"><span class="tf-rank">${i}</span><span class="tf-sym">${U.esc(r.tf)}</span>${dot}</div>
      <div class="meta" title="${abundant ? "mean log₂(TPM+1)" : "empirical-null FDR"}">${meta}</div></div>`;
  };

  async function selectVal(code) {
    cur = code;
    const abundant = call === "abundant";
    const mtfs = await DataLoader.loadTSV(file(axis));
    const rows = mtfs.rows.filter(r => r.cancer === code);
    curRows = rows;
    const keep = r => abundant ? (+r.top5e === 1) : (r.category === "specific");
    const groups = {};
    for (const r of rows) if (keep(r)) (groups[r.group] ||= []).push(r);
    const allGroups = [...new Set(rows.map(r => r.group))];
    const order = axis === "tumornormal" ? ["Tumor", "Normal"] : null;
    const sizeOf = g => { const any = rows.find(r => r.group === g); return any ? +any.group_size : 0; };
    const gnames = allGroups.sort((a, b) =>
      order ? order.indexOf(a) - order.indexOf(b) : sizeOf(b) - sizeOf(a) || a.localeCompare(b));

    const axisWord = axis === "subtype" ? "molecular subtypes" : "sample states (tumor vs adjacent normal)";
    U.el("within-desc").innerHTML =
      `<b>${U.esc(cancerLabel(code))}</b> scored across its own ${gnames.length} ${axisWord} (reference set = this cancer only). ` +
      (abundant
        ? `Each column lists that subgroup's <b>abundant</b> master TFs (top-5% expressed here), ordered by CaCTS specificity; <span class="pill spec tiny">S</span> marks those also FDR-significant. The value shown is expression (log₂TPM+1). Within one lineage this surfaces the highly-expressed canonical masters.`
        : `Each column lists that subgroup's <b>specific</b> master TFs (empirical-null FDR &lt; 0.10 AND mean ≥ 1 TPM), ordered by CaCTS specificity; the value shown is the FDR.`);

    const noun = abundant ? "abundant" : "specific";
    U.el("within-cards").innerHTML = gnames.map(g => {
      const gr = (groups[g] || []).slice().sort((a, b) => (+a.cacts_score) - (+b.cacts_score));
      const n = sizeOf(g), more = gr.length > CAP ? `<div class="empty">…and ${gr.length - CAP} more (in the TSV)</div>` : "";
      const body = gr.length
        ? gr.slice(0, CAP).map((r, i) => tfRow(r, i + 1, abundant)).join("") + more
        : `<div class="empty">None.</div>`;
      return `<div class="card"><div class="card-h"><h3 title="${U.esc(g)}">${U.esc(shortGroup(g))}</h3>` +
        `<span class="muted-s" title="${gr.length} ${noun} TFs · ${n} samples">${gr.length} · n=${n}</span></div>` +
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

  function setCall(c) {
    call = c;
    U.el("within-call").querySelectorAll("button").forEach(b =>
      b.setAttribute("aria-selected", String(b.dataset.call === c)));
    if (cur) selectVal(cur);
  }

  function download() {
    if (!cur || !curRows) return;
    const cols = ["group", "group_size", "tf", "category", "jsd_rank", "cacts_score", "group_expr_log2tpm",
                  "fdr_log10", "top5e"].map(k => ({ label: k, key: k }));
    U.downloadTSV(`pycacts_within_${axis}_${cur}.tsv`, cols, curRows);
  }

  const seg = (id, items, sel, attr, onclick) => {
    U.el(id).innerHTML = items.map(([k, lab, title]) =>
      `<button data-${attr}="${k}" aria-selected="${k === sel}"${title ? ` title="${title}"` : ""}>${lab}</button>`).join("");
    U.el(id).querySelectorAll("button").forEach(b => b.addEventListener("click", () => onclick(b.dataset[attr])));
  };

  async function init() {
    [manifest, typeDesc] = await Promise.all([
      DataLoader.loadJSON("data/tcga/within_manifest.json"),
      DataLoader.loadJSON("data/tcga/type_desc.json").catch(() => ({})),
    ]);
    seg("within-axis", AXES.map(([k, l]) => [k, l, `score each cancer across its ${k === "subtype" ? "molecular subtypes" : "tumor vs adjacent-normal samples"}`]), axis, "axis", setAxis);
    seg("within-call", [
      ["fdr", "Significant", "empirical-null FDR &lt; 0.10 AND mean ≥ 1 TPM (the tool-wide call)"],
      ["abundant", "Abundant", "top-5% expressed, ranked by specificity: the highly-expressed canonical masters (e.g. LumA → ESR1 / FOXA1 / GATA3)"],
    ], call, "call", setCall);
    combo = Combo.make(U.el("within-group"), key => selectVal(key));
    U.el("within-dl").addEventListener("click", download);
    setAxis(axis);
  }
  return { init };
})();
