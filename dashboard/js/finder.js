/* finder.js, search a TF, list every group where it is a specific MTF. */
const Finder = (() => {
  let index = null;
  const SHORT = Object.fromEntries(U.DIVS.map(([k, lab]) => [k, lab.toLowerCase()]));
  async function ensure() { index ||= await DataLoader.loadJSON("data/tf_index.json"); }
  function search(q) {
    const box = U.el("finder-res");
    q = q.trim().toUpperCase();
    if (!q) { box.innerHTML = `<div class="hint" style="padding:6px 2px">Start typing a gene symbol.</div>`; return; }
    const hits = Object.keys(index).filter(t => t.startsWith(q)).sort().slice(0, 14);
    if (!hits.length) {
      box.innerHTML = `<div class="card fr"><span class="hint">No group has <b class="mono">${U.esc(q)}</b> as a <em>specific</em> MTF. It may still rank highly by specificity without clearing the expression cutoff, check the TF scores tab.</span></div>`;
      return;
    }
    box.innerHTML = hits.map(tf => {
      const rows = index[tf].slice().sort((a, b) => a[2] - b[2]);
      const chips = rows.map(r => `<span class="chip"><span class="r">${SHORT[r[0]] || r[0]}</span> <b>${U.esc(r[1])}</b> <span class="r">#${r[2]}</span></span>`).join("");
      return `<div class="card fr"><span class="sym">${tf}</span>${U.tags(tf)}<div class="chips">${chips}</div></div>`;
    }).join("");
  }
  async function init() {
    await ensure();
    const inp = U.el("finder-input");
    inp.addEventListener("input", e => search(e.target.value));
    search("");
  }
  return { init };
})();
