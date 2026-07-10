/* about.js — what CaCTS/pyCaCTS computes, how an MTF is called, the DepMap levels, validation, credit. */
const About = (() => {
  function init() {
    U.el("about-body").innerHTML = `
      <h3>What pyCaCTS computes</h3>
      <p>For each transcription factor, pyCaCTS scores how <b>specifically</b> it is expressed in one group
      relative to the diversity of all groups — the <b>Jensen-Shannon divergence</b> between the factor's
      normalized cross-group expression profile and a one-hot "ideal" at the query group. A <b>lower</b>
      score means more group-specific (rank 1 = most specific). It uses expression only — no ChIP-seq.</p>

      <h3>The DepMap organizational levels</h3>
      <p>pyCaCTS is run here on the DepMap/CCLE panel (~1,450 cell lines) grouped at every level of the
      Oncotree disease hierarchy, so specificity can be read at whatever resolution is useful:</p>
      <ul>
        <li><b>Lineage</b> (<code>OncotreeLineage</code>) — 29 groups, e.g. <em>Skin</em>.</li>
        <li><b>Primary disease</b> (<code>OncotreePrimaryDisease</code>) — 79 groups, e.g. <em>Melanoma</em>.</li>
        <li><b>Subtype</b> (<code>OncotreeSubtype</code>) — 191 groups, e.g. <em>Cutaneous Melanoma</em>.</li>
        <li><b>Model type</b> (<code>DepmapModelType</code>) — 192 groups, the finest Oncotree code (e.g. <em>GB → Glioblastoma</em>).</li>
        <li><b>Cell line</b> — each of the ~1,450 individual DepMap models scored on its own (its expression profile vs. the whole panel).</li>
      </ul>
      <p>The same 1,651 expressed TFs (of the CaCTS 1,671-TF catalogue) are scored at every level; a group's
      "representative profile" is the mean expression across its member cell lines. Coarser levels pool more
      lines (more robust, less granular); finer levels resolve individual subtypes.</p>

      <h3>How an MTF is called</h3>
      <p>Following CaCTS, a factor is a <b>specific MTF</b> in a group if it is in the <b>top 5% by CaCTS score</b>
      <em>and</em> the <b>top 5% by mean expression</b> there. A <b>non-specific MTF</b> — CaCTS's
      <em>candidate ubiquitous (multi-cancer) master regulator</em> — has high expression (top 5%) but low
      lineage-specificity (CaCTS rank outside the top 5%). In a 1,651-TF universe, "top 5%" = the top ~83 TFs.</p>
      <div class="callout"><b>Two orthogonal axes.</b> Specificity (CaCTS score) and abundance (expression) are
      independent. The MTF call requires both; the <b>TF scores</b> tab exposes the pure specificity signal
      (score / rank) on its own, so you can see factors that are highly specific yet don't clear the abundance
      gate — which CaCTS's authors note "may also be MTFs."</div>

      <h3>The empirical-null FDR column</h3>
      <p>The <b>TF scores</b> tab also reports an <b>empirical-null FDR</b> per TF for the selected group, and lets
      you filter by it. The group's non-specific majority (the high-JSD side of its own score distribution) defines
      the null; each TF gets a left-tail p-value against it; Benjamini-Hochberg gives the FDR. It answers
      <em>"is this TF significantly more group-specific than the background?"</em> — a data-driven, non-arbitrary
      alternative to the top-5% cutoff (e.g. filter to FDR &lt; 0.10). It is recomputed for whichever group and
      level you select.</p>

      <h3>Validation &amp; performance</h3>
      <p>pyCaCTS reproduces the original CaCTS R output to floating-point precision
      (max |Δ| &lt; 1×10⁻¹⁵ on identical input) and runs roughly <b>3,000–4,300× faster</b> via an exact
      vectorization of the per-query loop — a genome-wide re-score across all groups is a milliseconds-to-seconds
      job rather than minutes. Input: DepMap/CCLE <code>OmicsExpressionProteinCodingGenesTPMLogp1</code>.</p>

      <h3>Credit</h3>
      <p>CaCTS is the method of the Lawrenson lab — <b>Reddy J, Fonseca MAS, Corona RI, et al.,
      &ldquo;Predicting master transcription factors from pan-cancer expression data,&rdquo; <em>Sci. Adv.</em>
      2021;7(48):eabf6123</b> (PMID 34818047; DOI 10.1126/sciadv.abf6123),
      <a href="https://github.com/lawrenson-lab/CaCTS" target="_blank" rel="noopener">github.com/lawrenson-lab/CaCTS</a>.
      pyCaCTS is an independent Python reimplementation; all credit for the method is theirs. Their original
      R is GPL and is not redistributed here.</p>`;
  }
  return { init };
})();
