/* about.js, what CaCTS/pyCaCTS computes, how an MTF is called, the DepMap levels, validation, credit. */
const About = (() => {
  function init() {
    U.el("about-body").innerHTML = `
      <h3>What pyCaCTS computes</h3>
      <p>For each transcription factor, pyCaCTS scores how <b>specifically</b> it is expressed in one group
      relative to the diversity of all groups, the <b>Jensen-Shannon divergence</b> between the factor's
      normalized cross-group expression profile and a one-hot "ideal" at the query group. A <b>lower</b>
      score means more group-specific (rank 1 = most specific). It uses expression only, no ChIP-seq.</p>

      <h3>The DepMap organizational levels</h3>
      <p>pyCaCTS is run here on the DepMap/CCLE panel (~1,450 cell lines) grouped at every level of the
      Oncotree disease hierarchy, so specificity can be read at whatever resolution is useful:</p>
      <ul>
        <li><b>Lineage</b> (<code>OncotreeLineage</code>): 29 groups, e.g. <em>Skin</em>.</li>
        <li><b>Primary disease</b> (<code>OncotreePrimaryDisease</code>): 79 groups, e.g. <em>Melanoma</em>.</li>
        <li><b>Subtype</b> (<code>OncotreeSubtype</code>): 191 groups, e.g. <em>Cutaneous Melanoma</em>.</li>
        <li><b>Model type</b> (<code>DepmapModelType</code>): 192 groups, the finest Oncotree code (e.g. <em>GB → Glioblastoma</em>).</li>
        <li><b>Cell line</b>, each of the ~1,450 individual DepMap models scored on its own (its expression profile vs. the whole panel).</li>
      </ul>
      <p>The same 1,651 expressed TFs (of the CaCTS 1,671-TF catalogue) are scored at every level; a group's
      "representative profile" is the mean expression across its member cell lines. Coarser levels pool more
      lines (more robust, less granular); finer levels resolve individual subtypes.</p>

      <h3>TCGA tumors</h3>
      <p>The <b>Data</b> toggle (on the MTF-atlas and TF-scores tabs) switches from the DepMap cell-line
      panel to <b>TCGA</b> tumors, scored the same way at three levels: <b>tumor type</b> (33, e.g. SKCM,
      BRCA, LUAD), <b>molecular subtype</b> (~94 curated TCGA subtype calls, e.g. BRCA.LumA, GBM.Proneural),
      and <b>sample type</b> (primary / metastatic / recurrent / normal, from the sample barcode).
      Expression is the UCSC Xena <b>Toil</b> RSEM <b>TPM</b> matrix (TCGA uniformly re-quantified), re-encoded
      to log₂(TPM+1) so it is in the same units as DepMap and the <b>1 TPM</b> abundance floor applies
      identically to both panels; tumor types follow the CaCTS 33-type map and subtypes the TCGA
      subtype working-group calls. Specificity is computed relative to each
      dataset's own reference set, so a TF's DepMap and TCGA scores are independent, and the difference is
      the point (the tumor vs its cell-line models). Two caveats: there is no CRISPR essentiality for tumors
      (that column is blank for TCGA), and the sample-type level is cross-cutting across all cancers, so it
      reads most cleanly within a single tumor type. The <b>Tumor vs model</b> tab pairs each tumor type
      with its matched DepMap group (a curated crosswalk) and splits their specific master TFs into
      <b>shared</b> (preserved in the model), <b>tumor-only</b>, and <b>model-only</b> (e.g. ovarian: PAX8
      shared, but WT1 / SOX17 tumor-only, lost in the cell-line models).</p>

      <h3>Within-cancer master TFs</h3>
      <p>The <b>Within-cancer</b> tab re-scores CaCTS with the reference set restricted to a <b>single
      cancer's own samples</b>, so specificity is read <em>inside</em> one disease rather than across all
      of TCGA. Two axes are offered: <b>molecular subtype</b> (the TCGA subtype working-group calls, e.g.
      BRCA into LumA / LumB / Basal / Her2 / Normal) and <b>tumor vs adjacent-normal</b> tissue (from the
      sample barcode). A <b>Call</b> toggle picks the definition, because within one shared lineage the two
      diverge. <b>Significant</b> (the tool-wide call: FDR &lt; 0.10 AND mean ≥ 1 TPM) favours the most
      subgroup-<em>discriminative</em> factors: breast LumA surfaces the progesterone receptor <b>PGR</b>, and
      some subgroups have many such TFs while some (e.g. lung adenocarcinoma normal tissue) have none, which is
      itself informative. <b>Abundant</b> (the top-5% expressed TFs ranked by specificity) instead surfaces the
      highly-expressed <em>canonical</em> masters that a low abundance floor lets ultra-specific but modestly
      expressed TFs outrank: breast LumA / LumB → <b>ESR1, FOXA1, GATA3</b>, Basal → FOXM1 / MYBL2. Metastatic
      and normal samples inherit their patient's cancer via the sample barcode. A stage axis is not shown: the
      AJCC-stage clinical table is not freely fetchable from the Xena mirror.</p>

      <h3>How an MTF is called</h3>
      <p>A factor is a <b>specific MTF</b> in a group if it is both <b>significantly group-specific</b>
      (empirical-null <b>FDR &lt; 0.10</b> on its CaCTS score, see below) <em>and</em> <b>expressed</b>
      (mean <b>≥ 1 TPM</b>, i.e. log₂(TPM+1) ≥ 1). A <b>non-specific MTF</b>, CaCTS's <em>candidate ubiquitous
      (multi-cancer) master regulator</em>, is highly expressed (top 5%) but not group-specific (FDR ≥ 0.10).
      This replaces CaCTS's original fixed cutoffs (top 5% by score ∩ top 5% by expression) with a data-driven
      significance threshold and a light 1 TPM abundance floor: the floor keeps genuinely expressed lineage
      factors that the aggressive top-5%-expression gate dropped (ovarian <b>SOX17 / WT1 / MECOM</b> are all
      recovered as specific) while still excluding near-silent JSD artifacts.</p>
      <div class="callout"><b>Two orthogonal axes.</b> Specificity (CaCTS score / FDR) and abundance
      (expression) are independent. The MTF call requires both; the <b>TF scores</b> tab shows the FDR and
      expression as their own sortable columns, so you can see factors that are highly specific yet sit below
      the abundance floor, which CaCTS's authors note "may also be MTFs."</div>

      <h3>The empirical-null FDR</h3>
      <p>The specificity gate is an <b>empirical-null FDR</b>, computed per TF for the selected group. The
      group's non-specific majority (the high-JSD side of its own score distribution) defines the null; each TF
      gets a left-tail p-value against it; Benjamini-Hochberg gives the FDR. It answers <em>"is this TF
      significantly more group-specific than the background?"</em>, a data-driven, non-arbitrary threshold that
      replaces CaCTS's fixed top-5%-by-score cutoff. The <b>TF scores</b> tab shows it per TF and lets you
      filter further (e.g. FDR &lt; 0.05); it is recomputed for whichever group and level you select.</p>

      <h3>The TF-scores table</h3>
      <p>The <b>TF scores</b> tab puts every TF in one table for the selected group. Alongside the CaCTS
      score and class it shows the <b>two gates</b> behind the class as their own columns (empirical-null
      FDR &lt; 0.10, and mean expression ≥ 1 TPM), each TF's <b>family</b> (DNA-binding domain, Lambert et al. 2018), <b>per-group
      CRISPR essentiality</b> (mean Chronos across the group's DepMap cell lines; lower = stronger
      dependency, ~ -1 = common-essential, 0 = neutral; not staged at the single-cell-line level),
      <b>cross-group breadth</b> (in how many lineage / disease / subtype / model-type groups the TF is a
      specific MTF), and out-links to NCBI Gene, GeneCards, and DepMap. The <b>&#8595; TSV</b> button exports
      the current sorted and filtered view; the MTF-atlas lists export the same way.</p>

      <h3>Validation &amp; performance</h3>
      <p>pyCaCTS reproduces the original CaCTS R output to floating-point precision
      (max |Δ| &lt; 1×10⁻¹⁵ on identical input) and runs roughly <b>~3,900× faster</b> via an exact
      vectorization of the per-query loop: scoring the whole DepMap panel at all five levels takes pyCaCTS
      <b>0.19 s</b> versus <b>12.5 min</b> for the original R. Input: DepMap/CCLE
      <code>OmicsExpressionProteinCodingGenesTPMLogp1</code>.</p>

      <h3>Credit</h3>
      <p>CaCTS is the method of the Lawrenson lab, <b>Reddy J, Fonseca MAS, Corona RI, et al.,
      &ldquo;Predicting master transcription factors from pan-cancer expression data,&rdquo; <em>Sci. Adv.</em>
      2021;7(48):eabf6123</b> (PMID 34818047; DOI 10.1126/sciadv.abf6123),
      <a href="https://github.com/lawrenson-lab/CaCTS" target="_blank" rel="noopener">github.com/lawrenson-lab/CaCTS</a>.
      pyCaCTS is an independent Python reimplementation; all credit for the method is theirs. Their original
      R is GPL and is not redistributed here.</p>`;
  }
  return { init };
})();
