# Follow-up Questions — Analysis Briefing

**Date:** 2026-09-02  
**Analysis context:** T30 (post-stress) samples; DESeq2 interaction model
`~ Tol + Treatment + Tol:Treatment` (LRT, reduced model `~ Tol + Treatment`),
contrast `TolST:TreatmentNaCl`.  
**DEGs:** padj < 0.05 and |log2FC| ≥ 1 → **2,551 genes**
(2,337 more NaCl-induced in ST; 214 more suppressed in ST).

---

## Q1 — Are PPR40, AS-1, and trmH differentially expressed?

| Gene | Phytozome ID | NCBI Gene ID | qPCR conditions |
|------|-------------|-------------|-----------------|
| PPR40 | AUR62042953 | 110697655 | CONT / NaCl |
| AS-1  | AUR62030471 | 110710750 | CONT only |
| trmH  | AUR62044515 | 110703716 | CONT / NaCl |

**NaCl main effect (T30, model `~ Tol + Treatment`):**

| Gene | log2FC | pvalue | padj | Significant |
|------|--------|--------|------|-------------|
| PPR40 | −0.275 | 0.589 | 0.732 | No |
| AS-1  | +0.110 | 0.833 | 0.905 | No |
| trmH  | +0.175 | 0.818 | 0.895 | No |

Effect sizes are negligible (|log2FC| < 0.3). None of the three genes are
among the significant interaction DEGs either.

**Conclusion:** All three genes are **stably expressed** across NaCl/CONT
conditions and ST/NT classes at T30. This is consistent with their use as
RT-qPCR reference/housekeeping genes rather than stress-responsive genes.

Note: the Phytozome expression profiles indicate these genes are *expressed*
under NaCl conditions — not that they are *differentially expressed*. Our
analysis confirms the former; the latter is not the case at T30.

**Potential follow-up (to align with qPCR data):**
Report ST-internal NaCl effect, NT-internal NaCl effect, and ST vs NT under
CONT separately, with normalised count boxplots (ST/NT × CONT/NaCl).

---

## Q2 — Can genes be grouped by expression profile across conditions?

Yes — we ran `degPatterns` (DEGreport) on the top 1,000 interaction DEGs
(T30 only). 947 genes were assigned to **4 clusters**:

| Cluster | n | Expression pattern (T30) |
|---------|---|--------------------------|
| 1 | 257 | Crossover: ST↑ / NT↓ under NaCl |
| 2 | 495 | ST-specific strong induction; NT largely unchanged |
| 3 | 178 | NT collapses under NaCl; ST stable |
| 4 | 17  | ST-specific downregulation; NT stable |

**Current limitation:** the x-axis covers only Treatment (CONT → NaCl);
T0 samples are not included in this clustering.

**Possible extension:** rebuild the VST matrix using all time points and run
`degPatterns` with `time = "Time"` (or a combined `Treatment:Time` factor on
the x-axis), faceted by ST/NT — this would show trajectories across
CONT-T0 / CONT-T30 / NaCl-T0 / NaCl-T30 for each cluster.

> Note: T0 precedes NaCl application, so by design CONT-T0 ≈ NaCl-T0.
> Including T0 is still informative: it distinguishes stress-induced divergence
> from pre-existing baseline differences between ST and NT.

### GO enrichment per cluster (genome-wide background)

Universe = ~8,218 Arabidopsis orthologs of all expressed quinoa genes
(not just DEGs — the statistically correct choice):

| Cluster | Terms | Key finding |
|---------|-------|-------------|
| 2 | 40 | **Thylakoid** (CC, p.adj 4×10⁻³³); **Photosynthesis** (BP, p.adj 2×10⁻¹⁷); light harvesting, photosynthetic membrane, plastid stroma → ST actively upregulates the chloroplast light-harvesting apparatus under NaCl |
| 3 | 23 | **Pigment metabolic process** (p.adj 3×10⁻⁵); tetrapyrrole/porphyrin metabolism; chlorophyll metabolism; lipid/wax biosynthesis; **cell wall** (CC) → NT loses chlorophyll biosynthesis capacity and cell wall integrity under NaCl |
| 1 | 2  | Phototropism (p.adj 0.019); response to lipid (p.adj 0.043) → divergent light and membrane signalling between ST and NT |
| 4 | — | < 5 orthologs; enrichment not run |

Clusters 2 and 3 together define opposite poles of the same physiological axis:
ST protects and induces photosynthetic capacity while NT simultaneously loses
pigment synthesis and structural integrity.

---

## Q3 — Could poorly clustering replicates affect peer review?

**Automated outlier detection result:** 0 samples flagged.

We ran PCA across all 73 samples and flagged any replicate > 2 SD from its
group centroid (groups: Tol × Treatment × Time). No samples met this threshold.

**Why the visual impression may differ from the statistical result:** the
grouping merges the three accessions within each ST/NT class. Visual clustering
issues observed at the accession level (accessions 14/227/287/460/577/670) can
be masked at this coarser grouping level.

**For peer review:** one mildly outlying replicate per group is common in
RNA-seq and is rarely fatal. With n = 3 per group, one noisy replicate mainly
reduces power (raises dispersion estimates) rather than biasing fold-change
estimates; DESeq2's dispersion shrinkage, Cook's distance, and independent
filtering handle extreme cases.

**Recommended strategy:**
- Include a QC figure (PCA + sample-correlation heatmap) as a supplementary
  figure — pre-empt the reviewer concern rather than waiting for it.
- Quantify whether visually atypical replicates are true outliers or simply
  noisier (within-group distances, silhouette score).
- Run a sensitivity analysis: exclude those replicates, re-run DESeq2, and
  show that key results are stable (high DEG overlap; photosynthesis GO terms
  still enriched). This is the strongest possible response to a reviewer.
- Do not silently remove samples; document the decision with justification.

**Suggested next step:** re-run the QC grouped by accession × Treatment × Time
(sample-correlation heatmap + hierarchical clustering) to pinpoint exactly which
accession's replicate is the concern, then proceed to the sensitivity analysis.

---

## Key output files

| Content | Path |
|---------|------|
| Q1 NaCl main effect | `results/04_cluster/Q1_qPCR_NaCl_main_effect.tsv` |
| Q2 cluster plot | `results/04_cluster/Q2_degPatterns_clusters.png` |
| Q2 cluster assignments | `results/04_cluster/Q2_cluster_assignments.tsv` |
| Q3 outlier PCA | `results/04_cluster/Q3_PCA_outliers.png` |
| Q3 outlier sample list | `results/04_cluster/Q3_outlier_samples.tsv` |
| Per-cluster GO (all) | `results/05_cluster_go/all_clusters_GO_corrected_combined.tsv` |
| Per-cluster GO dotplots | `results/05_cluster_go/cluster{1,2,3}_GO_corrected_dotplot.png` |
| Full analysis report | `results/analysis_report_en.md` |
