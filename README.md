# Step 1 — true-variant vs sequencing/mapping-error confidence model

Distinguishes a real variant from a **sequencing error** or **mapping error**
using UMI/co-barcode **molecule-linkage** features, validated against HG002 GIAB
truth. This is the compute-light (CPU, LightGBM) core of the cLFR/cWGS
molecule-linkage variant re-scorer. See `../memory/plan.md` for the full rationale.

## Why molecule features (the differentiator)

A standard per-read pileup caller sees read counts. It does **not** see whether
alt support is consistent *within independent molecules*. A true variant appears
across many independent molecules and consistently within each; a sequencing
error is scattered across single reads; a mapping error rides on low-MAPQ /
soft-clipped / repeat-adjacent reads. `02_extract_features.py` extracts exactly
these molecule-resolved signals (independent molecule count, within-molecule
agreement, MAPQ/BQ/strand/softclip/read-position).

## Truth (Step 1 vs Step 2)

- **Step 1 (this code):** HG002 gives two gold labels from one pure sample —
  hom-ref confident positions with alt reads = **error** (negatives); GIAB SNV
  sites = **true** (positives, at germline ~50/100% AF). Enough for error
  suppression + calibrated confidence. **No titration needed.**
- **Step 2 (later):** true variants at *low* AF need a GIAB **titration**
  (HG002 + HG003/HG004 at known ratios). Pure HG002 has no low-AF positives.
  Do not claim low-AF sensitivity from Step 1 data.

## Pipeline

```bash
# edit paths in run.sh, then:
bash run.sh
```

1. `01_make_candidates.py` — pileup within confident BED → labeled sites
   (true = GIAB SNV; error = confident hom-ref with alt reads).
2. `02_extract_features.py` — molecule-grouped features per site.
3. `03_train_eval.py` — LightGBM true-vs-error, **chrom-held-out test**,
   isotonic calibration fit on a TRAIN slice, PR/ROC/Brier + VAF-stratified PR-AUC
   + feature importance.

## The money chart (ablation)

`run.sh` trains twice: `--feature-set all` vs `--feature-set baseline`
(baseline = dp/alt_reads/ref_reads/vaf, i.e. what a plain pileup caller has).
The **PR-AUC delta** — largest among higher-VAF errors (mapping artifacts that
VAF alone can't flag) — is the value the molecule/mapping features add.

## Step 4 — from re-score to an isoform with correct SNPs

`04_apply_rescore.py` turns the scorer into a deliverable: run `02` on the
consensus-supporting reads (candidates = consensus-vs-reference diffs), then per
SNP decide `KEEP` (p_true ≥ thr → real SNP), `REVERT` (< thr → error, drop), or
`EDIT` (A>G/T>C at a REDIportal site → RNA editing, kept & annotated, never
reverted). The `KEEP` set → `bcftools consensus` → corrected isoform FASTA.
Validate on ERCC: per-base false-SNP rate before vs after correction. See
`readme_cn.md` for the full workflow. Note: the Step 1 model is DNA-trained;
recalibrate on ERCC for the RNA/isoform application.

## Honesty guardrails (baked in)

- Split by chromosome **before** anything else; calibration fit only on TRAIN.
- No scaling/SMOTE before the split (the leak that inflated the old `sqanti3.py`).
- On pure HG002, true=high-VAF / error=low-VAF, so VAF alone separates much of it;
  report the ablation and VAF strata honestly rather than one headline AUC.
- Step 1 = error suppression + calibration only. Low-AF sensitivity needs Step 2.

## Requires

`pysam, pandas, numpy, scikit-learn, lightgbm` (see `environment.yml`).
Truth VCF must be bgzip+tabix; reference needs `.fai`.