# MR_UNAGI

Core analysis scripts for a multiple myeloma study: two-sample cis-eQTL Mendelian randomization (MR) for JAK1/JAK2, and four-stage single-cell model training with UNAGI and iDREM. The scripts start from prepared inputs. This is not the complete preprocessing, perturbation, or figure-generation workflow.

## JAK1/JAK2 MR

Requires R, `data.table`, [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), PLINK 1.9, and a European LD reference in PLINK binary format. The checked R environment used R 4.5.3, data.table 1.18.2.1, and TwoSampleMR 0.7.4.

```sh
Rscript MR.R PROJECT_DIR PLINK_EXECUTABLE EUR_REFERENCE_PREFIX
```

`EUR_REFERENCE_PREFIX` is the shared path prefix of the `.bed`, `.bim`, and `.fam` files. Under `PROJECT_DIR`, supply these tab-separated, gzipped files:

| File | Columns used |
| --- | --- |
| `clean/exposure/JAK1_JAK2_eqtlgen_p5e8_candidate_rows.tsv.gz` | SNP, GeneSymbol, Pvalue, AssessedAllele, OtherAllele, Zscore, NrSamples |
| `clean/exposure/JAK1_JAK2_eqtlgen_af_for_clumped_snps.tsv.gz` | SNP, AlleleA, AlleleB, AlleleB_all (frequency of AlleleB) |
| `clean/outcome/finngen_R12_C3_MYELOMA_MONOCLONAL_for_JAK1_JAK2.tsv.gz` | SNP, beta (log odds), se, effect_allele, other_allele, pval, eaf, id, phenotype |

The [eQTLGen Phase I blood cis-eQTL data](https://www.eqtlgen.org/cis-eqtls.html) supply the exposures. The outcome is [FinnGen](https://www.finngen.fi/en/access_results) release 12 endpoint `C3_MYELOMA_MONOCLONAL`: MM/MGUS combined, not progression from MGUS to MM. The script prints gene-specific odds ratios, 95% confidence intervals, P values, and Benjamini-Hochberg q values across the two genes for this endpoint.

## Four-stage UNAGI training

Requires [UNAGI](https://github.com/mcgilldinglab/UNAGI), Java, and the [iDREM program and resources](https://github.com/phoenixding/idrem). The checked Python environment used Python 3.10.20 and scUNAGI 0.5.10.

```sh
python train_unagi.py STAGE_DIR IDREM_DIR SEED
```

Use a fresh stage directory with copies of `0.h5ad` (healthy), `1.h5ad` (MGUS), `2.h5ad` (SMM), and `3.h5ad` (MM). These study inputs share 30,232 genes in the same order and contain the observation fields `stage` and `marker_driven_celltype`. The underlying single-cell datasets are [GSE124310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124310) and [GSE271107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271107).

`IDREM_DIR` is the software directory containing `idrem.jar`, `example_settings.txt`, and the companion resources, not a previous results directory. Use paths without spaces. Pass an integer seed for the new run; the original run-specific seed was not recorded.

Training runs on CPU for five outer iterations. UNAGI calls iDREM within each iteration. Model checkpoints are written to `STAGE_DIR/model_save/` and final-iteration iDREM output to `STAGE_DIR/4/idremResults/`.

## License and data

The [MIT License](LICENSE) covers the original scripts and documentation in this repository. Third-party software and datasets are not included and retain their own licenses or access conditions. This repository does not distribute patient-level data, summary-statistic files, pretrained models, or LD panels.
