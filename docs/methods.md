# Current methods: exploratory tissue-group comparison

## Scope and question

What regional patterns do predefined rod- and cone-associated genes show in the published macular and peripheral retina group FPKM estimates from GSE94437? This replaces the proposed eight-pair inferential analysis because the available processed expression file lacks individual sample measurements. Sample metadata still documents the original 16 retinal libraries from eight eye IDs; these are not 16 expression columns in the current analysis.

## Data preparation

Run `python3 scripts/prepare_uploaded_retina.py`. The original gzip is preserved. Extract both neural-retina groups, retaining tracking IDs, gene annotation, original confidence bounds and estimation status. Exclude RPE/choroid/sclera columns. Do not infer replicate-level expression from group estimates or bounds. The upload audit contains all original column names and file hashes.

## Marker matching and QC

Run `python3 scripts/explore_groups.py`. Use the existing frozen panel in data/gene_sets.csv (10 rod-associated and 11 cone-associated genes), without selecting genes based on observed regional direction. Match gene_short_name exactly. Missing and nonunique matches remain in the output audit and are excluded from plots; no ambiguous entries are summed. Validate finite, nonnegative FPKM. Retain only unique matches with OK status in both regions in plots; flag FPKM below 1 without removing them. All 21 panel genes currently match uniquely, have OK status, and exceed FPKM 1 in both groups.

## Descriptive measures

The heatmap displays log2(FPKM + 1) on one shared scale. It shows group expression magnitudes, without row z scores or claims about cell abundance. The ratio chart displays log2(Macular FPKM / Peripheral FPKM), without a pseudocount. Ratios involving a zero are left undefined and omitted from that chart. The CSV separately provides log2(FPKM + 1) differences; these are distinct from the raw FPKM ratios. Positive ratios mean higher macular group estimates; negative ratios mean higher peripheral group estimates.

Do not compute paired tests, FDR, donor-level confidence intervals, standardized paired effects, donor consistency, sample PCA, or sample correlations. Genes are not biological replicates. Original Cuffdiff interval fields remain available as source annotations and are not converted into regional contrast intervals.

## Interpretation and limits

These are descriptive patterns in published tissue-group estimates, not statistically significant donor-level differences. Bulk tissue expression alone does not determine rod/cone cell counts, and FPKM comparisons can be affected by library composition and the original quantification procedure. Close ratios should not be treated as robust directionality. Opsin paralogs also warrant caution about short-read mapping specificity. The observed direction may contradict the initial hypothesis; report it without changing the marker panel. Further interpretation requires sample-level data and source-study context.

The original paired workflow remains available in scripts/analyze.py for a future verified sample-level input. It is not the current analysis. Its former plan is archived in docs/archive/paired-methods.md.

## Verification and post-hoc descriptive checks (20 November 2026)

Run scripts/verify_exploratory.py after generating the results. This separately reads the original gzip and verifies all retained fields, all 21 marker values and the regional log ratios. GEO sample titles are cross-checked against tissue-characteristic fields. No evidence of a column swap was found in the local processing; upstream biological labeling cannot be established by these checks.

Three opsins (OPN1LW, OPN1MW, OPN1SW) are excluded only in an additional sensitivity summary; the main panel remains unchanged. The separate log2(FPKM + 1) contrast preserves all signs, as mathematically expected. A post-hoc magnitude band from 1/1.1 to 1.1 flags four near-equal gene ratios; it is not an equivalence or significance test. Panel medians summarize gene ratios, not donor effects. See results/exploratory/verification.json for complete outputs and docs/source-context.md for source checks.
