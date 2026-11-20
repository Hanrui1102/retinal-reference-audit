> Historical data-access audit. The project now uses a descriptive tissue-group comparison; see [current results](exploratory/summary.md). The missing sample-level input still prevents the original paired analysis.

# Data verification summary

The GSE94437 records identify 32 biological samples. Exact retinal title matching retains 16 neural-retina libraries from eight eyes. Each eye has one macular and one peripheral library. For example, eye `11-1516` maps to macular `GSM2475287` and peripheral `GSM2475285`. See the [16-sample metadata](../data/metadata.csv) and [machine-readable audit](geo_audit.json).

The downloaded GEO FPKM file has 57,815 tracking rows and four group-level FPKM fields, two for retina and two for RPE/choroid/sclera. It has zero of the 16 retinal GSM columns. The downloaded Series Matrix has 32 GSM headers and zero expression rows. The original publication's supplementary tables cover sample mapping and sequencing metrics, enrichment, selected group summaries and primer information. They contain no gene-by-library FPKM matrix. [GEO series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94437) · [Li et al., with supplementary materials](https://pmc.ncbi.nlm.nih.gov/articles/PMC7297232/).

**Paired biological results are pending.** A verified GSM-level input is required for expression QC, tests, effect sizes and figures. The workflow ran end-to-end with temporary synthetic values, and the GEO group file was rejected as analysis input. No synthetic output is presented as a research result.
