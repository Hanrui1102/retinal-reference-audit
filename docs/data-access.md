> Current project: the published tissue-group file is sufficient for the revised descriptive comparison. Run `python3 scripts/prepare_uploaded_retina.py` then `python3 scripts/explore_groups.py`. The sample-level input described below is only required to restore the former paired analysis.

# Data access and required input

## GEO files you can reproduce today

- [GSE94437 series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94437)
- [Series Matrix download](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94437/matrix/GSE94437_series_matrix.txt.gz)
- [Cuffdiff FPKM tracking download](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94437/suppl/GSE94437_genes.fpkm_tracking.gz)

Run `python3 scripts/audit_geo.py --download`, or download the two files manually, keep their original names and place them in `data/raw/`. Then run `python3 scripts/audit_geo.py`.

The audit pairs samples by the eye ID and exact tissue label in the GEO sample titles. The 16 supporting tissue samples consist of RPE/choroid/sclera and are outside this study. The two retinal tissue-group estimates are `Macular_Retina_FPKM` and `Peripheral_Retina_FPKM`. A single group estimate cannot supply the eight within-eye differences required by the proposal.

## Obtain sample-level FPKM

The missing input is a gene-level FPKM matrix for the **16 retinal GSM accessions** in [metadata.csv](../data/metadata.csv). Ask the original study authors or the research supervisor whether the original per-library quantification or Cuffdiff replicate-level tracking output is available. Request all measured genes for sample QC and PCA, annotation build/version, gene IDs, and the library label to GSM correspondence. If provided with internal labels, use the authors' mapping and validate against the eye and region labels. [Original paper and supplementary material](https://pmc.ncbi.nlm.nih.gov/articles/PMC7297232/).

GEO also links raw sequencing reads. Using these to produce new FPKM values would require a revised data-acquisition method and a documented quantification workflow. RNA-seq counts, TPM, normalized compendium scores and inferred replicate values are outside the current FPKM input specification.

## Format after obtaining the matrix

- Save as CSV or TSV, with optional gzip compression. The first field is `gene_symbol` and the remaining 16 headers are the verified GSM accessions. The script aligns sample columns by ID, independent of their input order.
- Use one row per official gene symbol and numeric, nonnegative FPKM. Use an empty cell for missing values and zero only for a reported zero. Resolve duplicated and ambiguous gene identifiers with the source annotation before running. Preserve the raw file and annotation separately.
- Supply provenance in `--provenance`. Include your source/version, confirmation of gene-level FPKM, and how sample identifiers were matched. The script records a checksum for the input and validates its structure. The units and mapping require source documentation.
- Keep input files in `data/raw/`. Generated results go in a new output directory. The source matrix is excluded from Git.
