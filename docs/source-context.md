# Source context and evidence boundaries

Checked 20 November 2026. This note supports interpretation, not an independent replication.

## Evidence used

- [GEO GSE94437](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94437): official study record and sample identities. The local Series Matrix contains 32 sample labels, including 16 neural-retina libraries, but zero expression rows. Titles and tissue-characteristic fields agree for the retained samples.
- [Li et al. 2014](https://doi.org/10.1093/hmg/ddu114), Methods and Discussion: broad biopsy collection and anatomical heterogeneity caution against equating the macular sample with the cone-dominated foveal center. The publication describes RefSeq quantification, whereas the deposited tracking file uses versioned Ensembl IDs. These outputs should not be assumed numerically identical without the original annotation and analysis files.
- [Curcio et al. 1990, Human photoreceptor topography](https://pubmed.ncbi.nlm.nih.gov/2324310/): primary anatomical mapping demonstrates strongly varying rod/cone densities across retinal locations. Anatomical expectations depend on sampling extent.
- [Christopher et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3743458/): primary phototransduction study used for pathway membership in the frozen panel; pathway membership is not proof of exclusive cell-type expression.
- [Lukowski et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6745503/): primary single-cell study supports selected rod/cone marker associations. Its warning about highly homologous opsins motivates an additional descriptive check excluding all three opsins. The 3-prime assay limitation in that paper does not itself demonstrate misassignment in GSE94437.

The frozen panel also retains a 2025 review as its original supplementary functional reference. It is not an independent validation dataset. No marker was added or removed from the main analysis during this review.

## What was verified locally

scripts/verify_exploratory.py independently reads the original compressed tracking table, compares every retained expression/annotation field to the cleaned table, checks all 21 marker FPKM values and recomputes the log ratios as a difference of logs. It cross-checks sample titles against tissue characteristics. It records descriptive checks excluding opsins and using a +1 transformation.

## What remains unknown

No BAM/FASTQ or independent sample-level matrix was used. Internal label consistency does not establish the biological correctness of all upstream labels. Group estimates are not asserted to be arithmetic donor means. The annotation-version provenance of the deposited tracking output is incomplete. Neither the observed regional directions nor any explanation for them have been independently replicated here.
