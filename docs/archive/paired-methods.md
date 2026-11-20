# Methods and figure guide

## 1. Pairing and quality control

Use eight donor eyes as the units of analysis. Each contributes matched macular and peripheral neural retina. Transform every observed FPKM by log₂(FPKM + 1). Flag missing values, duplicate symbols, zeros and low expression (FPKM below 1). This threshold is descriptive. Duplicate IDs stop the workflow. A frozen marker with a missing value in any sample is documented and omitted from the complete-pair calculations. All-zero and constant markers remain visible in QC and descriptive gene-level contrasts. Constant markers have undefined signature z scores and are omitted from signatures and the heatmap.

Check the sample distributions, Spearman correlations and exploratory PCA before interpreting any tests. The PCA uses all complete, variable genes in the supplied matrix, centered across samples on the transformed scale. Review unusually low correlation or extreme scores with the underlying sample information. The workflow never drops an eye automatically.

## 2. Signatures and paired tests

Freeze the panel before testing. It contains 10 rod-associated and 11 cone-associated genes. Inclusion is based on subtype-associated phototransduction function or established marker evidence in the cited literature. Each row stores its rationale and source. The terms rod-associated and cone-associated allow overlapping expression between cell types.

For each complete, variable marker, compute a z score across all 16 samples using the sample standard deviation. For each sample, average its marker z scores within each gene set, giving each included gene equal weight. Report the genes retained in each signature.

For each gene and signature calculate eight differences, each defined as macular minus peripheral. Positive values indicate higher macular expression or signature score. The mean difference has units of log₂(FPKM + 1) for genes and mean z-score units for signatures. Cohen's paired d_z divides that mean by the standard deviation of the eight differences. The t-based 95% confidence intervals are pointwise.

- **Primary comparison.** A two-sided paired t-test, implemented as a one-sample t-test of the differences against zero. Inspect the individual differences for asymmetry and influential eyes. The eight eyes are assumed independent and the t-test assumes approximately normal differences. At this sample size, a normality test has limited ability to diagnose departures.
- **Sensitivity check.** Two-sided Wilcoxon signed-rank test. The workflow excludes zero differences, uses average ranks for absolute ties and enumerates all sign assignments for the nonzero ranks. Differences are rounded to 12 decimals for tie handling. A fully zero difference vector receives sensitivity p = 1. The sign-flip calculation assumes symmetry of differences under the null.
- **Multiplicity.** Apply Benjamini–Hochberg to the primary t-test p-values across the full 21-gene family. A missing or undefined test occupies a slot with p = 1 during adjustment and displays blank in the result. Signature tests and Wilcoxon sensitivity p-values are reported as nominal.
- **Consistency.** Read the count of eyes in either direction and the difference SD alongside the mean and interval. Where the difference vector is constant, the workflow reports the descriptive mean and leaves the t-test, interval and d_z blank.

## 3. Read the figures

The workflow generates PNG (300 dpi) and SVG copies of these figures when the required input is available.

1. **Paired signatures.** Each line connects the two samples from one eye. Look for shifts that repeat across lines.
2. **Marker heatmap.** Each row has its own standardized color scale values on a common z-score legend. Adjacent M and P columns belong to one eye.
3. **Gene effects.** Points show mean paired differences and bars show the 95% intervals when estimable. Zero marks equal mean regional expression on the transformed scale.
4. **Exploratory PCA.** Look for sample structure and check possible outliers. Region colors and donor connectors are descriptive. At least two complete, variable genes are required.
5. **Sample distributions.** Compare the medians and spreads of transformed expression across samples.
6. **Sample correlations.** Look for samples that have unusually low rank correlation with the rest of the set.

## 4. Build the research story once expression is available

Suggested 5–8 page report structure: background and hypothesis, data provenance and pairing, prespecified markers and methods, sample QC, paired results, biological interpretation and limitations, references. The one-page abstract can use Objective, Methods, Results and Conclusion. For a poster or short talk, start with the paired study design, show the signature and effect-size figures, then discuss uncertainty and the interpretation of bulk RNA-seq.

Eight eyes give limited precision and coverage of donor variation. Bulk tissue expression reflects a mixture of cell types and levels of gene expression within those cells. Regional association alone cannot determine their relative contribution. Keep the interpretation anchored to the sampled tissue and observed paired differences.
