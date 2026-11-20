"""Build the six-page English research report from verified project results."""
import csv
import json
from pathlib import Path
from xml.sax.saxutils import escape
import reportlab
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image, Table, TableStyle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/pdf/retina_exploratory_report.pdf'
FONT_DIR=Path(reportlab.__file__).parent/'fonts'
for name,file in [('Vera','Vera.ttf'),('Vera-Bold','VeraBd.ttf'),('Vera-Italic','VeraIt.ttf'),('Vera-BoldItalic','VeraBI.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(FONT_DIR/file)))
pdfmetrics.registerFontFamily('Vera',normal='Vera',bold='Vera-Bold',italic='Vera-Italic',boldItalic='Vera-BoldItalic')
INK=colors.HexColor('#19313b')
ACCENT=colors.HexColor('#276678')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='ReportTitle',fontName='Vera-Bold',fontSize=24,leading=28,textColor=INK,spaceAfter=17))
styles.add(ParagraphStyle(name='Section',fontName='Vera-Bold',fontSize=14,leading=18,textColor=INK,spaceBefore=13,spaceAfter=8))
styles.add(ParagraphStyle(name='Text',fontName='Vera',fontSize=10,leading=14,spaceAfter=10,textColor=INK))
styles.add(ParagraphStyle(name='CaptionText',fontName='Vera',fontSize=9,leading=12,textColor=INK,spaceAfter=9))
styles.add(ParagraphStyle(name='Eyebrow',fontName='Vera-Bold',fontSize=9,leading=12,textColor=ACCENT,spaceAfter=12))
styles.add(ParagraphStyle(name='RefText',fontName='Vera',fontSize=8.4,leading=12,spaceAfter=8,textColor=INK))


def para(text,style='Text'):
    return Paragraph(text,styles[style])

def heading(text):
    return para(text,'Section')

def footer(canvas,doc):
    canvas.saveState()
    w,h=doc.pagesize
    canvas.setStrokeColor(colors.HexColor('#ced8dc'))
    canvas.line(48,41,w-48,41)
    canvas.setFillColor(INK);canvas.setFont('Vera',8)
    canvas.drawString(48,28,'GSE94437 | Exploratory secondary analysis | 20 November 2026')
    canvas.drawRightString(w-48,28,str(doc.page))
    canvas.restoreState()

def main():
    audit=json.loads((ROOT/'results/exploratory/verification.json').read_text())
    assert audit['raw_to_cleaned_all_cells_equal'] and audit['raw_to_marker_values_and_ratios_verified']==21
    with (ROOT/'results/exploratory/marker_expression.csv').open() as f: rows=list(csv.DictReader(f))
    by_symbol={r['symbol']:r for r in rows}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    story=[]
    def add(*items): story.extend(items)
    add(para('UNDERGRADUATE COMPUTATIONAL BIOLOGY PROJECT','Eyebrow'),
        para('Rod/cone-associated expression across human retinal regions','ReportTitle'),
        para('A reproducible exploratory comparison of published tissue-group FPKM estimates from GSE94437.'),
        para('Version 1.0 | Research report | 20 November 2026','CaptionText'),
        heading('Abstract'),
        para('<b>Objective.</b> Describe regional patterns in an existing panel of rod- and cone-associated genes without treating group estimates as individual biological replicates.'),
        para('<b>Methods.</b> The original compressed tracking file was audited directly. Its two neural-retina FPKM columns were extracted, preserving gene identifiers and source estimation status. Twenty-one predefined genes were matched by exact symbol. Expression was visualized as log2(FPKM + 1), and regional contrasts as log2(macular FPKM / peripheral FPKM). A separate verification script compared the cleaned data and marker values back to the original file.'),
        para('<b>Results.</b> All 21 genes matched uniquely and had OK status in both groups. All 10 rod-associated genes had higher macular estimates. Ten of 11 cone-associated genes had higher peripheral estimates; OPN1LW was nearly equal, with a macular/peripheral ratio of 1.021. Excluding all three opsins left eight cone-associated genes, all higher in the peripheral group.'),
        para('<b>Conclusion.</b> The supplied group values do not show the initially expected combination of peripheral rod enrichment and macular cone enrichment. This is a descriptive finding about a deposited processed table, not evidence of a statistically significant population-level reversal, altered cell counts, or a mechanism.'),
        heading('Research question and scope'),
        para('What regional expression patterns are visible in the published group estimates for a fixed rod/cone-associated panel? The original plan called for eight eye-level paired comparisons. That plan was revised when inspection showed that individual-sample expression was unavailable in the supplied file. This revision and all subsequent sensitivity checks are exploratory, not preregistered.'))
    add(PageBreak(),para('01 / DATA AND METHODS','Eyebrow'),heading('Source data and provenance'),
        para('GSE94437 is the official GEO series for the source study [1,2]. The local Series Matrix contains 32 sample records. Sixteen have neural-retina titles, yielding eight complete eye-ID pairs, each with one macular and one peripheral record. Titles and tissue-characteristic fields agree. The Series Matrix has zero expression rows.'),
        para('The actual genes.fpkm_tracking.gz file contains 57,815 tracking records and 25 fields: nine annotation fields plus four tissue groups, each with FPKM, lower bound, upper bound and status. Only Macular_Retina_FPKM and Peripheral_Retina_FPKM enter the comparison. The other tissue groups are excluded. No individual retinal GSM expression columns are present.'),
        para('The source SHA-256 is recorded in the project audits. Gene IDs, symbols, loci, confidence bounds and status fields are retained in the cleaned table. Group estimates are not assumed to be simple donor averages. No missing expression values are filled, no duplicate symbols are merged, and no source values are re-normalized.'),
        heading('Fixed marker panel and quality checks'),
        para('The main panel contains 10 rod-associated and 11 cone-associated genes, retained from the project plan before the first regional marker extraction. Its rationale and literature references are stored in data/gene_sets.csv. Primary phototransduction and single-cell studies support pathway or cell-type associations [3,4]. These associations are not claims of exclusive expression in one cell type.'),
        para('An exact gene_short_name match is required. Multiple matches or absent matches would be reported and excluded from plotting. Every current marker matches exactly one tracking record, is finite and nonnegative, exceeds FPKM 1 in both groups, and has source status OK. Source status does not replace sample-level QC.'),
        heading('Transforms and units'),
        para('Heatmap values are log2(FPKM + 1), using one shared color scale. The contrast is log2(M/P), where M and P denote macular and peripheral FPKM. A value of +1 corresponds to twice the macular estimate; -1 to half. Ratios involving zeros would be undefined and omitted, but none occur in this panel. The results CSV also stores the distinct quantity log2(M + 1) - log2(P + 1).'),
        para('No paired tests, p-values, FDR, donor confidence intervals, donor consistency, sample PCA or sample correlations are estimated. Marker genes are not treated as independent donor replicates.'))
    add(PageBreak(),para('02 / REGIONAL CONTRASTS','Eyebrow'),heading('Direction and magnitude of group estimates'),
        Image(str(ROOT/'figures/group_regional_ratios.png'),width=6.1*inch,height=6.8625*inch),
        para('<b>Figure 1.</b> Regional ratios for the fixed panel. Positive values mean higher macular FPKM; negative values mean higher peripheral FPKM. Colors indicate rod-associated and cone-associated panel membership. Bars are genes, not replicates. No uncertainty interval for this contrast is available.','CaptionText'),
        para('The rod-associated panel is uniformly positive in this file. Cone-associated contrasts are mostly negative, often modest in size. OPN1LW is slightly positive, so its direction should not be emphasized as robust.'))
    add(PageBreak(),para('03 / EXPRESSION MAGNITUDES','Eyebrow'),heading('A shared scale preserves expression differences'),
        Image(str(ROOT/'figures/group_expression_heatmap.png'),width=5.1*inch,height=6.5571*inch),
        para('<b>Figure 2.</b> Both retinal group estimates on a shared log2(FPKM + 1) scale. Rows are not standardized separately. The plot displays measured transcript abundance estimates and does not constitute a rod/cone abundance score.','CaptionText'),
        para('For example, RHO is highly expressed in both groups, whereas several cone-associated genes have much lower FPKM. Cross-gene magnitude is shown for context; neither high FPKM nor a regional ratio directly reports the number of cells expressing that gene.'))
    add(PageBreak(),para('04 / VERIFICATION AND INTERPRETATION','Eyebrow'),heading('Independent calculation checks'),
        para('The verification script reads the original gzip separately from the plotting workflow. It checks every retained annotation/expression field against the cleaned file, checks all 21 marker measurements, and recalculates each log ratio as log2(M) - log2(P). It also verifies agreement between GEO sample titles and tissue characteristics and confirms eight complete eye-ID pairs. These checks pass. They establish internal processing consistency, not external biological validation.'),
        heading('Post-hoc descriptive sensitivity checks'),
        para('Removing OPN1LW, OPN1MW and OPN1SW leaves all eight remaining cone-associated genes higher in the peripheral group. The 10 rod genes are unchanged. This shows that the panel direction is not entirely driven by the three opsin entries; it does not validate the upstream quantification.'),
        para('Replacing log2(M/P) by log2(M + 1) - log2(P + 1) preserves every marker direction. Direction preservation is mathematically expected for a monotonic transform of positive values and is a consistency check, not independent evidence. The original ratio remains the primary descriptive measure.'),
        para('An illustrative, post-hoc near-equality band of 1/1.1 to 1.1 contains OPN1LW, GNGT2, PDE6H and ARR3. This is a magnitude aid, not a biological equivalence or significance threshold. Median gene ratios are 1.748 for the rod panel and 0.868 for the cone panel; these summarize genes, not a donor distribution.'),
        heading('Why the initial hypothesis was too simple'),
        para('The source paper describes broad biopsy sampling and notes that macular tissue contains many rods [2]. A macular biopsy must not be equated with a pure cone-rich foveal sample. Primary anatomical mapping shows strong spatial variation in photoreceptor density [5]. Sampling extent therefore matters when forming regional expectations.'),
        para('This context weakens the simple expectation that every cone-associated transcript must be higher in a macular tissue sample. It does not establish the cause of the observed ratios. Tissue composition, within-cell expression and original normalization remain possible contributors that cannot be separated using two group estimates.'))
    add(PageBreak(),para('05 / CONCLUSIONS AND REPRODUCIBILITY','Eyebrow'),heading('What this project establishes'),
        para('A documented and reproducible descriptive pattern is present in the supplied GSE94437 group table. All rod-panel estimates are higher in the macular column; most cone-panel estimates are higher in the peripheral column. The main panel was kept intact despite disagreement with the initial directional expectation. The result is suitable as an undergraduate secondary-analysis case study, without a claim of novel biological discovery.'),
        heading('What remains unresolved'),
        para('Donor variability and statistical reliability cannot be measured. Internal label checks do not independently establish upstream tissue identities. The source paper reports RefSeq quantification, while the deposited table uses versioned Ensembl identifiers [2]; exact annotation and processing equivalence is unresolved. The original study findings should not be presented as numerically reproduced here. Bulk RNA-seq ratios cannot separate cell proportions from cell-specific regulation.'),
        para('The primary single-cell atlas reports difficulty separating highly similar L/M opsins in its own assay [4]. This motivates caution and the opsin-exclusion check, but does not prove an opsin mapping error in this dataset. No independent expression dataset or laboratory experiment was used to validate the current observation.'),
        heading('Reproduce and inspect'),
        para('From the project root: install requirements-exploratory.txt; run scripts/audit_geo.py --download, scripts/prepare_uploaded_retina.py, scripts/explore_groups.py and scripts/verify_exploratory.py. Full inputs and SHA-256 hashes are documented. Figures and the 21-gene results are regenerated locally. scripts/build_report.py rebuilds this report with the optional report dependencies. See README.md and docs/methods.md.'),
        heading('References'),
        para('[1] NCBI GEO. GSE94437. <link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94437" color="#276678">Official series record</link>. Accessed 20 November 2026.','RefText'),
        para('[2] Li M et al. Comprehensive analysis of gene expression in human retina and supporting tissues. Human Molecular Genetics. 2014;23:4001-4014. <link href="https://doi.org/10.1093/hmg/ddu114" color="#276678">doi:10.1093/hmg/ddu114</link>.','RefText'),
        para('[3] Christopher M et al. Selection of Phototransduction Genes in Homo sapiens. IOVS. 2013;54:5489-5496. <link href="https://doi.org/10.1167/iovs.12-11454" color="#276678">doi:10.1167/iovs.12-11454</link>.','RefText'),
        para('[4] Lukowski SW et al. A single-cell transcriptome atlas of the adult human retina. EMBO Journal. 2019;38:e100811. <link href="https://doi.org/10.15252/embj.2018100811" color="#276678">doi:10.15252/embj.2018100811</link>.','RefText'),
        para('[5] Curcio CA et al. Human photoreceptor topography. Journal of Comparative Neurology. 1990;292:497-523. <link href="https://doi.org/10.1002/cne.902920402" color="#276678">doi:10.1002/cne.902920402</link>.','RefText'))
    doc=SimpleDocTemplate(str(OUT),pagesize=(612,792),rightMargin=48,leftMargin=48,topMargin=48,bottomMargin=54,
                          title='Rod/cone-associated expression across human retinal regions',author='Independent research project')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    print(OUT)

if __name__=='__main__': main()
