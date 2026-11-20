"""Independent raw-to-result verification and post-hoc descriptive sensitivity checks."""
import csv
import gzip
import hashlib
import json
import math
import platform
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))

def main():
    raw = ROOT / 'data/raw/GSE94437_genes.fpkm_tracking.gz'
    with gzip.open(raw, 'rt') as handle:
        source = list(csv.DictReader(handle, delimiter='\t'))
    raw_by_id = {row['tracking_id']: row for row in source}
    assert len(raw_by_id)==len(source)
    cleaned = read_csv(ROOT / 'data/retinal_group_expression.csv')
    assert len(cleaned)==len(source) and len({r['tracking_id'] for r in cleaned})==len(source)
    for row in cleaned:
        assert all(raw_by_id[row['tracking_id']][key]==value for key,value in row.items())
    assert len(raw_by_id) == len(source) == len(cleaned) == 57815
    for row in cleaned:
        assert all(raw_by_id[row['tracking_id']][key] == value for key, value in row.items())
    results = read_csv(ROOT / 'results/exploratory/marker_expression.csv')
    panel = read_csv(ROOT / 'data/gene_sets.csv')
    assert [r['symbol'] for r in results] == [r['symbol'] for r in panel]
    for row in results:
        matches = [r for r in source if r['gene_short_name'] == row['symbol']]
        assert len(matches)==1
        original = matches[0]
        m,p = [float(original[f'{region}_Retina_FPKM']) for region in ('Macular','Peripheral')]
        assert float(row['macular_fpkm']) == m and float(row['peripheral_fpkm']) == p
        assert math.isclose(float(row['log2_macular_over_peripheral']), math.log2(m)-math.log2(p), abs_tol=1e-12)
        assert math.isclose(float(row['log2_plus1_difference']), math.log2((m+1)/(p+1)), abs_tol=1e-12)
    with gzip.open(ROOT / 'data/raw/GSE94437_series_matrix.txt.gz','rt') as handle:
        matrix = list(csv.reader(handle,delimiter='\t'))
    accession = next(r[1:] for r in matrix if r and r[0]=='!Sample_geo_accession')
    titles = next(r[1:] for r in matrix if r and r[0]=='!Sample_title')
    tissues = next(r[1:] for r in matrix if r and r[0]=='!Sample_characteristics_ch1' and r[1].startswith('tissue:'))
    mapping = dict(zip(accession,zip(titles,tissues),strict=True))
    metadata = read_csv(ROOT / 'data/metadata.csv')
    assert len(metadata)==16 and len({r['sample_id'] for r in metadata})==16
    for row in metadata:
        title,tissue = mapping[row['sample_id']]
        assert row['geo_title']==title==f"{row['eye_id']} {row['region'].title()} Retina"
        assert tissue == f"tissue: {row['region'].title()} Retina"
    for eye in {r['eye_id'] for r in metadata}:
        assert sorted(r['region'] for r in metadata if r['eye_id']==eye)==['macular','peripheral']
    sensitivity = []
    for label, subset in [('full_panel', results), ('excluding_three_opsins', [r for r in results if not r['symbol'].startswith('OPN1')])]:
        for group in ('rod','cone'):
            rows = [r for r in subset if r['gene_set']==group]
            ratios = [float(r['macular_over_peripheral']) for r in rows]
            sensitivity.append({'panel':label,'association':group,'genes':len(rows),
                'higher_macular':sum(x>1 for x in ratios),'higher_peripheral':sum(x<1 for x in ratios),
                'median_macular_over_peripheral':statistics.median(ratios),
                'direction_unchanged_with_plus1':all((float(r['log2_macular_over_peripheral'])>0)==(float(r['log2_plus1_difference'])>0) for r in rows),
                'within_posthoc_1_over_1_1_to_1_1_band':[r['symbol'] for r in rows if 1/1.1<=float(r['macular_over_peripheral'])<=1.1]})
    target = ROOT / 'results/exploratory'
    report = {'raw_to_cleaned_all_cells_equal':True,'raw_to_marker_values_and_ratios_verified':len(results),
              'sample_titles_and_tissue_fields_consistent':True,'complete_eye_pairs':8,
              'raw_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
              'python_version':platform.python_version(), 'sensitivity':sensitivity,
              'interpretation':'Post-hoc descriptive checks, not statistical validation. Median ratios summarize genes, not donors. The 1/1.1 to 1.1 band is an illustrative magnitude band, not a significance threshold.',
              'provenance_limit':'File labels and metadata are internally consistent. Upstream biological label correctness and exact quantification cannot be independently verified from this file.'}
    (target / 'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
