"""Inspect the original upload and extract group-level retina data, never donor replicates."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from audit_geo import ROOT, audit


def prepare(source):
    source = Path(source)
    canonical = ROOT / 'data/raw/GSE94437_genes.fpkm_tracking.gz'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != hashlib.sha256(canonical.read_bytes()).hexdigest():
        raise ValueError('Uploaded file differs from audited GEO input; review before proceeding.')
    audit(False)
    with gzip.open(source, 'rt', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        columns = reader.fieldnames
        rows = list(reader)
    annotation = ['tracking_id', 'gene_id', 'gene_short_name', 'locus']
    measures = [f'{region}_Retina_{field}' for region in ('Macular', 'Peripheral')
                for field in ('FPKM', 'conf_lo', 'conf_hi', 'status')]
    with (ROOT / 'data/retinal_group_expression.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=annotation + measures)
        writer.writeheader()
        writer.writerows({key: row[key] for key in annotation + measures} for row in rows)
    with (ROOT / 'data/metadata.csv').open() as handle:
        metadata = list(csv.DictReader(handle))
    pair_counts = Counter((row['eye_id'], row['region']) for row in metadata)
    report = {
        'uploaded_filename': source.name, 'sha256': digest,
        'identical_to_local_geo_download': True,
        'columns': columns, 'column_count': len(columns), 'tracking_rows': len(rows),
        'retinal_samples_in_metadata': len(metadata),
        'complete_pairs': len({row['eye_id'] for row in metadata}),
        'one_sample_per_eye_region': all(n == 1 for n in pair_counts.values()),
        'sample_expression_columns': [r['sample_id'] for r in metadata if r['sample_id'] in columns],
        'retina_fpkm_columns': [c for c in measures if c.endswith('_FPKM')],
        'missing_retina_fpkm_cells': sum(row[c] in ('', '-', 'NA', 'nan') for row in rows for c in measures if c.endswith('_FPKM')),
        'retina_status_counts': {region: dict(Counter(row[f'{region}_Retina_status'] for row in rows)) for region in ('Macular', 'Peripheral')},
        'unique_tracking_ids': len({row['tracking_id'] for row in rows}),
        'paired_analysis_ready': False,
        'limitation': 'Two tissue-group estimates do not provide sixteen individual sample measurements. Confidence bounds are not donor replicates.'
    }
    for row in rows:
        for col in report['retina_fpkm_columns']:
            value = float(row[col])
            if not 0 <= value < float('inf'):
                raise ValueError(f'Invalid FPKM: {row["tracking_id"]}, {col}')
    (ROOT / 'results/upload_inspection.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'columns'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default=str(ROOT / 'data/raw/GSE94437_genes.fpkm_tracking.gz'))
    prepare(parser.parse_args().source)
