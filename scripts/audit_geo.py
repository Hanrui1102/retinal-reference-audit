"""Download GEO records, verify the eight pairs, and audit expression granularity."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from shutil import copyfileobj
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94437/"
FILES = {
    "matrix/GSE94437_series_matrix.txt.gz": "GSE94437_series_matrix.txt.gz",
    "suppl/GSE94437_genes.fpkm_tracking.gz": "GSE94437_genes.fpkm_tracking.gz",
}


def audit(download):
    raw = ROOT / "data/raw"
    raw.mkdir(parents=True, exist_ok=True)
    paths = [raw / name for name in FILES.values()]
    for remote, path in zip(FILES, paths, strict=True):
        if not path.exists():
            if not download:
                raise SystemExit(f"Missing {path.name}. Run with --download or see docs/data-access.md.")
            temporary = path.with_suffix(".part")
            with urlopen(BASE + remote, timeout=30) as response, temporary.open("wb") as out:
                copyfileobj(response, out)
            temporary.replace(path)

    with gzip.open(paths[0], "rt", encoding="utf-8") as handle:
        matrix = list(csv.reader(handle, delimiter="\t"))
    fields = {row[0]: row[1:] for row in matrix if row and row[0].startswith("!Sample_")}
    entries = list(zip(fields["!Sample_geo_accession"], fields["!Sample_title"], strict=True))
    if len(entries) != 32 or len({gsm for gsm, _ in entries}) != 32:
        raise ValueError("Expected 32 uniquely labeled GEO samples.")
    metadata = []
    for gsm, title in entries:
        match = re.fullmatch(r"(\d+-\d+) (Macular|Peripheral) Retina", title)
        if match:
            eye, region = match.groups()
            metadata.append((gsm, eye, region.lower(), title))
    metadata.sort(key=lambda row: (row[1], row[2]))
    pairs = Counter((eye, region) for _, eye, region, _ in metadata)
    eyes = sorted({eye for _, eye, _, _ in metadata})
    if len(eyes) != 8 or set(pairs) != {(eye, region) for eye in eyes
                                       for region in ("macular", "peripheral")} or set(pairs.values()) != {1}:
        raise ValueError("Expected eight complete, unique retinal pairs.")

    with gzip.open(paths[1], "rt", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        columns = next(rows)
        n_rows = sum(bool(row) for row in rows)
    start = matrix.index(["!series_matrix_table_begin"])
    end = matrix.index(["!series_matrix_table_end"])
    matrix_rows = matrix[start+2:end]
    required = {row[0] for row in metadata}
    observed = sorted(required.intersection(columns))
    report = {
        "total_geo_samples": len(entries), "retina_samples": len(metadata),
        "complete_eye_pairs": len(eyes), "series_matrix_expression_rows": len(matrix_rows),
        "tracking_rows": n_rows, "tracking_fpkm_columns": [c for c in columns if c.endswith("_FPKM")],
        "tracking_retinal_gsm_columns": observed,
        "per_sample_expression_available_in_downloads": len(observed) == 16 or bool(matrix_rows),
        "sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    }
    with (ROOT / "data/metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "eye_id", "region", "geo_title"])
        writer.writerows(metadata)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results/geo_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Download missing GEO files.")
    audit(parser.parse_args().download)
