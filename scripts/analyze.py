"""Analyze a verified GSE94437 gene-by-sample FPKM matrix."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
COLORS = {"rod": "#326899", "cone": "#AA5419", "macular": "#AA5419", "peripheral": "#326899"}
REGIONS = ["macular", "peripheral"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path, samples):
    delimiter = "," if path.name.removesuffix(".gz").endswith(".csv") else "\t"
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        headers = next(csv.reader(handle, delimiter=delimiter), [])
    if any(name.endswith("_Retina_FPKM") for name in headers):
        raise ValueError("GEO Cuffdiff file contains tissue-group estimates. Need 16 GSM-level FPKM columns.")
    if len(headers) != 17 or set(headers) != {"gene_symbol", *samples}:
        raise ValueError("Expected gene_symbol and exactly the 16 retinal GSM columns from data/metadata.csv.")
    frame = pd.read_csv(path, sep=delimiter, index_col="gene_symbol", encoding="utf-8-sig")
    if frame.empty or frame.index.hasnans or not frame.index.is_unique:
        raise ValueError("Provide nonempty, unique official gene symbols. Resolve ambiguous IDs at the source.")
    values = frame.loc[:, samples].apply(pd.to_numeric, errors="raise")
    if (values < 0).any().any() or np.isinf(values).any().any():
        raise ValueError("FPKM must be finite when observed and cannot be negative.")
    return values


def contrast(frame, pair_ids):
    macular = frame.loc[:, pair_ids.macular].to_numpy()
    peripheral = frame.loc[:, pair_ids.peripheral].to_numpy()
    return pd.DataFrame(macular - peripheral, index=frame.index, columns=pair_ids.index)


def wilcoxon_p(differences):
    # Enumerate all sign flips of nonzero differences, including tied ranks.
    nonzero = np.round(differences, 12)
    nonzero = nonzero[nonzero != 0]
    if not len(nonzero):
        return 1.0
    ranks = stats.rankdata(np.abs(nonzero))
    observed = abs(ranks @ np.sign(nonzero))
    sums = np.asarray(list(itertools.product((-1, 1), repeat=len(ranks)))) @ ranks
    return np.mean(np.abs(sums) >= observed)


def summarize(delta):
    records = []
    for name, row in delta.iterrows():
        values = row.to_numpy()
        mean = values.mean()
        sd = values.std(ddof=1)
        if sd > 0 and not np.all(values == values[0]):
            test = stats.ttest_1samp(values, 0)
            ci = test.confidence_interval()
            low, high, p_t, dz = ci.low, ci.high, test.pvalue, mean / sd
        else:
            low = high = p_t = dz = np.nan
        records.append({
            "feature": name, "n_pairs": len(values), "mean_delta": mean, "sd_delta": sd,
            "ci95_low": low, "ci95_high": high, "cohen_dz": dz, "p_paired_t": p_t,
            "p_wilcoxon": wilcoxon_p(values), "macular_higher": (values > 0).sum(),
            "peripheral_higher": (values < 0).sum(), "ties": (values == 0).sum(),
        })
    return pd.DataFrame(records).set_index("feature")


def save(fig, folder, name):
    fig.savefig(folder / f"{name}.png", dpi=300, facecolor="white")
    fig.savefig(folder / f"{name}.svg", facecolor="white")
    plt.close(fig)


def plot_results(log, z, signatures, effects, groups, metadata, pairs, folder):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "svg.fonttype": "none"})
    samples = metadata.index
    sample_labels = [f"{row.eye_id} {row.region[0].upper()}" for row in metadata.itertuples()]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True, layout="constrained")
    for ax, group in zip(axes, ["rod", "cone"], strict=True):
        for pair in pairs.itertuples():
            scores = signatures.loc[group, [pair.peripheral, pair.macular]]
            ax.plot([0, 1], scores, color=COLORS[group], marker="o", alpha=0.55)
        ax.set(xticks=[0, 1], xticklabels=["Peripheral", "Macular"], title=f"{group.title()} signature")
        ax.axhline(0, color="#777777", lw=0.6, ls="--")
    axes[0].set_ylabel("Mean within-gene z score")
    fig.suptitle("Photoreceptor signatures within donor eyes")
    save(fig, folder, "01_paired_signatures")

    fig, ax = plt.subplots(figsize=(10, 7.5), layout="constrained")
    values = z.loc[:, samples]
    limit = max(1, np.abs(values.to_numpy()).max())
    image = ax.imshow(values, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set(xticks=range(len(samples)), xticklabels=sample_labels, yticks=range(len(z)),
           yticklabels=[f"{g}  ({groups[g]})" for g in z.index], title="Photoreceptor gene expression across paired samples")
    ax.tick_params(axis="x", labelrotation=60)
    for boundary in np.arange(1.5, len(samples) - 0.5, 2):
        ax.axvline(boundary, color="white", lw=1.25)
    fig.colorbar(image, ax=ax, shrink=0.72, label="Within-gene z score")
    save(fig, folder, "02_marker_heatmap")

    fig, ax = plt.subplots(figsize=(8.5, 7.5), layout="constrained")
    for pos, (gene, row) in enumerate(effects.iterrows()):
        color = COLORS[groups[gene]]
        ax.plot(row.mean_delta, pos, "o", color=color)
        if pd.notna(row.ci95_low):
            ax.hlines(pos, row.ci95_low, row.ci95_high, color=color, lw=1.5)
    ax.set(yticks=range(len(effects)), yticklabels=effects.index,
           xlabel="Mean paired difference in log₂(FPKM + 1) with 95% CI",
           title="Gene-level regional differences")
    ax.axvline(0, color="#555555", ls="--", lw=1)
    ax.invert_yaxis()
    for group in ("rod", "cone"):
        ax.plot([], [], "o-", label=group.title(), color=COLORS[group])
    ax.legend(frameon=False)
    save(fig, folder, "03_gene_effects")

    # Mean-center complete, variable genes across samples before PCA.
    complete = log.dropna()
    variable = complete.loc[complete.var(axis=1) > 0]
    if len(variable) >= 2:
        centered = variable.T - variable.mean(axis=1)
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        components = pd.DataFrame(u[:, :2] * singular[:2], index=log.columns, columns=["PC1", "PC2"])
        explained = singular[:2] ** 2 / (singular ** 2).sum()
        fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
        for pair in pairs.itertuples():
            points = components.loc[[pair.peripheral, pair.macular]]
            ax.plot(points.PC1, points.PC2, color="#b7b7b7", lw=0.9)
            ax.annotate(pair.Index, points.loc[pair.macular], xytext=(4, 4), textcoords="offset points", fontsize=8)
        for region in REGIONS:
            points = components.loc[metadata.region == region]
            ax.scatter(points.PC1, points.PC2, label=region.title(), color=COLORS[region], s=36)
        ax.set(xlabel=f"PC1 ({explained[0]:.1%})", ylabel=f"PC2 ({explained[1]:.1%})",
               title="Exploratory PCA of retinal expression")
        ax.legend(frameon=False)
        save(fig, folder, "04_exploratory_pca")

    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    boxes = ax.boxplot([log[sample].dropna() for sample in samples], tick_labels=sample_labels,
                       showfliers=False, patch_artist=True)
    for box, region in zip(boxes["boxes"], metadata.region, strict=True):
        box.set_facecolor(COLORS[region])
        box.set_alpha(0.40)
    ax.set(ylabel="log₂(FPKM + 1)", title="Expression distributions by sample")
    ax.tick_params(axis="x", labelrotation=60)
    save(fig, folder, "05_sample_distributions")

    fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
    image = ax.imshow(log.corr(method="spearman").loc[samples, samples], vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set(xticks=range(len(samples)), xticklabels=sample_labels, yticks=range(len(samples)),
           yticklabels=sample_labels, title="Sample similarity across genes")
    ax.tick_params(axis="x", labelrotation=60)
    fig.colorbar(image, ax=ax, shrink=0.72, label="Spearman correlation")
    save(fig, folder, "06_sample_correlations")


def analyze(path, provenance, output):
    metadata_path = ROOT / "data/metadata.csv"
    genes_path = ROOT / "data/gene_sets.csv"
    metadata = pd.read_csv(metadata_path).set_index("sample_id")
    pairs = metadata.reset_index().pivot(index="eye_id", columns="region", values="sample_id")
    if pairs.shape != (8, 2) or set(pairs.columns) != set(REGIONS) or pairs.isna().any().any():
        raise ValueError("Metadata must have exactly one sample per region for each of eight eyes.")
    groups = pd.read_csv(genes_path).set_index("symbol").gene_set
    fpkm = load_matrix(path, metadata.index.tolist())
    log = np.log2(fpkm + 1)
    panel = log.reindex(groups.index)
    panel_qc = pd.DataFrame({"group": groups, "absent": ~groups.index.isin(log.index),
                             "missing_samples": panel.isna().sum(axis=1),
                             "all_zero": panel.eq(0).all(axis=1),
                             "constant": panel.nunique(axis=1, dropna=False).eq(1),
                             "max_fpkm_below_1": fpkm.reindex(groups.index).max(axis=1).lt(1)})
    subgroup = panel.dropna()
    if subgroup.empty:
        raise ValueError("No frozen markers have complete expression across the 16 samples.")
    z_input = subgroup.loc[subgroup.var(axis=1) > 0]
    z = z_input.sub(z_input.mean(axis=1), axis=0).div(z_input.std(axis=1, ddof=1), axis=0)
    signatures = z.groupby(groups).mean()
    if set(signatures.index) != {"rod", "cone"}:
        raise ValueError("Both signatures require at least one complete, variable frozen marker.")
    delta = contrast(subgroup, pairs)
    effects = summarize(delta)
    effects.insert(0, "group", groups)
    # Preserve one BH family across all frozen candidates, including missing and untestable genes.
    effects = effects.reindex(groups.index)
    effects["group"] = groups
    valid = effects.p_paired_t.notna()
    effects["q_bh"] = np.where(valid, stats.false_discovery_control(effects.p_paired_t.fillna(1)), np.nan)
    signature_delta = contrast(signatures, pairs)
    signature_effects = summarize(signature_delta)

    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()
    sample_qc = pd.DataFrame({"missing_genes": fpkm.isna().sum(), "zero_fpkm_genes": fpkm.eq(0).sum(),
                              "fpkm_below_1": fpkm.lt(1).sum(), "median_log2_fpkm_plus_1": log.median()})
    for name, frame in {
        "gene_effects": effects, "gene_differences": delta, "signature_scores": signatures,
        "signature_differences": signature_delta, "signature_effects": signature_effects,
        "marker_qc": panel_qc, "sample_qc": sample_qc, "sample_correlations": log.corr(method="spearman"),
    }.items():
        frame.to_csv(output / f"{name}.csv")
    plot_results(log, z, signatures, effects.loc[subgroup.index], groups, metadata, pairs, figures)
    manifest = {"input_fpkm_sha256": digest(path), "metadata_sha256": digest(metadata_path),
                "gene_sets_sha256": digest(genes_path), "input_provenance": provenance,
                "n_pairs": len(pairs), "signature_markers": {
                    group: z.index[groups.reindex(z.index).eq(group)].tolist() for group in ("rod", "cone")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    lines = ["# Paired expression summary", "", f"Input provenance: {provenance}", "",
             f"{len(pairs)} eye pairs and {len(groups)} frozen markers. "
             f"{len(subgroup)} markers had complete expression. "
             f"{valid.sum()} had estimable paired t-tests, and {(effects.q_bh < 0.05).sum()} had BH q < 0.05.", ""]
    for group, row in signature_effects.iterrows():
        lines.append(f"- {group.title()} signature mean macular-minus-peripheral difference "
                     f"{row.mean_delta:+.3f} z-score units (95% CI {row.ci95_low:+.3f} to {row.ci95_high:+.3f}). "
                     f"Macular higher in {int(row.macular_higher)}/{len(pairs)} pairs.")
    (output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fpkm", type=Path, required=True)
    parser.add_argument("--provenance", required=True, help="Document source, GSM mapping, and FPKM units.")
    parser.add_argument("--out", type=Path, default=ROOT / "analysis-output", help="New output directory.")
    args = parser.parse_args()
    try:
        analyze(args.fpkm, args.provenance, args.out)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Analysis stopped: {error}\n")
