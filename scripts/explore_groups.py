"""Descriptive comparison of published tissue-group FPKM estimates; no inference."""
import csv
import hashlib
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def main():
    source = ROOT / 'data/retinal_group_expression.csv'
    with source.open() as f:
        expression = list(csv.DictReader(f))
    with (ROOT / 'data/gene_sets.csv').open() as f:
        panel = list(csv.DictReader(f))
    results = []
    for marker in panel:
        symbol = marker['symbol']
        hits = [r for r in expression if r['gene_short_name'] == symbol]
        entry = dict(marker, match_count=len(hits), match_status='unique' if len(hits)==1 else 'missing' if not hits else 'ambiguous')
        if len(hits) == 1:
            row = hits[0]
            m, p = [float(row[f'{region}_Retina_FPKM']) for region in ('Macular', 'Peripheral')]
            if not all(math.isfinite(x) and x >= 0 for x in (m, p)):
                raise ValueError(f'Invalid expression: {symbol}')
            entry.update(tracking_id=row['tracking_id'], macular_fpkm=m, peripheral_fpkm=p,
                         macular_status=row['Macular_Retina_status'], peripheral_status=row['Peripheral_Retina_status'],
                         macular_below_1=m<1, peripheral_below_1=p<1,
                         log2_macular_plus1=math.log2(m+1), log2_peripheral_plus1=math.log2(p+1),
                         macular_over_peripheral=m/p if p>0 else '',
                         log2_macular_over_peripheral=math.log2(m/p) if m>0 and p>0 else '',
                         log2_plus1_difference=math.log2(m+1)-math.log2(p+1))
        results.append(entry)
    out = ROOT / 'results/exploratory'
    figures = ROOT / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(exist_ok=True)
    keys = list(dict.fromkeys(k for r in results for k in r))
    with (out / 'marker_expression.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader(); writer.writerows(results)
    valid = [r for r in results if r['match_status']=='unique' and r['macular_status']==r['peripheral_status']=='OK']
    ratio_rows = [r for r in valid if r['log2_macular_over_peripheral']!='']
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'svg.fonttype':'none'})
    labels = [f"{r['symbol']} ({r['gene_set']})" for r in valid]
    fig, ax = plt.subplots(figsize=(7,9))
    values = np.array([[r['log2_macular_plus1'],r['log2_peripheral_plus1']] for r in valid])
    image = ax.imshow(values, cmap='viridis', aspect='auto')
    ax.set(xticks=[0,1], xticklabels=['Macular retina','Peripheral retina'], yticks=range(len(valid)), yticklabels=labels,
           title='Rod/cone panel: tissue-group expression', xlabel='Published tissue group', ylabel='Predefined gene (association)')
    fig.colorbar(image, ax=ax, label='log2(FPKM + 1)')
    fig.text(.5,.012,'Source: GSE94437 processed FPKM | Group estimates; no donor-level inference',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,1])
    for ext in ('png','svg'): fig.savefig(figures / f'group_expression_heatmap.{ext}', dpi=200)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8,9))
    y = np.arange(len(ratio_rows))
    x = [r['log2_macular_over_peripheral'] for r in ratio_rows]
    colors = ['#276678' if r['gene_set']=='rod' else '#ad653c' for r in ratio_rows]
    ax.barh(y,x,color=colors,height=.65)
    ax.axvline(0,color='0.3',linewidth=.8)
    ax.set(yticks=y,yticklabels=[f"{r['symbol']} ({r['gene_set']})" for r in ratio_rows],
           xlabel='log2(Macular FPKM / Peripheral FPKM)', ylabel='Predefined gene (association)',
           title='Regional ratios of tissue-group estimates')
    ax.invert_yaxis()
    ax.set_title('Regional ratios of tissue-group estimates', pad=32)
    ax.text(.01,1.005,'Negative: higher peripheral',transform=ax.transAxes,fontsize=9)
    ax.text(.99,1.005,'Positive: higher macular',transform=ax.transAxes,ha='right',fontsize=9)
    fig.text(.5,.012,'Source: GSE94437 processed FPKM | No p-values or donor variability estimates',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,1])
    for ext in ('png','svg'): fig.savefig(figures / f'group_regional_ratios.{ext}',dpi=200)
    plt.close(fig)
    summary = {'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
               'predefined_markers':len(panel),'uniquely_matched':len(valid),
               'missing':[r['symbol'] for r in results if r['match_status']=='missing'],
               'ambiguous':[r['symbol'] for r in results if r['match_status']=='ambiguous'],
               'group_summary':{group:{'genes':sum(r['gene_set']==group for r in ratio_rows),
                                      'higher_macular':sum(r['gene_set']==group and r['macular_fpkm']>r['peripheral_fpkm'] for r in ratio_rows),
                                      'higher_peripheral':sum(r['gene_set']==group and r['macular_fpkm']<r['peripheral_fpkm'] for r in ratio_rows)} for group in ('rod','cone')},
               'inference':'Descriptive only. Genes are not donor replicates; no hypothesis tests performed.'}
    (out / 'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines = ['# Exploratory results', '', 'Descriptive comparison of published tissue-group estimates; no donor-level inference.', '',
             f"Matched {len(valid)} of {len(panel)} predefined markers uniquely with OK status in both groups.", '']
    for group, stats in summary['group_summary'].items():
        lines.append(f"- {group}: {stats['higher_macular']} higher macular; {stats['higher_peripheral']} higher peripheral, among {stats['genes']} genes.")
    lines += ['', '## Examples', '']
    for symbol in ('RHO','GNAT1','OPN1LW','OPN1MW','GNAT2'):
        r = next(r for r in valid if r['symbol']==symbol)
        lines.append(f"- {symbol}: macular {r['macular_fpkm']:.5g} FPKM; peripheral {r['peripheral_fpkm']:.5g} FPKM; macular/peripheral ratio {r['macular_over_peripheral']:.3f}.")
    lines += ['', 'The panel pattern differs from the initial directional hypothesis. These observations describe the supplied group estimates only; they do not establish donor consistency, statistical significance, cell counts or causality.', '',
              'Inspect figures/group_regional_ratios.png and figures/group_expression_heatmap.png. Full values and QC are in marker_expression.csv. Original bounds remain in data/retinal_group_expression.csv; they are not confidence intervals for this contrast.', '']
    (out / 'summary.md').write_text('\n'.join(lines))
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
