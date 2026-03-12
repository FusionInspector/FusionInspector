#!/usr/bin/env python3
"""
plot_reciprocal_fusions.py

Visualize reciprocal gene fusion pairs from FusionInspector output.
Shows stacked gene models (one per gene, same bp scale) with:
  - Breakpoint lines from both fusion directions
  - Cross-panel connector lines linking paired breakpoints per fusion event
  - Chromosome identifiers and parsed annotation summary

Uses the GTF's local intron-compressed coordinate system for display.

Usage:
    python plot_reciprocal_fusions.py \\
        --gtf  finspector.gtf \\
        --tsv  finspector.FusionInspector.fusions.abridged.tsv \\
        --out  reciprocal_plots/
"""

import re
import sys
import os
import argparse
import collections
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch

# ── visual constants (halved from original to keep image compact) ─────────────
UTR_H        = 0.13   # half-height of UTR exon rectangles
CDS_H        = 0.22   # half-height of CDS overlay rectangles
ROW_SPACING  = 0.60   # vertical spacing between transcript rows
FIRST_ROW_Y  = 0.75   # y of lowest transcript row
AXIS_Y       = 0.0    # y of the 5'→3' direction arrow

GENE_A_COL   = '#3A6FC4'
GENE_A_CDS   = '#1E3F7A'
GENE_B_COL   = '#E07B35'
GENE_B_CDS   = '#8B4A18'

BP_AB_COL    = '#D62728'   # red   – GeneA→GeneB breakpoints
BP_BA_COL    = '#2CA02C'   # green – GeneB→GeneA breakpoints

CONN_ALPHA   = 0.35
CONN_LW      = 2.0
INTRON_LW    = 0.7


# ── GTF parsing ───────────────────────────────────────────────────────────────

def _attr(s, key):
    m = re.search(r'\b' + re.escape(key) + r'\s+"([^"]*)"', s)
    return m.group(1) if m else None


def parse_gtf(gtf_file):
    feats = []
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            contig, _, feat, start, end, _, strand, _, attrs = cols
            if feat not in ('exon', 'CDS'):
                continue
            feats.append({
                'contig':     contig,
                'feature':    feat,
                'start':      int(start),
                'end':        int(end),
                'strand':     strand,
                'tid':        _attr(attrs, 'transcript_id'),
                'gene_name':  _attr(attrs, 'gene_name'),
                'tx_type':    _attr(attrs, 'transcript_type'),
                'orig_coord': _attr(attrs, 'orig_coord_info'),
            })
    return feats


# ── fusion TSV parsing ────────────────────────────────────────────────────────

def parse_fusion_tsv(tsv_file):
    fusions = []
    with open(tsv_file) as fh:
        header = fh.readline().lstrip('#').rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            row['LeftLocalBreakpoint']  = int(row['LeftLocalBreakpoint'])
            row['RightLocalBreakpoint'] = int(row['RightLocalBreakpoint'])
            row['JunctionReadCount']    = int(row['JunctionReadCount'])
            row['SpanningFragCount']    = int(row['SpanningFragCount'])
            fusions.append(row)
    return fusions


# ── annotation parsing ────────────────────────────────────────────────────────

# Known fusion databases to look for in the annots field (in display order)
_KNOWN_DBS = ['ChimerKB', 'ChimerPub', 'ChimerSeq', 'DEEPEST2019',
              'TumorFusionsNAR2018', 'Cosmic', 'TCGA_StarF2019']

# Short display aliases
_DB_ALIAS = {'TumorFusionsNAR2018': 'NAR2018', 'TCGA_StarF2019': 'TCGA'}


def parse_annots(annots_str):
    """
    Return a dict with keys:
        chromosomal  – 'Intrachromosomal [chrX:Y.YYMb]' or 'Interchromosomal [chrX--chrY]'
        databases    – list of recognised DB names
        oncogenes    – list of genes flagged as Oncogene
        flags        – list of structural flags (NEIGHBORS_OVERLAP, LOCAL_REARRANGEMENT, …)
    """
    info = {'chromosomal': None, 'databases': [], 'oncogenes': [], 'flags': []}

    m = re.search(r'INTRACHROMOSOMAL\[([^\]]+)\]', annots_str)
    if m:
        info['chromosomal'] = f'Intrachromosomal  [{m.group(1)}]'
    else:
        m = re.search(r'INTERCHROMOSOMAL\[([^\]]+)\]', annots_str)
        if m:
            info['chromosomal'] = f'Interchromosomal  [{m.group(1)}]'

    dbs = [db for db in _KNOWN_DBS if db in annots_str]
    info['databases'] = list(dict.fromkeys(dbs))

    info['oncogenes'] = list(dict.fromkeys(
        re.findall(r'"(\w+):Oncogene"', annots_str)))

    for flag in ('NEIGHBORS_OVERLAP', 'LOCAL_REARRANGEMENT', 'ONLY_REF_SPLICE',
                 'INCL_NON_REF_SPLICE'):
        if flag in annots_str:
            info['flags'].append(flag)

    return info


def format_annots_line(all_fusions):
    """Build a compact annotation summary string from a list of fusion rows."""
    # Merge annots across all fusions in the pair
    chromosomal, databases, oncogenes, flags = None, [], [], []
    for f in all_fusions:
        a = parse_annots(f.get('annots', ''))
        if chromosomal is None and a['chromosomal']:
            chromosomal = a['chromosomal']
        databases  += [d for d in a['databases']  if d not in databases]
        oncogenes  += [g for g in a['oncogenes']  if g not in oncogenes]
        flags      += [fl for fl in a['flags']    if fl not in flags]

    parts = []
    if chromosomal:
        parts.append(chromosomal)
    if oncogenes:
        parts.append('Oncogenes: ' + ', '.join(oncogenes))
    if databases:
        parts.append('Databases: ' + ', '.join(_DB_ALIAS.get(d, d) for d in databases))
    if flags:
        parts.append(' | '.join(flags))
    return '   ·   '.join(parts) if parts else '(no annotation)'


def chrom_from_breakpoint(bp_str):
    """Extract chromosome from 'chr12:57093598:+' → 'chr12'."""
    return bp_str.split(':')[0] if bp_str else '?'


# ── reciprocal pair detection ─────────────────────────────────────────────────

def find_reciprocal_pairs(fusions):
    by_pair = collections.defaultdict(lambda: {'ab': [], 'ba': []})
    for f in fusions:
        lg = f['LeftGene'].split('^')[0]
        rg = f['RightGene'].split('^')[0]
        key = tuple(sorted([lg, rg]))
        if lg == key[0]:
            by_pair[key]['ab'].append(f)
        else:
            by_pair[key]['ba'].append(f)

    result = []
    for (gene_a, gene_b), d in by_pair.items():
        if d['ab'] and d['ba']:
            result.append({'gene_a': gene_a, 'gene_b': gene_b,
                           'ab_fusions': d['ab'], 'ba_fusions': d['ba']})
    return result


# ── feature / coordinate helpers ──────────────────────────────────────────────

def gene_features_from_contig(all_feats, contig, gene_name):
    return [f for f in all_feats
            if f['contig'] == contig and f['gene_name'] == gene_name]


def _unique_exon_maps(gene_feats):
    seen = {}
    for f in gene_feats:
        if f['feature'] != 'exon' or not f['orig_coord']:
            continue
        parts = f['orig_coord'].split(',')
        if len(parts) < 4:
            continue
        g_start, g_end = int(parts[1]), int(parts[2])
        key = (g_start, g_end)
        if key not in seen:
            seen[key] = (g_start, g_end, parts[3], f['start'], f['end'])
    return sorted(seen.values(), key=lambda x: x[3])


def genomic_to_local(genomic_pos, gene_feats):
    emaps = _unique_exon_maps(gene_feats)
    if not emaps:
        return None
    strand = emaps[0][2]

    for g_start, g_end, g_strand, l_start, l_end in emaps:
        if g_start <= genomic_pos <= g_end:
            if g_strand == '+':
                return l_start + (genomic_pos - g_start)
            else:
                return l_start + (g_end - genomic_pos)

    if strand == '+':
        by_g = sorted(emaps, key=lambda x: x[0])
        for e1, e2 in zip(by_g, by_g[1:]):
            if e1[1] < genomic_pos < e2[0]:
                frac = (genomic_pos - e1[1]) / (e2[0] - e1[1])
                return e1[4] + frac * (e2[3] - e1[4])
    else:
        by_g = sorted(emaps, key=lambda x: x[0], reverse=True)
        for e_hi, e_lo in zip(by_g, by_g[1:]):
            if e_lo[1] < genomic_pos < e_hi[0]:
                frac = (e_hi[0] - genomic_pos) / (e_hi[0] - e_lo[1])
                return e_hi[4] + frac * (e_lo[3] - e_hi[4])
    return None


# ── transcript selection ──────────────────────────────────────────────────────

TX_PRIORITY = {'protein_coding': 0, 'processed_transcript': 1,
               'retained_intron': 2, 'nonsense_mediated_decay': 3}


def build_transcripts(gene_feats):
    tx = collections.defaultdict(lambda: {'exons': [], 'cds': [], 'type': None})
    for f in gene_feats:
        tid = f['tid']
        if f['feature'] == 'exon':
            tx[tid]['exons'].append((f['start'], f['end']))
            if tx[tid]['type'] is None:
                tx[tid]['type'] = f['tx_type']
        elif f['feature'] == 'CDS':
            tx[tid]['cds'].append((f['start'], f['end']))
    for t in tx.values():
        t['exons'].sort()
        t['cds'].sort()
    return dict(tx)


def select_transcripts(transcripts, max_n=7):
    def sort_key(tid):
        t = transcripts[tid]
        return (TX_PRIORITY.get(t['type'] or '', 10), -len(t['exons']))

    seen, selected = set(), []
    for tid in sorted(transcripts, key=sort_key):
        struct = tuple(transcripts[tid]['exons'])
        if struct not in seen and transcripts[tid]['exons']:
            seen.add(struct)
            selected.append(tid)
        if len(selected) >= max_n:
            break
    return selected


# ── drawing helpers ───────────────────────────────────────────────────────────

def draw_transcript_row(ax, row_y, exons, cds_list, offset, exon_col, cds_col):
    if not exons:
        return
    ax.plot([exons[0][0] - offset, exons[-1][1] - offset],
            [row_y, row_y], color=exon_col, lw=INTRON_LW, zorder=1)
    for s, e in exons:
        ax.add_patch(mpatches.Rectangle(
            (s - offset, row_y - UTR_H), e - s, 2 * UTR_H,
            lw=0.3, ec='black', fc=exon_col, zorder=2))
    for s, e in cds_list:
        ax.add_patch(mpatches.Rectangle(
            (s - offset, row_y - CDS_H), e - s, 2 * CDS_H,
            lw=0.4, ec='black', fc=cds_col, zorder=3))


def _stagger_labels(breakpoints, offset, max_span, n_levels=3, prox_frac=0.04):
    """Assign stagger level to each breakpoint to avoid overlapping labels."""
    threshold = max_span * prox_frac
    levels = [0] * len(breakpoints)
    indexed = sorted(enumerate(breakpoints), key=lambda x: x[1][0])
    prev_x, cur_level = -1e18, 0
    for rank, (orig_i, (pos, col, lbl)) in enumerate(indexed):
        x = pos - offset
        cur_level = (cur_level + 1) % n_levels if (rank > 0 and abs(x - prev_x) < threshold) else 0
        levels[orig_i] = cur_level
        prev_x = x
    return [(pos, col, lbl, levels[i]) for i, (pos, col, lbl) in enumerate(breakpoints)]


def draw_gene_panel(ax, gene_name, chrom, sel_tids, transcripts,
                    breakpoints, offset, max_span, exon_col, cds_col):
    """Draw transcripts + breakpoints. Returns (y_top, y_bottom) in data coords."""
    n_rows = len(sel_tids)
    y_top = FIRST_ROW_Y + n_rows * ROW_SPACING + 0.30

    for i, tid in enumerate(sel_tids):
        row_y = FIRST_ROW_Y + i * ROW_SPACING
        t = transcripts[tid]
        draw_transcript_row(ax, row_y, t['exons'], t['cds'],
                            offset, exon_col, cds_col)

    # 5'→3' arrow
    ax.annotate('', xy=(max_span * 1.01, AXIS_Y), xytext=(0, AXIS_Y),
                arrowprops=dict(arrowstyle='->', color='#444', lw=1.2))

    # Staggered breakpoint lines + labels
    label_ys = [y_top * 0.96, y_top * 0.80, y_top * 0.64]
    for local_pos, color, label, lvl in _stagger_labels(breakpoints, offset, max_span):
        x = local_pos - offset
        ax.axvline(x=x, color=color, ls='--', lw=1.5, alpha=0.85, zorder=10)
        ax.text(x + max_span * 0.005, label_ys[min(lvl, 2)], label,
                color=color, fontsize=7, va='top', ha='left', zorder=11,
                bbox=dict(fc='white', ec='none', alpha=0.65, pad=0.8))

    margin = max_span * 0.02
    ax.set_xlim(-margin, max_span + margin)
    ax.set_ylim(AXIS_Y - 0.35, y_top)
    ax.set_title(f'{gene_name}   ({chrom},  5′ → 3′)', fontsize=12,
                 fontweight='bold', loc='left', pad=3)
    ax.set_xlabel('Local contig position  [bp, introns compressed]', fontsize=9)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.yaxis.set_visible(False)
    ax.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f'{int(x + offset):,}'))

    return y_top, AXIS_Y


def draw_connections(fig, ax_a, ax_b, connections, a_lmin, b_lmin,
                     y_bottom_a, y_top_b):
    """Draw one ConnectionPatch per fusion event between the two gene panels."""
    for a_local, b_local, color, label in connections:
        con = ConnectionPatch(
            xyA=(a_local - a_lmin, y_bottom_a), coordsA='data', axesA=ax_a,
            xyB=(b_local - b_lmin, y_top_b),   coordsB='data', axesB=ax_b,
            color=color, lw=CONN_LW, ls='-', alpha=CONN_ALPHA,
            clip_on=False, zorder=20)
        fig.add_artist(con)


# ── main plotting function ────────────────────────────────────────────────────

def plot_reciprocal_pair(pair, all_feats, out_dir, max_tx=7):
    gene_a, gene_b = pair['gene_a'], pair['gene_b']
    ab_fusions = pair['ab_fusions']
    ba_fusions = pair['ba_fusions']

    contig_ab = f'{gene_a}--{gene_b}'
    contig_ba = f'{gene_b}--{gene_a}'

    feats_a = gene_features_from_contig(all_feats, contig_ab, gene_a)
    feats_b = gene_features_from_contig(all_feats, contig_ba, gene_b)

    for gene, feats, contig in [(gene_a, feats_a, contig_ab),
                                 (gene_b, feats_b, contig_ba)]:
        if not feats:
            print(f'  WARNING: no features for {gene} in {contig}', file=sys.stderr)
            return

    # Chromosomes
    chrom_a = chrom_from_breakpoint(ab_fusions[0]['LeftBreakpoint'])
    chrom_b = chrom_from_breakpoint(ab_fusions[0]['RightBreakpoint'])

    a_lmin = min(f['start'] for f in feats_a)
    a_lmax = max(f['end']   for f in feats_a)
    b_lmin = min(f['start'] for f in feats_b)
    b_lmax = max(f['end']   for f in feats_b)
    max_span = max(a_lmax - a_lmin, b_lmax - b_lmin)

    tx_a, tx_b = build_transcripts(feats_a), build_transcripts(feats_b)
    sel_a = select_transcripts(tx_a, max_tx)
    sel_b = select_transcripts(tx_b, max_tx)

    # ── breakpoints & connections ─────────────────────────────────────────────
    bp_a, bp_b, conns = [], [], []

    for f in ab_fusions:
        j, s = f['JunctionReadCount'], f['SpanningFragCount']
        lbl = f'{gene_a}→{gene_b}\nJ={j} S={s}'
        a_pos = f['LeftLocalBreakpoint']
        g_pos_b = int(f['RightBreakpoint'].split(':')[1])
        b_pos = genomic_to_local(g_pos_b, feats_b)
        bp_a.append((a_pos, BP_AB_COL, lbl))
        if b_pos is not None:
            bp_b.append((b_pos, BP_AB_COL, lbl))
            conns.append((a_pos, b_pos, BP_AB_COL, lbl))
        else:
            print(f'  WARNING: could not map {gene_b} breakpoint {g_pos_b}',
                  file=sys.stderr)

    for f in ba_fusions:
        j, s = f['JunctionReadCount'], f['SpanningFragCount']
        lbl = f'{gene_b}→{gene_a}\nJ={j} S={s}'
        b_pos = f['LeftLocalBreakpoint']
        g_pos_a = int(f['RightBreakpoint'].split(':')[1])
        a_pos = genomic_to_local(g_pos_a, feats_a)
        bp_b.append((b_pos, BP_BA_COL, lbl))
        if a_pos is not None:
            bp_a.append((a_pos, BP_BA_COL, lbl))
            conns.append((a_pos, b_pos, BP_BA_COL, lbl))
        else:
            print(f'  WARNING: could not map {gene_a} breakpoint {g_pos_a}',
                  file=sys.stderr)

    def dedup(lst):
        seen, out = set(), []
        for pos, col, lbl in lst:
            k = (round(pos), col)
            if k not in seen:
                seen.add(k); out.append((pos, col, lbl))
        return out

    bp_a = dedup(bp_a)
    bp_b = dedup(bp_b)

    # ── annotation summary ────────────────────────────────────────────────────
    annots_line = format_annots_line(ab_fusions + ba_fusions)

    # ── figure ────────────────────────────────────────────────────────────────
    h_a     = FIRST_ROW_Y + len(sel_a) * ROW_SPACING + 0.60
    h_b     = FIRST_ROW_Y + len(sel_b) * ROW_SPACING + 0.60
    h_title = 1.6   # inches reserved for the title / annotation block

    fig = plt.figure(figsize=(14, h_title + h_a + h_b + 0.5))
    gs  = fig.add_gridspec(3, 1,
                           height_ratios=[h_title, h_a, h_b],
                           hspace=0.55)
    ax_t = fig.add_subplot(gs[0])   # hidden title panel
    ax_a = fig.add_subplot(gs[1])
    ax_b = fig.add_subplot(gs[2])
    ax_t.set_axis_off()

    y_top_a, y_bot_a = draw_gene_panel(
        ax_a, gene_a, chrom_a, sel_a, tx_a, bp_a,
        a_lmin, max_span, GENE_A_COL, GENE_A_CDS)
    y_top_b, y_bot_b = draw_gene_panel(
        ax_b, gene_b, chrom_b, sel_b, tx_b, bp_b,
        b_lmin, max_span, GENE_B_COL, GENE_B_CDS)

    draw_connections(fig, ax_a, ax_b, conns, a_lmin, b_lmin, y_bot_a, y_top_b)

    # ── legend ────────────────────────────────────────────────────────────────
    legend_elems = [
        mpatches.Patch(fc=GENE_A_COL, ec='black', lw=0.3, label=f'{gene_a}  exon/UTR'),
        mpatches.Patch(fc=GENE_A_CDS, ec='black', lw=0.3, label=f'{gene_a}  CDS'),
        mpatches.Patch(fc=GENE_B_COL, ec='black', lw=0.3, label=f'{gene_b}  exon/UTR'),
        mpatches.Patch(fc=GENE_B_CDS, ec='black', lw=0.3, label=f'{gene_b}  CDS'),
        plt.Line2D([0], [0], color=BP_AB_COL, ls='--', lw=1.5,
                   label=f'{gene_a}→{gene_b}  breakpoints'),
        plt.Line2D([0], [0], color=BP_BA_COL, ls='--', lw=1.5,
                   label=f'{gene_b}→{gene_a}  breakpoints'),
        plt.Line2D([0], [0], color='grey', ls='-', lw=CONN_LW, alpha=CONN_ALPHA,
                   label='fusion connector (per event)'),
    ]
    fig.legend(handles=legend_elems, loc='upper right', fontsize=9,
               frameon=True, bbox_to_anchor=(0.99, 0.99))

    # ── titles ────────────────────────────────────────────────────────────────
    ab_jr = sum(f['JunctionReadCount'] for f in ab_fusions)
    ba_jr = sum(f['JunctionReadCount'] for f in ba_fusions)
    n_ev  = len(conns)

    title_text = (
        f'Reciprocal fusion:  {gene_a}--{gene_b}  ↔  {gene_b}--{gene_a}\n'
        f'{gene_a}→{gene_b}: {ab_jr} junc reads   |   '
        f'{gene_b}→{gene_a}: {ba_jr} junc reads   '
        f'({n_ev} fusion event{"s" if n_ev != 1 else ""} connected)')

    # Place bold title in the upper portion of the hidden title panel
    ax_t.text(0.5, 0.90, title_text,
              ha='center', va='top', fontsize=12, fontweight='bold',
              linespacing=1.6, transform=ax_t.transAxes)

    # Place annotation in the lower portion of the title panel
    wrapped = textwrap.fill(annots_line, width=110)
    ax_t.text(0.5, 0.15, wrapped,
              ha='center', va='bottom', fontsize=9, color='#444444',
              style='italic', transform=ax_t.transAxes)

    out_path = os.path.join(out_dir, f'{gene_a}__{gene_b}.reciprocal_fusion.pdf')
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  → {out_path}')


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True, help='finspector.gtf')
    ap.add_argument('--tsv', required=True,
                    help='finspector.FusionInspector.fusions.abridged.tsv')
    ap.add_argument('--out', default='.', help='Output directory')
    ap.add_argument('--max_transcripts', type=int, default=7,
                    help='Max unique isoforms per gene (default 7)')
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print('Parsing GTF …', file=sys.stderr)
    all_feats = parse_gtf(args.gtf)

    print('Parsing fusion TSV …', file=sys.stderr)
    fusions = parse_fusion_tsv(args.tsv)

    print('Finding reciprocal pairs …', file=sys.stderr)
    pairs = find_reciprocal_pairs(fusions)

    if not pairs:
        print('No reciprocal pairs found.', file=sys.stderr)
        return

    for pair in pairs:
        print(f"Plotting {pair['gene_a']}--{pair['gene_b']} "
              f"↔ {pair['gene_b']}--{pair['gene_a']} …", file=sys.stderr)
        plot_reciprocal_pair(pair, all_feats, args.out,
                             max_tx=args.max_transcripts)

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
