#!/usr/bin/env python
"""
Standalone script to build phylogenetic pie plots from a tree file and a cluster
abundances file. Use this when you want to visualize branch proportions with a
custom or edited tree (e.g. one you don't agree with from BuildTree).

Inputs:
  1) Tree file: either
     - BuildTree posteriors TSV (e.g. indiv_build_tree_posteriors.tsv): columns
       n_iter, likelihood, edges. Uses the first data row (best tree). Edges
       column format: "None-1,1-2,2-4,..."
     - Simple edge list TSV: header "Parent_ID\tChild_ID", one edge per line.
       Root must be cluster 1. Example:
         Parent_ID	Child_ID
         1	2
         1	3
         2	4

  2) Abundances file: either
     - Constrained CCF style (e.g. indiv_constrained_ccf.tsv): columns
       Patient_ID, Sample_ID, Cell_population, Constrained_CCF. Cluster is
       inferred from Cell_population (e.g. CL1_CL2 -> cluster 2).
     - Simple format: columns Sample_ID, Cluster_ID, Abundance (or
       Constrained_CCF). One row per (sample, cluster).

Output:
  - One pie plot SVG per sample in the requested output directory (and
    optionally base64 PNGs if used as a library).
"""

from __future__ import print_function

import argparse
import os
import re
import sys

# Allow running from repo root or from scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from io import BytesIO
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import BuildTree.Tree
from BuildTree.Node import Node
from output.PhylogicOutput import PhylogicOutput, ClusterColors


def parse_tree_tsv_buildtree(tree_path):
    """Parse BuildTree posteriors TSV. Returns list of (parent_id, child_id); root edge (None,1) excluded."""
    with open(tree_path, 'r') as f:
        header = f.readline().strip().split('\t')
        if 'edges' not in header:
            raise ValueError('Tree file must have an "edges" column (BuildTree format).')
        idx = header.index('edges')
        first_line = f.readline()
        if not first_line.strip():
            raise ValueError('Tree file has no data rows.')
        edges_str = first_line.strip().split('\t')[idx]
    edges = PhylogicOutput.reformat_edges_for_input(edges_str)
    # Exclude (None, 1) so we don't need node None when building tree
    return [(p, c) for (p, c) in edges if p is not None]


def parse_tree_simple(tree_path):
    """Parse simple edge list: Parent_ID, Child_ID. Returns list of (parent_id, child_id)."""
    edges = []
    with open(tree_path, 'r') as f:
        header = f.readline().strip().lower().split('\t')
        # Accept various column names
        pcol = 'parent_id' if 'parent_id' in header else (header[0] if header else 'parent_id')
        ccol = 'child_id' if 'child_id' in header else (header[1] if len(header) > 1 else 'child_id')
        try:
            pi, ci = header.index(pcol), header.index(ccol)
        except (ValueError, IndexError):
            pi, ci = 0, 1
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            try:
                p = int(parts[pi].strip())
                c = int(parts[ci].strip())
            except ValueError:
                continue
            edges.append((p, c))
    return edges


def load_tree(tree_path, tree_format='auto'):
    """
    Load tree from file. tree_format one of 'auto', 'buildtree', 'simple'.
    Returns BuildTree.Tree.Tree instance with root set to cluster 1.
    """
    if tree_format == 'auto':
        with open(tree_path, 'r') as f:
            header = f.readline().strip().split('\t')
        if 'edges' in header and 'n_iter' in header:
            edges = parse_tree_tsv_buildtree(tree_path)
        else:
            edges = parse_tree_simple(tree_path)
    elif tree_format == 'buildtree':
        edges = parse_tree_tsv_buildtree(tree_path)
    else:
        edges = parse_tree_simple(tree_path)

    if not edges:
        raise ValueError('No edges found in tree file.')

    node_ids = set()
    for (p, c) in edges:
        node_ids.add(p)
        node_ids.add(c)
    if 1 not in node_ids:
        raise ValueError('Tree must include cluster 1 (clonal/root).')

    tree = BuildTree.Tree.Tree()
    for nid in node_ids:
        tree.add_node(nid, data=None)
    tree.add_edges(edges)
    tree.set_root(tree.get_node_by_id(1))
    return tree


def parse_abundances_constrained_ccf(abundances_path):
    """Parse constrained CCF style file. Returns {sample_id: {cluster_id: value}}."""
    out = {}
    with open(abundances_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = dict(zip(header, line.strip().split('\t')))
            sample = fields['Sample_ID']
            cell_pop = fields['Cell_population']
            # Last cluster in path (e.g. CL1_CL2 -> 2)
            m = re.search(r'(\d{1,2})$', cell_pop)
            if not m:
                continue
            c = int(m.group(1))
            val = float(fields['Constrained_CCF'])
            out.setdefault(sample, {})[c] = val
    return out


def parse_abundances_simple(abundances_path):
    """Parse simple Sample_ID, Cluster_ID, Abundance file."""
    out = {}
    with open(abundances_path, 'r') as f:
        header = f.readline().strip().split('\t')
        header_lower = [h.strip().lower() for h in header]
        # Column indices
        sidx = header_lower.index('sample_id') if 'sample_id' in header_lower else 0
        cidx = header_lower.index('cluster_id') if 'cluster_id' in header_lower else 1
        if 'abundance' in header_lower:
            vidx = header_lower.index('abundance')
        elif 'constrained_ccf' in header_lower:
            vidx = header_lower.index('constrained_ccf')
        else:
            vidx = 2 if len(header_lower) > 2 else 1
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < max(sidx, cidx, vidx) + 1:
                continue
            sample = parts[sidx].strip()
            try:
                c = int(parts[cidx].strip())
                val = float(parts[vidx].strip())
            except (ValueError, IndexError):
                continue
            out.setdefault(sample, {})[c] = val
    return out


def load_abundances(abundances_path, format='auto'):
    """
    Load abundances. format one of 'auto', 'constrained_ccf', 'simple'.
    Returns {sample_id: {cluster_id: abundance}}.
    """
    if format == 'auto':
        with open(abundances_path, 'r') as f:
            header = f.readline().strip().lower().split('\t')
        if 'cell_population' in header and 'constrained_ccf' in header:
            return parse_abundances_constrained_ccf(abundances_path)
        return parse_abundances_simple(abundances_path)
    elif format == 'constrained_ccf':
        return parse_abundances_constrained_ccf(abundances_path)
    else:
        return parse_abundances_simple(abundances_path)


def safe_cluster_color(cluster_id):
    """Hex color for cluster; safe for large cluster IDs (wraps with modulo)."""
    n = len(ClusterColors.color_list)
    return ClusterColors.get_hex_string(cluster_id % n)


def make_pie_plot_standalone(tree, cluster_abundances, outdir='', sample='', dpi=150,
                             negative_tolerance=1e-6):
    """
    Build pie plot for one sample (same logic as PhylogicOutput.make_pie_plot but
    uses safe color for any cluster id).
    Returns base64-encoded PNG string.

    Raises ValueError if a cluster in the tree has no abundance entry for this sample, or if the
    abundances imply a negative pie slice for some cluster -- the latter usually means a
    hand-edited tree topology is inconsistent with the original CCF data for this particular
    sample (e.g. a child reparented such that its cumulative abundance now exceeds its new
    parent's, violating the pigeonhole principle the tree was originally built under). Both are
    reported with enough detail to identify the offending cluster/sample rather than surfacing a
    raw matplotlib error.
    """
    plt.figure(figsize=(1, 1))
    try:
        ax = plt.gca()
        pie_slices = {}
        node_order = list(tree.traverse_by_branch())
        for level, clusters in enumerate(tree.traverse_by_level()):
            for node in clusters:
                nid = node.identifier
                if nid not in cluster_abundances:
                    raise ValueError('Cluster {} missing in abundances for sample "{}".'.format(nid, sample))
                pie_slices[node] = cluster_abundances[nid]
                if node.parent:
                    pie_slices[node.parent] -= cluster_abundances[nid]

            negative = [(node.identifier, pie_slices[node]) for node in node_order
                       if node in pie_slices and pie_slices[node] < -negative_tolerance]
            if negative:
                detail = ', '.join('cluster {} = {:.4g}'.format(nid, val) for nid, val in negative)
                raise ValueError(
                    'Negative pie slice(s) for sample "{}": {}. This usually means the tree '
                    "topology is inconsistent with this sample's abundance data -- e.g. a "
                    "hand-edited parent/child relationship where a child's abundance exceeds its "
                    'new parent\'s, violating the pigeonhole principle the original tree was built '
                    'under.'.format(sample, detail))

            x = []
            colors = []
            for node in node_order:
                if node in pie_slices:
                    # Clip any within-tolerance negative floating-point noise to exactly 0, since
                    # genuinely negative values (beyond tolerance) were already raised above.
                    x.append(max(pie_slices[node], 0.))
                    colors.append(safe_cluster_color(node.identifier))
            ax.pie(x, colors=colors, radius=.9 - (.1 * level))
        ax.set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        bytes_io = BytesIO()
        plt.savefig(bytes_io, bbox_inches='tight', pad_inches=0., transparent=True, dpi=dpi)
        if outdir and sample:
            os.makedirs(outdir, exist_ok=True)
            plt.savefig(os.path.join(outdir, '{}.pie_plot.svg'.format(sample)), format='svg')
        return base64.b64encode(bytes_io.getvalue()).decode('UTF-8')
    finally:
        # Always close the figure, even on a raised ValueError -- otherwise, now that main() no
        # longer aborts on the first per-sample failure, repeated failures across many samples
        # would leak an unbounded number of open matplotlib figures in one process.
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Build phylogenetic pie plots from a tree file and cluster abundances file.'
    )
    parser.add_argument('tree_file', help='Path to tree file (BuildTree TSV or simple edge list)')
    parser.add_argument('abundances_file', help='Path to cluster abundances / constrained CCF file')
    parser.add_argument('-o', '--outdir', default='pie_plots', help='Output directory for SVG pie plots (default: pie_plots)')
    parser.add_argument('--tree-format', choices=['auto', 'buildtree', 'simple'], default='auto',
                        help='Tree file format (default: auto-detect)')
    parser.add_argument('--abundances-format', choices=['auto', 'constrained_ccf', 'simple'], default='auto',
                        help='Abundances file format (default: auto-detect)')
    parser.add_argument('--samples', nargs='*', default=None,
                        help='Only plot these sample IDs (default: all samples in abundances file)')
    args = parser.parse_args()

    tree = load_tree(args.tree_file, tree_format=args.tree_format)
    abundances = load_abundances(args.abundances_file, format=args.abundances_format)
    if not abundances:
        print('No abundances found in', args.abundances_file, file=sys.stderr)
        sys.exit(1)

    sample_ids = args.samples if args.samples is not None else sorted(abundances.keys())
    failed_samples = []
    for s in sample_ids:
        if s not in abundances:
            print('Sample "{}" not in abundances file, skipping.'.format(s), file=sys.stderr)
            failed_samples.append(s)
            continue
        try:
            make_pie_plot_standalone(tree, abundances[s], outdir=args.outdir, sample=s)
            print('Wrote {}.pie_plot.svg'.format(s))
        except ValueError as e:
            # Don't let one bad sample -- e.g. a hand-edited tree that's only inconsistent with
            # this particular sample's abundance data -- block every other sample from being
            # plotted. Report it and keep going; still exit non-zero at the end so scripted/CI
            # callers can detect a partial failure.
            print('Error for sample "{}": {}'.format(s, e), file=sys.stderr)
            failed_samples.append(s)
            continue

    if failed_samples:
        print('Done with errors. {}/{} samples failed: {}. Plots for the rest are in {}'.format(
            len(failed_samples), len(sample_ids), ', '.join(map(str, failed_samples)), args.outdir),
            file=sys.stderr)
        sys.exit(1)

    print('Done. Plots in {}'.format(args.outdir))


if __name__ == '__main__':
    main()
