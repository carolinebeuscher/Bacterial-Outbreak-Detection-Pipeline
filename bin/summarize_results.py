#!/usr/bin/env python3
"""Combine one sample's genotyping, QC, and identification results into a single summary row."""

import argparse
import re
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--genotyping', required=True, help='mlst output tsv')
    parser.add_argument('--fastq_metrics', required=True, help='seqkit stats -a -b output')
    parser.add_argument('--genomeQA', required=True, help='BUSCO short_summary txt')
    parser.add_argument('--contigQA', required=True, help='GUNC *maxCSS_level.tsv')
    parser.add_argument('--phylogenetic', required=True, help='shared parsnp tree svg')
    parser.add_argument('--species', required=True)
    parser.add_argument('--genus', required=True)
    parser.add_argument('--output', required=True)
    return parser.parse_args()


def find_col(columns, *keywords):
    """Return the first column whose name contains all keywords, ignoring case/punctuation.

    Tool output column names (seqkit, GUNC) drift slightly between versions
    (eg. 'Q30(%)' vs 'Q30_pct'), so match loosely instead of by exact string.
    """
    for col in columns:
        norm = re.sub(r'[^a-z0-9]', '', col.lower())
        if all(re.sub(r'[^a-z0-9]', '', kw.lower()) in norm for kw in keywords):
            return col
    return None


def warn(path, exc):
    print(f"WARNING: could not parse {path}: {exc}", file=sys.stderr)


def parse_mlst(path):
    # mlst prints one tab-separated line per queried assembly, no header:
    # file  scheme  ST  gene1(allele)  gene2(allele)  ...
    try:
        df = pd.read_csv(path, sep='\t', header=None)
        row = df.iloc[0]
        alleles = ';'.join(str(v) for v in row[3:].dropna())
        return {'mlst_scheme': row[1], 'mlst_st': row[2], 'mlst_alleles': alleles}
    except Exception as exc:
        warn(path, exc)
        return {'mlst_scheme': 'NA', 'mlst_st': 'NA', 'mlst_alleles': 'NA'}


def parse_fastq_metrics(path):
    # seqkit stats -a -b writes one row per input fastq file (R1, R2); combine into totals.
    try:
        df = pd.read_csv(path, sep=r'\s+', engine='python')
        num_seqs_col = find_col(df.columns, 'num', 'seqs')
        sum_len_col = find_col(df.columns, 'sum', 'len')
        q20_col = find_col(df.columns, 'q20')
        q30_col = find_col(df.columns, 'q30')
        gc_col = find_col(df.columns, 'gc')

        def to_num(col):
            # seqkit formats large integers with thousands separators (eg. "1,301,326"),
            # which pd.to_numeric silently turns into NaN unless they're stripped first
            if not col:
                return None
            cleaned = df[col].astype(str).str.rstrip('%').str.replace(',', '', regex=False)
            return pd.to_numeric(cleaned, errors='coerce')

        total_reads = to_num(num_seqs_col).sum() if num_seqs_col else 'NA'
        total_bases = to_num(sum_len_col).sum() if sum_len_col else 'NA'
        avg_q20 = to_num(q20_col).mean() if q20_col else 'NA'
        avg_q30 = to_num(q30_col).mean() if q30_col else 'NA'
        avg_gc = to_num(gc_col).mean() if gc_col else 'NA'

        return {
            'total_reads': total_reads,
            'total_bases': total_bases,
            'avg_q20_pct': avg_q20,
            'avg_q30_pct': avg_q30,
            'avg_gc_pct': avg_gc,
        }
    except Exception as exc:
        warn(path, exc)
        return {
            'total_reads': 'NA', 'total_bases': 'NA',
            'avg_q20_pct': 'NA', 'avg_q30_pct': 'NA', 'avg_gc_pct': 'NA',
        }


def parse_busco(path):
    # BUSCO's short_summary file embeds one line like:
    # C:98.5%[S:98.0%,D:0.5%],F:0.5%,M:1.0%,n:124
    try:
        text = open(path).read()
        match = re.search(
            r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)',
            text,
        )
        if not match:
            raise ValueError('BUSCO result line not found')
        complete, single, dup, frag, missing, n = match.groups()
        return {
            'busco_complete_pct': complete,
            'busco_single_pct': single,
            'busco_duplicated_pct': dup,
            'busco_fragmented_pct': frag,
            'busco_missing_pct': missing,
            'busco_n_groups': n,
        }
    except Exception as exc:
        warn(path, exc)
        return {
            'busco_complete_pct': 'NA', 'busco_single_pct': 'NA',
            'busco_duplicated_pct': 'NA', 'busco_fragmented_pct': 'NA',
            'busco_missing_pct': 'NA', 'busco_n_groups': 'NA',
        }


def parse_gunc(path):
    # GUNC's *maxCSS_level.tsv is a real tab-delimited table with a header row.
    try:
        df = pd.read_csv(path, sep='\t')
        row = df.iloc[0]
        pass_col = find_col(df.columns, 'pass', 'gunc')
        contam_col = find_col(df.columns, 'contamination')
        return {
            'gunc_pass': row[pass_col] if pass_col else 'NA',
            'gunc_contamination_portion': row[contam_col] if contam_col else 'NA',
        }
    except Exception as exc:
        warn(path, exc)
        return {'gunc_pass': 'NA', 'gunc_contamination_portion': 'NA'}


def main():
    args = parse_args()

    row = {'sample_id': args.sample_id, 'genus': args.genus, 'species': args.species}
    row.update(parse_mlst(args.genotyping))
    row.update(parse_fastq_metrics(args.fastq_metrics))
    row.update(parse_busco(args.genomeQA))
    row.update(parse_gunc(args.contigQA))
    row['phylogenetic_tree'] = args.phylogenetic

    pd.DataFrame([row]).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
