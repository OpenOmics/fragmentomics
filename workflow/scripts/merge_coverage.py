#!/usr/bin/env python
"""Merge per-sample finaletoolkit coverage BEDs into a single Excel workbook.

Each input BED is produced by the `coverage` rule and is a headerless,
tab-delimited file with five columns:

    contig  start  stop  name  coverage

The workbook contains two sheets:

  * summary  - one row per sample with descriptive statistics over its
               intervals (interval count, total/mean/median/sd/min/max
               coverage). This is the sheet most users want.
  * coverage - the merged coverage matrix, one row per interval and one
               column per sample, outer-joined on (contig, start, stop,
               name) so intervals missing from a sample become blanks.

Excel caps a worksheet at 1,048,576 rows. When the merged matrix exceeds
that, the `coverage` sheet is skipped and the matrix is written next to the
workbook as a `.tsv` instead; the summary sheet is always written.
"""
import argparse
import os
import sys

import pandas as pd


RED     = '\033[31m'
GREEN   = '\033[32m'
YELLOW  = '\033[33m'
RESET   = '\033[0m'

# Maximum rows in an xlsx worksheet, minus one for the header row.
EXCEL_MAX_ROWS = 1048576 - 1
# Columns shared by every coverage BED, used as the join key.
INTERVAL_COLUMNS = ['contig', 'start', 'stop', 'name']
BED_COLUMNS = INTERVAL_COLUMNS + ['coverage']


def get_sample_name(bed_file, suffix='_coverage.bed'):
    """Derives a sample name from a coverage BED path. Mirrors the naming
    of the `coverage` rule, i.e. coverage/{sid}_coverage.bed -> {sid}.
    :param bed_file <str>: path to a per-sample coverage BED
    :param suffix <str>: filename suffix added by the coverage rule
    :return sample <str>: sample identifier
    """
    name = os.path.basename(bed_file)
    if name.endswith(suffix):
        name = name[:-len(suffix)]
    return name


def read_coverage_bed(bed_file):
    """Reads one coverage BED into a dataframe with a sample-named
    coverage column. Empty files yield an empty frame rather than raising,
    so one uncovered sample cannot fail the whole merge.
    :param bed_file <str>: path to a per-sample coverage BED
    :return (sample <str>, frame <pd.DataFrame>)
    """
    sample = get_sample_name(bed_file)
    try:
        frame = pd.read_csv(
            bed_file,
            sep='\t',
            header=None,
            names=BED_COLUMNS,
            dtype={'contig': str, 'name': str},
        )
    except pd.errors.EmptyDataError:
        print(
            f'{YELLOW}Warning: {bed_file} is empty, '
            f'sample {sample} will have no intervals{RESET}',
            file=sys.stderr,
        )
        frame = pd.DataFrame(columns=BED_COLUMNS)

    frame = frame.rename(columns={'coverage': sample})
    return sample, frame


def merge_coverage(bed_files):
    """Outer-joins every per-sample coverage BED on its interval columns.
    Rows are ordered by contig first-appearance across the inputs, then by
    coordinate, so the merged matrix keeps the genomic contig order of the
    source BEDs rather than a lexical one (chr2 before chr10).
    :param bed_files list[<str>]: per-sample coverage BEDs
    :return (samples list[<str>], matrix <pd.DataFrame>)
    """
    samples, matrix, contigs = [], None, []

    for bed_file in bed_files:
        sample, frame = read_coverage_bed(bed_file)
        if sample in samples:
            raise ValueError(
                'Duplicate sample name derived from coverage BEDs: '
                f'{sample}. Please rename the inputs to have distinct '
                'basenames.'
            )
        samples.append(sample)
        print(f'\t > Read {GREEN}{bed_file}{RESET} ({len(frame)} intervals)')

        for contig in frame['contig'].tolist():
            if contig not in contigs:
                contigs.append(contig)

        if matrix is None:
            matrix = frame
        else:
            matrix = matrix.merge(frame, on=INTERVAL_COLUMNS, how='outer')

    if matrix is None:
        matrix = pd.DataFrame(columns=INTERVAL_COLUMNS)

    if len(matrix):
        contig_rank = {contig: i for i, contig in enumerate(contigs)}
        matrix = (
            matrix
            .assign(_contig_rank=matrix['contig'].map(contig_rank))
            .sort_values(['_contig_rank'] + INTERVAL_COLUMNS[1:])
            .drop(columns='_contig_rank')
            .reset_index(drop=True)
        )
    return samples, matrix


def summarize(matrix, samples):
    """Builds per-sample descriptive statistics from the merged matrix.
    :param matrix <pd.DataFrame>: merged coverage matrix
    :param samples list[<str>]: sample column names in the matrix
    :return summary <pd.DataFrame>: one row per sample
    """
    rows = []
    for sample in samples:
        coverage = pd.to_numeric(matrix[sample], errors='coerce').dropna()
        rows.append({
            'sample': sample,
            'intervals': int(coverage.size),
            'total_coverage': coverage.sum() if coverage.size else 0.0,
            'mean_coverage': coverage.mean() if coverage.size else float('nan'),
            'median_coverage': coverage.median() if coverage.size else float('nan'),
            'sd_coverage': coverage.std() if coverage.size > 1 else float('nan'),
            'min_coverage': coverage.min() if coverage.size else float('nan'),
            'max_coverage': coverage.max() if coverage.size else float('nan'),
        })
    return pd.DataFrame(rows)


def main(args):
    print(f'- Merging coverage for {len(args.beds)} sample(s)')
    samples, matrix = merge_coverage(args.beds)
    summary = summarize(matrix, samples)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    too_many_rows = len(matrix) > EXCEL_MAX_ROWS
    if too_many_rows:
        sidecar = os.path.splitext(args.output)[0] + '.tsv'
        print(
            f'{YELLOW}Warning: merged matrix has {len(matrix)} rows, which '
            f'exceeds the {EXCEL_MAX_ROWS} row xlsx limit. Writing the matrix '
            f'to {sidecar} instead of a worksheet.{RESET}',
            file=sys.stderr,
        )
        matrix.to_csv(sidecar, sep='\t', index=False)

    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        summary.to_excel(writer, sheet_name='summary', index=False)
        if not too_many_rows:
            matrix.to_excel(writer, sheet_name='coverage', index=False)

    print(f'- Wrote {GREEN}{args.output}{RESET}')
    print(f'  samples: {len(samples)}, intervals: {len(matrix)}')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Merge per-sample finaletoolkit coverage BEDs into an '
                    'Excel workbook'
    )
    parser.add_argument(
        '--beds',
        nargs='+',
        help='per-sample coverage BED files from the coverage rule',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='output xlsx workbook',
        required=True,
    )
    main(parser.parse_args())
