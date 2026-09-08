#!/usr/bin/env python
"""Merge per-sample cfDNAanalyzer feature matrices into project-level ones.

cfDNAanalyzer walks its BAM list serially, so the `cda_extract` rule runs it
once per sample and this script performs the gather. Each per-sample directory
holds one CSV per feature matrix:

    <sample dir>/<matrix>.csv

shaped exactly like the matrix a single all-samples cfDNAanalyzer run would
have written, but with a single data row: a leading `sample,label` column pair
followed by one column per measurement. Merging is therefore a row
concatenation aligned on column name.

Column sets are expected to agree across samples, because every feature derives
its columns from an input shared by all samples (the region BED, the bundled
site lists, the 256 end motifs, the fixed genomic bins). Where they do not, the
union is kept and the missing cells are left blank rather than dropping either
side, and the mismatch is reported.

A sample can be legitimately absent from a matrix: cfDNAanalyzer drops any
sample that fails a feature's quality control (PFE requires deep coverage over
the promoters it scores), leaving that sample's CSV with only its header. Those
samples are simply not in the merged rows, exactly as they would not have been.

With --site-lists, the named subdirectory of each sample directory is merged the
same way, file by file. This is how the nucleosome profile (NP) feature reports
per-site-list coverage: one table per Griffin transcription factor site list,
each holding one row per sample.
"""
import argparse
import os
import sys

import pandas as pd


RED     = '\033[31m'
GREEN   = '\033[32m'
YELLOW  = '\033[33m'
RESET   = '\033[0m'

# Columns every cfDNAanalyzer feature matrix leads with. The label is a
# placeholder written by the workflow (cfDNAanalyzer will not assemble its
# matrices without a label file), not a real annotation.
INDEX_COLUMNS = ['sample', 'label']


def read_matrix(csv_file):
    """Reads one per-sample feature matrix CSV. A missing or header-only file
    yields an empty frame rather than raising, so a sample that cfDNAanalyzer
    dropped for failing a feature's quality control cannot fail the merge.
    :param csv_file <str>: path to a per-sample <matrix>.csv
    :return frame <pd.DataFrame>: the sample's row(s), possibly empty
    """
    try:
        frame = pd.read_csv(csv_file, dtype={'sample': str, 'label': str})
    except pd.errors.EmptyDataError:
        print(
            f'{YELLOW}Warning: {csv_file} is empty{RESET}',
            file=sys.stderr,
        )
        return pd.DataFrame(columns=INDEX_COLUMNS)

    if frame.empty:
        print(
            f'{YELLOW}Warning: {csv_file} has no rows, cfDNAanalyzer dropped '
            f'this sample from the matrix{RESET}',
            file=sys.stderr,
        )
    return frame


def concatenate(frames, sources, label):
    """Concatenates per-sample frames into one matrix, aligning on column name
    and keeping the first frame's column order. Mismatched column sets are
    reported and unioned rather than silently trimmed, since a mismatch means
    two samples were measured over different features and the caller needs to
    know that before using the matrix.
    :param frames list[<pd.DataFrame>]: per-sample frames, in sample order
    :param sources list[<str>]: the file each frame was read from
    :param label <str>: what is being merged, for the warning messages
    :return matrix <pd.DataFrame>: one row per sample that had data
    """
    populated = [
        (frame, source) for frame, source in zip(frames, sources)
        if not frame.empty
    ]
    if not populated:
        return pd.DataFrame(columns=INDEX_COLUMNS)

    reference_columns = list(populated[0][0].columns)
    for frame, source in populated[1:]:
        if list(frame.columns) != reference_columns:
            missing = [c for c in reference_columns if c not in frame.columns]
            extra = [c for c in frame.columns if c not in reference_columns]
            print(
                f'{YELLOW}Warning: {label} columns in {source} do not match '
                f'{populated[0][1]} ({len(missing)} missing, {len(extra)} '
                f'additional). Keeping the union.{RESET}',
                file=sys.stderr,
            )

    matrix = pd.concat([frame for frame, _ in populated], ignore_index=True, sort=False)
    # Restore the leading sample,label pair in case a union reordered things.
    ordered = [c for c in INDEX_COLUMNS if c in matrix.columns]
    ordered += [c for c in matrix.columns if c not in INDEX_COLUMNS]

    return matrix[ordered]


def merge_matrix(sample_dirs, matrix_name, output_dir):
    """Merges one feature matrix across every sample directory.
    :param sample_dirs list[<str>]: per-sample cfDNAanalyzer output directories
    :param matrix_name <str>: matrix to merge, i.e. 'EM_motifs_frequency'
    :param output_dir <str>: directory the merged CSV is written to
    :return matrix <pd.DataFrame>: the merged matrix
    """
    sources = [
        os.path.join(sample_dir, f'{matrix_name}.csv')
        for sample_dir in sample_dirs
    ]
    frames = [read_matrix(source) for source in sources]
    matrix = concatenate(frames, sources, matrix_name)

    output = os.path.join(output_dir, f'{matrix_name}.csv')
    matrix.to_csv(output, index=False)
    print(
        f'\t > {GREEN}{output}{RESET} '
        f'({len(matrix)} sample(s), {max(len(matrix.columns) - len(INDEX_COLUMNS), 0)} '
        'measurement(s))'
    )

    return matrix


def merge_site_lists(sample_dirs, site_lists, output_dir):
    """Merges the per-site-list nucleosome profile tables across samples. Every
    sample scores the same site lists, so the tables are matched by filename;
    any table only some samples have is still merged, over those samples, with
    the gap reported.
    :param sample_dirs list[<str>]: per-sample cfDNAanalyzer output directories
    :param site_lists <str>: name of the site list subdirectory in each of them
    :param output_dir <str>: directory the merged tables are written to
    """
    merged_dir = os.path.join(output_dir, site_lists)
    os.makedirs(merged_dir, exist_ok=True)

    # Table names in first-appearance order, so the merged directory is
    # ordered the way the first sample's tables were.
    tables = []
    for sample_dir in sample_dirs:
        source_dir = os.path.join(sample_dir, site_lists)
        if not os.path.isdir(source_dir):
            print(
                f'{YELLOW}Warning: {source_dir} does not exist, its sample '
                f'will be absent from the merged {site_lists} tables{RESET}',
                file=sys.stderr,
            )
            continue
        for name in sorted(os.listdir(source_dir)):
            if name.endswith('.txt') and name not in tables:
                tables.append(name)

    print(f'- Merging {len(tables)} {site_lists} table(s)')
    for name in tables:
        sources = [
            os.path.join(sample_dir, site_lists, name)
            for sample_dir in sample_dirs
        ]
        present = [source for source in sources if os.path.isfile(source)]
        if len(present) != len(sources):
            print(
                f'{YELLOW}Warning: {name} is missing for '
                f'{len(sources) - len(present)} sample(s){RESET}',
                file=sys.stderr,
            )

        frames = [read_matrix(source) for source in present]
        table = concatenate(frames, present, name)
        table.to_csv(os.path.join(merged_dir, name), index=False)

    print(f'- Wrote {GREEN}{merged_dir}{RESET}')


def main(args):
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    print(
        f'- Merging {len(args.matrices)} cfDNAanalyzer feature matrix/matrices '
        f'for {len(args.sample_dirs)} sample(s)'
    )
    for matrix_name in args.matrices:
        merge_matrix(args.sample_dirs, matrix_name, args.output)

    if args.site_lists:
        merge_site_lists(args.sample_dirs, args.site_lists, args.output)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Merge per-sample cfDNAanalyzer feature matrices into '
                    'project-level feature matrices'
    )
    parser.add_argument(
        '--sample-dirs',
        nargs='+',
        help='per-sample cfDNAanalyzer output directories, in sample order',
        required=True,
    )
    parser.add_argument(
        '--matrices',
        nargs='+',
        help='feature matrices to merge, named without the .csv suffix',
        required=True,
    )
    parser.add_argument(
        '--site-lists',
        default=None,
        help='name of the per-site-list table subdirectory to merge as well, '
             'i.e. NP_site_list for the nucleosome profile feature',
    )
    parser.add_argument(
        '--output',
        help='output directory for the merged feature matrices',
        required=True,
    )
    main(parser.parse_args())
