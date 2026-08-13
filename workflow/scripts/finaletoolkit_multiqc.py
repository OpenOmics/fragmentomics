#!/usr/bin/env python
"""Turn the finaletoolkit results into MultiQC sections.

MultiQC ships a module for every tool the QC rules run - samtools, FastQC,
fastp - but none for finaletoolkit, so left alone it walks straight past the
fragmentomics results and reports only the alignment QC. This script converts
those results into MultiQC *custom content*, which MultiQC does parse
natively, so the fragmentation features end up in the same report as the
alignment metrics instead of in a directory nobody opens.

Each output is a `*_mqc.yaml` document holding its own plot configuration
alongside its data, which is what lets the aggregate report be produced
without a project-level `multiqc_config.yaml`. YAML rather than JSON because
the line plots are keyed by a numeric x position (fragment length, distance
from the TSS) and YAML preserves numeric mapping keys, while JSON would turn
every one of them into a string. All of them share a `parent_id`, so they are
grouped under a single "Fragmentomics" heading in the report rather than
scattered among the samtools sections.

The sections written, from the outputs of the rule named in brackets:

  * fragmentomics_general_stats  [several]  headline metrics added as columns
                                           of MultiQC's general statistics
                                           table
  * fragmentomics_metrics        [several]  the same per-sample metrics, plus
                                           the ones too detailed for the
                                           general statistics table, as a
                                           standalone table
  * fragment_length_distribution [frag_length_bins]  fragment length histogram
                                           as a percentage of each sample's
                                           fragments, so libraries of
                                           different depth overlay
  * interval_fragment_length     [frag_length_intervals]  distribution of the
                                           per-interval median fragment
                                           length, i.e. how uniform the
                                           fragmentation is along the genome
  * end_motif_frequency          [end_motifs]        the most frequent 5'
                                           end motifs
  * interval_end_motif_frequency [interval_end_motifs]    the same, over the
                                           target intervals only
  * delfi_fragmentation_profile  [delfi]    short/long ratio along the genome
  * tss_wps_profile              [agg_wps]  aggregate WPS around the TSS
  * tss_adjusted_wps_profile     [agg_adjust_wps]    the same, after
                                           adjustment
  * tss_cleavage_profile         [agg_cleavage_profile]  aggregate cleavage
                                           proportion around the TSS

Only the inputs it is given are read, so the caller (the
`finaletoolkit_multiqc` rule) decides which sections exist: the rules that
need a reference artifact the selected genome build does not ship never run,
and their outputs are then not passed here.

Sample names are taken from the output filenames rather than from any column
inside them, which is what keeps them identical to the names MultiQC derives
for the samtools and FastQC reports of the same sample - if they differed,
the general statistics table would show two half-empty rows per sample
instead of one.
"""
import argparse
import os
import re
import sys

import pandas as pd
import yaml


RED     = '\033[31m'
GREEN   = '\033[32m'
YELLOW  = '\033[33m'
RESET   = '\033[0m'

# Suffixes the finaletoolkit rules append to a sample name. The fragment
# length histogram carries the bin size in its name (frag_length_bins writes
# {sample}_frag_bin{bin_size}.tsv), so that one has to be matched as a
# pattern rather than compared literally.
FRAG_LENGTH_BINS_SUFFIX = re.compile(r'_frag_bin\d+\.tsv$')
FRAG_LENGTH_INT_SUFFIX  = '_frag_interval.bed'
MDS_SUFFIX              = '_mds.tsv'
END_MOTIFS_SUFFIX       = '_endmotif.tsv'
INTERVAL_END_MOTIFS_SUFFIX = '_endmotif_interval.tsv'
COVERAGE_SUFFIX         = '_coverage.bed'
DELFI_SUFFIX            = '_delfi.bed'
WPS_AGGR_SUFFIX         = '_wps_out_tss_aggr.wig'
ADJUSTED_WPS_SUFFIX     = '_wps_out_tss_adj_aggr.wig'
CLEAVAGE_AGGR_SUFFIX    = '_cleavage_profile_aggr.wig'

# Rows read at a time from the per-interval end motif tables. Those hold 256
# frequency columns per interval and run to over a gigabyte per sample, which
# is more than the step's memory allocation, but the summary is a running
# weighted sum, so they are read in chunks instead of all at once.
INTERVAL_CHUNK_ROWS = 100000

# Width of the bins the per-interval median fragment lengths are counted into,
# in bp. The per-interval medians span a range of tens of bp, so binning them
# finely enough to see the shape without one point per possible half-base is
# what makes the interval distribution readable.
INTERVAL_LENGTH_BIN = 5

# Points kept per line plot. An aggregate wig is one value per base over the
# whole TSS window (4,000 of them) and a DELFI bed is one row per 5 kb bin
# (several hundred thousand), which is far more than a plot can resolve and
# enough to make the report slow to open. Both are averaged down into this
# many evenly sized groups. 500 is chosen so the DELFI profile lands near the
# ~5 Mb bin scale the method is normally read at, and so a TSS profile keeps
# single-nucleosome detail (one point per 8 bp over a 4 kb window).
MAX_PLOT_POINTS = 500

# End motifs shown in the bar plot, ranked by mean frequency across samples.
# finaletoolkit reports all 256 4-mers; the full set is unreadable as a bar
# plot and is already on disk for anyone who wants it.
TOP_END_MOTIFS = 12

# Fragments shorter than this are counted as "short" for the reported short
# fraction. 150 bp is the split DELFI uses between short and long fragments.
SHORT_FRAGMENT_MAX = 150

# Shared by every section written, so MultiQC groups them together under one
# heading instead of interleaving them with the samtools sections.
PARENT_ID = 'fragmentomics'
PARENT_NAME = 'Fragmentomics'
PARENT_DESCRIPTION = (
    'Fragmentation features computed by '
    '<a href="https://finaletoolkit.readthedocs.io" target="_blank">'
    'finaletoolkit</a> from the analysis BAM of each sample. These are '
    'pipeline results rather than quality metrics, but they are summarized '
    'here so that a sample whose fragmentation profile is an outlier can be '
    'spotted next to the alignment QC that would explain it.'
)


def sample_name(path, suffix):
    """Derives a sample name from a finaletoolkit output path by removing the
    suffix that produced it, i.e. mds/{sample}_mds.tsv -> {sample}.
    :param path <str>: path to a per-sample finaletoolkit output
    :param suffix <str|re.Pattern>: literal suffix, or pattern matching it
    :return sample <str>: sample identifier
    """
    name = os.path.basename(path)
    if isinstance(suffix, re.Pattern):
        return suffix.sub('', name)
    if name.endswith(suffix):
        return name[:-len(suffix)]
    return name


def warn(message):
    """Reports a recoverable problem with one input file. Every reader below
    keeps going after one, so a single unreadable or empty output costs the
    report that sample's section rather than failing the run at its last step.
    :param message <str>: what could not be used, and why
    """
    print(f'{YELLOW}Warning: {message}{RESET}', file=sys.stderr)


def finite(value):
    """Tests whether a value can go into a report. NaN and +/-inf are dropped
    rather than written out, because neither survives YAML round-tripping into
    MultiQC as a number.
    :param value <float>: value to test
    :return usable <bool>: True when the value is finite
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return False
    return value == value and value not in (float('inf'), float('-inf'))


def weighted_median(values, weights):
    """Median of a distribution given as (value, count) pairs, i.e. the
    smallest value whose cumulative count reaches half of the total.
    :param values <pd.Series>: distinct values, e.g. fragment lengths
    :param weights <pd.Series>: how many observations each value has
    :return median <float>: weighted median, NaN if there are no observations
    """
    frame = pd.DataFrame({'value': values, 'weight': weights}).sort_values('value')
    total = frame['weight'].sum()
    if not total:
        return float('nan')
    reached = frame['weight'].cumsum() >= total / 2.0
    return float(frame.loc[reached, 'value'].iloc[0])


def downsample(values, start=0, step=1):
    """Averages a dense profile into at most MAX_PLOT_POINTS groups, keyed by
    the position at the centre of each group.
    :param values <list[float]>: profile values, in position order
    :param start <int>: position of the first value
    :param step <int>: distance between consecutive values
    :return profile <dict>: position -> mean of the values around it
    """
    series = pd.Series(values, dtype=float).reset_index(drop=True)
    if series.empty:
        return {}

    # Ceiling division: the smallest group size that fits the series into
    # MAX_PLOT_POINTS groups. Groups are taken in order, so the last one is
    # short whenever the length does not divide evenly.
    group_size = max(1, -(-len(series) // MAX_PLOT_POINTS))
    profile = {}
    for group_id, group in series.groupby(series.index // group_size):
        mean = group.mean()
        if not finite(mean):
            continue
        centre = start + (group_id * group_size + (len(group) - 1) / 2.0) * step
        profile[int(round(centre))] = round(float(mean), 6)
    return profile


def read_frag_length_bins(paths):
    """Reads the fragment length histograms written by `frag_length_bins`,
    each a tab-delimited file of (min, max, count) rows over one sample's
    fragments.

    The distribution is returned as a percentage of the sample's fragments
    rather than as raw counts, so that samples sequenced to different depths
    can be read off the same axis - the shape of the curve is the point here,
    not its height. Lengths are reported at the midpoint of their bin, which
    is the length itself at the default bin size of 1.
    :param paths list[<str>]: per-sample histogram TSVs
    :return (distribution <dict>, metrics <dict>): sample -> length ->
        percentage of fragments, and sample -> summary metrics
    """
    distribution, metrics = {}, {}

    for path in paths:
        sample = sample_name(path, FRAG_LENGTH_BINS_SUFFIX)
        try:
            frame = pd.read_csv(path, sep='\t')
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read fragment lengths from {path}: {err}')
            continue

        if not {'min', 'max', 'count'}.issubset(frame.columns):
            warn(f'{path} is not a frag-length-bins TSV, skipping')
            continue

        lengths = (frame['min'] + frame['max']) / 2.0
        counts = pd.to_numeric(frame['count'], errors='coerce').fillna(0.0)
        total = float(counts.sum())
        if not total:
            warn(f'{path} counted no fragments, skipping sample {sample}')
            continue

        distribution[sample] = {
            int(round(length)): round(count / total * 100.0, 6)
            for length, count in zip(lengths, counts)
        }
        short = float(counts[lengths < SHORT_FRAGMENT_MAX].sum())
        metrics[sample] = {
            'fragments': int(total),
            'mean_fragment_length': round(float((lengths * counts).sum() / total), 2),
            'median_fragment_length': round(weighted_median(lengths, counts), 1),
            'short_fragment_percent': round(short / total * 100.0, 2),
        }
        print(f'\t > Read {GREEN}{path}{RESET} ({len(frame)} bins, {int(total)} fragments)')

    return distribution, metrics


def read_frag_length_intervals(paths):
    """Reads the per-interval fragment length statistics written by
    `frag_length_intervals`, one row per interval with the mean, median,
    spread and count of the fragments it caught and the fraction of them
    shorter than 150 bp.

    Two things are taken from it. The first is the distribution of the
    per-interval *median* length across intervals, which says something the
    genome-wide histogram cannot: whether the fragmentation is uniform along
    the genome or varies from region to region. The second is the cohort-level
    summary in the metrics table, weighted by each interval's fragment count
    so that a sparsely covered interval does not carry the same weight as a
    deeply covered one.

    Intervals that caught no fragments are dropped. finaletoolkit writes -1 in
    every statistic column for those, which would otherwise be read as a
    fragment length.
    :param paths list[<str>]: per-sample per-interval fragment length BEDs
    :return (distribution <dict>, metrics <dict>): sample -> median length ->
        percentage of intervals, and sample -> summary metrics
    """
    distribution, metrics = {}, {}

    for path in paths:
        sample = sample_name(path, FRAG_LENGTH_INT_SUFFIX)
        try:
            frame = pd.read_csv(path, sep='\t')
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read interval fragment lengths from {path}: {err}')
            continue

        if not {'median', 'count', 's150'}.issubset(frame.columns):
            warn(f'{path} is not a frag-length-intervals BED, skipping')
            continue

        intervals = len(frame)
        counts = pd.to_numeric(frame['count'], errors='coerce').fillna(0.0)
        frame = frame[counts > 0]
        counts = counts[counts > 0]
        if frame.empty:
            warn(f'{path} has no interval with fragments, skipping sample {sample}')
            continue

        medians = pd.to_numeric(frame['median'], errors='coerce')
        short = pd.to_numeric(frame['s150'], errors='coerce').fillna(0.0)
        # An unparsable median is dropped rather than carried as NaN: a NaN
        # reaching the weighted median below would make the reported length
        # NaN as well, which is not a value MultiQC can put in a table.
        usable = medians.dropna()
        if usable.empty:
            warn(f'{path} holds no usable interval medians, skipping sample {sample}')
            continue

        # Count the intervals into fixed-width length bins, keyed by the
        # centre of the bin, and report each bin as a percentage of the
        # sample's covered intervals so samples with different numbers of
        # covered intervals overlay.
        binned = (usable / INTERVAL_LENGTH_BIN).round() * INTERVAL_LENGTH_BIN
        tally = binned.value_counts().sort_index()
        distribution[sample] = {
            int(round(length)): round(intervals_at_length / len(usable) * 100.0, 6)
            for length, intervals_at_length in tally.items()
        }

        total = float(counts.sum())
        metrics[sample] = {
            'fragment_intervals': int(len(frame)),
            'fragment_interval_percent': round(len(frame) / intervals * 100.0, 2),
            'interval_median_fragment_length': round(
                weighted_median(usable, counts[usable.index]), 1
            ),
            'interval_short_fragment_percent': round(
                float((short * counts).sum() / total * 100.0), 2
            ),
        }
        print(
            f'\t > Read {GREEN}{path}{RESET} '
            f'({len(frame)} of {intervals} intervals with fragments)'
        )

    return distribution, metrics


def read_mds(paths):
    """Reads the motif diversity scores written by `mds`, one two-row file
    per sample.
    :param paths list[<str>]: per-sample MDS TSVs
    :return metrics <dict>: sample -> {'mds': score}
    """
    metrics = {}

    for path in paths:
        sample = sample_name(path, MDS_SUFFIX)
        try:
            frame = pd.read_csv(path, sep='\t')
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read MDS from {path}: {err}')
            continue

        if 'MDS_score' not in frame.columns or frame.empty:
            warn(f'{path} holds no MDS score, skipping')
            continue

        score = pd.to_numeric(frame['MDS_score'], errors='coerce').iloc[0]
        if not finite(score):
            warn(f'{path} holds a non-numeric MDS score, skipping')
            continue

        metrics[sample] = {'mds': round(float(score), 5)}
        print(f'\t > Read {GREEN}{path}{RESET} (MDS {float(score):.5f})')

    return metrics


def top_motifs(frequencies):
    """Keeps the TOP_END_MOTIFS motifs with the highest mean frequency across
    the cohort, and keeps the same set for every sample so that the bars stay
    comparable between them.
    :param frequencies <dict>: sample -> motif -> frequency, all 256 4-mers
    :return frequencies <dict>: sample -> motif -> frequency, the top motifs
    """
    if not frequencies:
        return {}

    ranked = pd.DataFrame(frequencies).mean(axis=1).sort_values(ascending=False)
    top = list(ranked.index[:TOP_END_MOTIFS])
    return {
        sample: {
            motif: round(float(sample_frequencies[motif]), 6)
            for motif in top if motif in sample_frequencies
        }
        for sample, sample_frequencies in frequencies.items()
    }


def read_end_motifs(paths):
    """Reads the 5' end motif frequencies written by `end_motifs`, a
    headerless (motif, frequency) file covering all 256 4-mers per sample.

    Only the TOP_END_MOTIFS motifs with the highest mean frequency across the
    cohort are kept, and the same set is kept for every sample so the bars
    stay comparable.
    :param paths list[<str>]: per-sample end motif TSVs
    :return frequencies <dict>: sample -> motif -> frequency
    """
    frequencies = {}

    for path in paths:
        sample = sample_name(path, END_MOTIFS_SUFFIX)
        try:
            frame = pd.read_csv(
                path,
                sep='\t',
                header=None,
                names=['motif', 'frequency'],
                dtype={'motif': str},
            )
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read end motifs from {path}: {err}')
            continue

        frame['frequency'] = pd.to_numeric(frame['frequency'], errors='coerce')
        frame = frame.dropna(subset=['motif', 'frequency'])
        if frame.empty:
            warn(f'{path} holds no end motif frequencies, skipping')
            continue

        frequencies[sample] = dict(zip(frame['motif'], frame['frequency']))
        print(f'\t > Read {GREEN}{path}{RESET} ({len(frame)} motifs)')

    return top_motifs(frequencies)


def read_interval_end_motifs(paths):
    """Reads the per-interval 5' end motif frequencies written by
    `interval_end_motifs`, one row per interval holding the fragment count and
    the frequency of each of the 256 4-mers within that interval.

    Reduced to one frequency per motif per sample, weighted by each interval's
    fragment count, which makes it the end motif profile of the fragments that
    fell in the target intervals - the counterpart of the genome-wide profile
    from `end_motifs`, over the regions the interval file selects. Comparing
    the two sections is the point: a difference between them is a difference
    between the intervals and the rest of the genome.

    Read in chunks, since these tables are the largest output the pipeline
    produces. Intervals with no fragments carry NaN in every motif column and
    contribute nothing, being weighted by a count of zero.
    :param paths list[<str>]: per-sample per-interval end motif TSVs
    :return (frequencies <dict>, metrics <dict>): sample -> motif -> weighted
        mean frequency, and sample -> summary metrics
    """
    frequencies, metrics = {}, {}
    positions = ['contig', 'start', 'stop', 'name', 'count']

    for path in paths:
        sample = sample_name(path, INTERVAL_END_MOTIFS_SUFFIX)
        try:
            chunks = pd.read_csv(
                path,
                sep='\t',
                chunksize=INTERVAL_CHUNK_ROWS,
                dtype={'contig': str, 'name': str},
            )
            weighted, weight, intervals, covered = None, 0.0, 0, 0
            for chunk in chunks:
                if 'count' not in chunk.columns:
                    raise ValueError('no count column')
                counts = pd.to_numeric(chunk['count'], errors='coerce').fillna(0.0)
                motifs = chunk.drop(columns=positions, errors='ignore')
                intervals += len(chunk)
                covered += int((counts > 0).sum())
                weight += float(counts.sum())
                # Frequencies are per-interval fractions, so each has to be
                # scaled by the fragments behind it before the intervals can
                # be summed.
                totals = motifs.apply(pd.to_numeric, errors='coerce') \
                               .fillna(0.0) \
                               .multiply(counts, axis=0) \
                               .sum()
                weighted = totals if weighted is None else weighted.add(totals)
        except (pd.errors.EmptyDataError, OSError, ValueError) as err:
            warn(f'could not read interval end motifs from {path}: {err}')
            continue

        if weighted is None or not weight:
            warn(f'{path} has no interval with fragments, skipping sample {sample}')
            continue

        frequencies[sample] = (weighted / weight).to_dict()
        metrics[sample] = {
            'motif_intervals': covered,
            'motif_interval_percent': round(covered / intervals * 100.0, 2),
        }
        print(
            f'\t > Read {GREEN}{path}{RESET} '
            f'({covered} of {intervals} intervals with fragment ends)'
        )

    return top_motifs(frequencies), metrics


def read_coverage(paths):
    """Reads the per-interval coverage BEDs written by `coverage`, headerless
    (contig, start, stop, name, coverage) files.

    Only summary statistics are taken. The full interval x sample matrix is
    already merged into `coverage/coverage_summary.xlsx` by
    `merge_coverage_excel`, and it is far too large to put in a report.
    :param paths list[<str>]: per-sample coverage BEDs
    :return metrics <dict>: sample -> coverage summary metrics
    """
    metrics = {}

    for path in paths:
        sample = sample_name(path, COVERAGE_SUFFIX)
        try:
            frame = pd.read_csv(
                path,
                sep='\t',
                header=None,
                names=['contig', 'start', 'stop', 'name', 'coverage'],
                dtype={'contig': str, 'name': str},
            )
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read coverage from {path}: {err}')
            continue

        coverage = pd.to_numeric(frame['coverage'], errors='coerce').dropna()
        if coverage.empty:
            warn(f'{path} holds no coverage values, skipping')
            continue

        metrics[sample] = {
            'coverage_intervals': int(coverage.size),
            'mean_coverage': round(float(coverage.mean()), 4),
            'median_coverage': round(float(coverage.median()), 4),
        }
        print(f'\t > Read {GREEN}{path}{RESET} ({coverage.size} intervals)')

    return metrics


def read_delfi(paths):
    """Reads the DELFI fragmentation profiles written by `delfi`, one row per
    genomic bin with short and long fragment counts and their ratio, both raw
    and GC-corrected.

    The GC-corrected ratio is the one reported; it is the profile the method
    is interpreted from, since the raw ratio still carries the GC bias of the
    library. Rows arrive in genome order, which is what makes averaging
    consecutive bins into a plottable profile meaningful.
    :param paths list[<str>]: per-sample DELFI BEDs
    :return (profiles <dict>, metrics <dict>): sample -> bin -> mean ratio,
        and sample -> summary metrics
    """
    profiles, metrics = {}, {}

    for path in paths:
        sample = sample_name(path, DELFI_SUFFIX)
        try:
            frame = pd.read_csv(path, sep='\t')
        except (pd.errors.EmptyDataError, OSError) as err:
            warn(f'could not read DELFI ratios from {path}: {err}')
            continue

        column = 'ratio_corrected' if 'ratio_corrected' in frame.columns else 'ratio'
        if column not in frame.columns:
            warn(f'{path} is not a delfi BED, skipping')
            continue

        ratio = pd.to_numeric(frame[column], errors='coerce')
        ratio = ratio.replace([float('inf'), float('-inf')], float('nan'))
        if ratio.dropna().empty:
            warn(f'{path} holds no usable {column} values, skipping')
            continue

        profiles[sample] = downsample(ratio.tolist(), start=1)
        metrics[sample] = {
            'delfi_bins': int(ratio.dropna().size),
            'delfi_median_ratio': round(float(ratio.median()), 4),
        }
        print(f'\t > Read {GREEN}{path}{RESET} ({len(frame)} bins, {column})')

    return profiles, metrics


def read_wig(paths, suffix):
    """Reads the aggregate profiles written by the `agg-bw` rules: a fixedStep
    wig holding one value per base over the window around the TSS.

    The header line gives the position of the first value and the distance
    between them, so the profile can be keyed by distance from the TSS rather
    than by line number.
    :param paths list[<str>]: per-sample aggregate wigs
    :param suffix <str>: filename suffix to strip for the sample name
    :return profiles <dict>: sample -> distance from TSS -> mean value
    """
    profiles = {}

    for path in paths:
        sample = sample_name(path, suffix)
        try:
            with open(path) as handle:
                lines = handle.read().split('\n')
        except OSError as err:
            warn(f'could not read profile from {path}: {err}')
            continue

        header, values = lines[0], lines[1:]
        if not header.startswith('fixedStep'):
            warn(f'{path} is not a fixedStep wig, skipping')
            continue

        # e.g. "fixedStep\tchrom=.\tstart=-2000\tstep=1\tspan=4000"
        fields = dict(
            field.split('=', 1) for field in header.split()[1:] if '=' in field
        )
        try:
            start = int(fields.get('start', 0))
            step = int(fields.get('step', 1))
        except ValueError:
            warn(f'{path} has an unparsable fixedStep header, skipping')
            continue

        profile = downsample(
            [value for value in values if value.strip()],
            start=start,
            step=step,
        )
        if not profile:
            warn(f'{path} holds no usable values, skipping')
            continue

        profiles[sample] = profile
        print(f'\t > Read {GREEN}{path}{RESET} ({len(values)} positions)')

    return profiles


def merge_metrics(*metrics):
    """Merges the per-sample metric dictionaries the readers return into one
    row per sample. A sample missing from one of them keeps the metrics it
    does have rather than being dropped.
    :param metrics <dict>: sample -> metric -> value mappings
    :return merged <dict>: sample -> all of its metrics
    """
    merged = {}
    for source in metrics:
        for sample, values in source.items():
            merged.setdefault(sample, {}).update(values)
    return merged


def write_section(output_dir, section_id, section):
    """Writes one MultiQC custom-content document. The `_mqc.yaml` suffix is
    what makes MultiQC parse the file as custom content, and the parent keys
    are what group every section this script writes under a single report
    heading.
    :param output_dir <str>: directory MultiQC will scan
    :param section_id <str>: identifier and filename stem of the section
    :param section <dict>: MultiQC custom-content document, minus its id
    :return path <str>: file written
    """
    document = {
        'id': section_id,
        'parent_id': PARENT_ID,
        'parent_name': PARENT_NAME,
        'parent_description': PARENT_DESCRIPTION,
    }
    document.update(section)

    path = os.path.join(output_dir, f'{section_id}_mqc.yaml')
    with open(path, 'w') as handle:
        yaml.safe_dump(document, handle, default_flow_style=False, sort_keys=False)
    print(f'- Wrote {GREEN}{path}{RESET}')
    return path


def general_stats_section(metrics):
    """Builds the section that adds the headline fragmentation metrics as
    columns of MultiQC's general statistics table, where they sit on the same
    row as the sample's alignment metrics.

    Deliberately a subset of the metrics table below: the general statistics
    table is the first thing in the report and stays readable only if each
    tool contributes a few columns.
    :param metrics <dict>: sample -> metric -> value
    :return section <dict>: MultiQC custom-content document
    """
    columns = [
        ('median_fragment_length', {
            'title': 'Median frag',
            'description': f'Median fragment length ({PARENT_NAME.lower()})',
            'suffix': ' bp',
            'format': '{:,.0f}',
            'scale': 'BuPu',
        }),
        ('short_fragment_percent', {
            'title': '% short',
            'description': f'Fragments shorter than {SHORT_FRAGMENT_MAX} bp',
            'suffix': '%',
            'format': '{:,.1f}',
            'min': 0,
            'max': 100,
            'scale': 'RdYlGn-rev',
        }),
        ('mds', {
            'title': 'MDS',
            'description': 'Motif diversity score of the 5\' end motifs',
            'format': '{:,.4f}',
            'min': 0,
            'max': 1,
            'scale': 'GnBu',
        }),
        ('mean_coverage', {
            'title': 'Mean cov',
            'description': 'Mean normalized coverage over the target intervals',
            'format': '{:,.3f}',
            'scale': 'YlGn',
        }),
    ]
    reported = [name for name, _ in columns
                if any(name in values for values in metrics.values())]

    return {
        'plot_type': 'generalstats',
        # Attributes the columns to this pipeline step in the table's tooltips
        # and column-visibility controls, next to the samtools and FastQC ones.
        'namespace': PARENT_NAME,
        # For general statistics columns, MultiQC takes the configuration as
        # a list of single-key mappings, one per column, rather than as the
        # single mapping the other plot types use.
        'pconfig': [{name: config} for name, config in columns if name in reported],
        'data': {
            sample: {name: values[name] for name in reported if name in values}
            for sample, values in metrics.items()
        },
    }


def metrics_section(metrics):
    """Builds the standalone per-sample metrics table. This is the full set,
    including the metrics left out of the general statistics table above.
    :param metrics <dict>: sample -> metric -> value
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'Fragmentation metrics',
        'description': (
            'Per-sample summary of the fragmentation features. Fragment '
            'counts and lengths come from '
            '<code>finaletoolkit frag-length-bins</code>, the motif '
            'diversity score from <code>finaletoolkit mds</code>, coverage '
            'from <code>finaletoolkit coverage</code> (normalized and '
            'scaled, so it is comparable between samples but is not a read '
            'depth), and the short/long ratio from '
            '<code>finaletoolkit delfi</code>, GC-corrected and taken as the '
            'median over all genomic bins. The columns marked '
            '<em>intervals</em> are the equivalents restricted to the '
            'genome build\'s target intervals, from '
            '<code>finaletoolkit frag-length-intervals</code> and '
            '<code>finaletoolkit interval-end-motifs</code>, and are weighted '
            'by the number of fragments in each interval.'
        ),
        'plot_type': 'table',
        'pconfig': {
            'id': 'fragmentomics_metrics_table',
            'title': 'Fragmentation metrics',
            'namespace': PARENT_NAME,
        },
        # Column titles, order and number formats. MultiQC drops any header
        # with no data behind it, so listing all of them here is what lets a
        # run that produced only some of the features still get a table.
        'headers': {
            'fragments': {
                'title': 'Fragments',
                'description': 'Fragments counted in the length histogram',
                'format': '{:,.0f}',
                'scale': 'Blues',
            },
            'mean_fragment_length': {
                'title': 'Mean length',
                'description': 'Mean fragment length',
                'suffix': ' bp',
                'format': '{:,.1f}',
                'scale': 'BuPu',
            },
            'median_fragment_length': {
                'title': 'Median length',
                'description': 'Median fragment length',
                'suffix': ' bp',
                'format': '{:,.0f}',
                'scale': 'BuPu',
            },
            'short_fragment_percent': {
                'title': f'% < {SHORT_FRAGMENT_MAX} bp',
                'description': f'Fragments shorter than {SHORT_FRAGMENT_MAX} bp',
                'suffix': '%',
                'format': '{:,.2f}',
                'min': 0,
                'max': 100,
                'scale': 'RdYlGn-rev',
            },
            'mds': {
                'title': 'MDS',
                'description': 'Motif diversity score of the 5\' end motifs',
                'format': '{:,.4f}',
                'min': 0,
                'max': 1,
                'scale': 'GnBu',
            },
            'coverage_intervals': {
                'title': 'Intervals',
                'description': 'Intervals with a coverage value',
                'format': '{:,.0f}',
                'scale': 'Greys',
            },
            'mean_coverage': {
                'title': 'Mean coverage',
                'description': 'Mean normalized coverage over the intervals',
                'format': '{:,.3f}',
                'scale': 'YlGn',
            },
            'median_coverage': {
                'title': 'Median coverage',
                'description': 'Median normalized coverage over the intervals',
                'format': '{:,.3f}',
                'scale': 'YlGn',
            },
            'fragment_intervals': {
                'title': 'Frag intervals',
                'description': 'Target intervals that caught at least one fragment',
                'format': '{:,.0f}',
                'scale': 'Greys',
            },
            'fragment_interval_percent': {
                'title': '% intervals',
                'description': 'Target intervals with fragments, of all of them',
                'suffix': '%',
                'format': '{:,.1f}',
                'min': 0,
                'max': 100,
                'scale': 'Greens',
            },
            'interval_median_fragment_length': {
                'title': 'Median length (intervals)',
                'description': 'Median fragment length over the target '
                               'intervals, weighted by interval fragment count',
                'suffix': ' bp',
                'format': '{:,.0f}',
                'scale': 'BuPu',
            },
            'interval_short_fragment_percent': {
                'title': f'% < {SHORT_FRAGMENT_MAX} bp (intervals)',
                'description': f'Fragments shorter than {SHORT_FRAGMENT_MAX} bp '
                               'over the target intervals',
                'suffix': '%',
                'format': '{:,.2f}',
                'min': 0,
                'max': 100,
                'scale': 'RdYlGn-rev',
            },
            'motif_intervals': {
                'title': 'Motif intervals',
                'description': 'Target intervals with fragment ends in them',
                'format': '{:,.0f}',
                'scale': 'Greys',
            },
            'motif_interval_percent': {
                'title': '% intervals (motifs)',
                'description': 'Target intervals with fragment ends, of all of them',
                'suffix': '%',
                'format': '{:,.1f}',
                'min': 0,
                'max': 100,
                'scale': 'Greens',
            },
            'delfi_bins': {
                'title': 'DELFI bins',
                'description': 'Genomic bins with a short/long ratio',
                'format': '{:,.0f}',
                'scale': 'Greys',
            },
            'delfi_median_ratio': {
                'title': 'DELFI short/long',
                'description': 'Median GC-corrected short/long fragment ratio',
                'format': '{:,.3f}',
                'scale': 'PuRd',
            },
        },
        'data': metrics,
    }


def fragment_length_section(distribution):
    """Builds the fragment length distribution line plot.
    :param distribution <dict>: sample -> fragment length -> percentage
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'Fragment length distribution',
        'description': (
            'Fragment length histogram from '
            '<code>finaletoolkit frag-length-bins</code>, as a percentage of '
            'each sample\'s fragments so that samples of differing depth can '
            'be compared. cfDNA is expected to peak near 167 bp - one '
            'nucleosome plus its linker - with a shoulder near 320 bp from '
            'dinucleosomal fragments and a 10 bp periodicity below the main '
            'peak. A profile shifted long, or missing the peak, points to '
            'genomic DNA contamination from lysed white cells.'
        ),
        'plot_type': 'linegraph',
        'pconfig': {
            'id': 'fragment_length_distribution_plot',
            'title': 'Fragmentomics: fragment length distribution',
            'xlab': 'Fragment length (bp)',
            'ylab': '% of fragments',
            'ymin': 0,
            'tt_decimals': 3,
            'tt_suffix': '%',
        },
        'data': distribution,
    }


def interval_fragment_length_section(distribution):
    """Builds the per-interval fragment length line plot. Unlike the
    genome-wide histogram this is a distribution *of intervals*, not of
    fragments: how many of the target intervals have a given median fragment
    length.
    :param distribution <dict>: sample -> median length -> % of intervals
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'Fragment length across intervals',
        'description': (
            'Distribution of the per-interval median fragment length from '
            '<code>finaletoolkit frag-length-intervals</code>, counted into '
            f'{INTERVAL_LENGTH_BIN} bp bins and shown as a percentage of the '
            'intervals that caught fragments. Where the section above shows '
            'the length of the fragments, this shows how consistent that '
            'length is along the genome: a narrow peak means every region '
            'fragments alike, while a wide or multi-peaked distribution means '
            'some regions differ from the rest. Intervals with no fragments '
            'are left out, and their share of the total is in the metrics '
            'table.'
        ),
        'plot_type': 'linegraph',
        'pconfig': {
            'id': 'interval_fragment_length_plot',
            'title': 'Fragmentomics: fragment length across intervals',
            'xlab': 'Median fragment length of an interval (bp)',
            'ylab': '% of intervals',
            'ymin': 0,
            'tt_decimals': 3,
            'tt_suffix': '%',
        },
        'data': distribution,
    }


def end_motif_section(frequencies):
    """Builds the end motif frequency bar plot.
    :param frequencies <dict>: sample -> motif -> frequency
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'End motif frequency',
        'description': (
            f'The {TOP_END_MOTIFS} most frequent 5\' end motifs, ranked by '
            'mean frequency across the samples in this run, from '
            '<code>finaletoolkit end-motifs</code>. Frequencies are '
            'fractions of all fragment ends, so the bars for one sample sum '
            'to well under 1 - the remaining 4-mers of the 256 are in '
            '<code>end_motifs/{sample}_endmotif.tsv</code>. End motif usage '
            'reflects nuclease activity and is summarized as a single number '
            'by the motif diversity score.'
        ),
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'end_motif_frequency_plot',
            'title': 'Fragmentomics: end motif frequency',
            'ylab': 'Frequency of fragment ends',
            # These are frequencies of a much larger set, not parts of a
            # whole, so offer neither a count/percentage toggle nor stacking.
            'cpswitch': False,
            'stacking': None,
            'tt_decimals': 5,
        },
        'data': frequencies,
    }


def interval_end_motif_section(frequencies):
    """Builds the interval end motif frequency bar plot, the counterpart of
    end_motif_section over the target intervals.
    :param frequencies <dict>: sample -> motif -> weighted mean frequency
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'End motif frequency across intervals',
        'description': (
            f'The {TOP_END_MOTIFS} most frequent 5\' end motifs within the '
            'genome build\'s target intervals, from '
            '<code>finaletoolkit interval-end-motifs</code>, ranked by mean '
            'frequency across the samples in this run. Each interval\'s '
            'frequencies are weighted by the fragment ends in it before being '
            'combined, so this is the profile of the fragments that fell in '
            'the intervals rather than an average over intervals. Read it '
            'against the genome-wide section above: the two differing means '
            'the intervals fragment differently from the rest of the genome. '
            'The per-interval values are in '
            '<code>interval_end_motifs/{sample}_endmotif_interval.tsv</code>.'
        ),
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'interval_end_motif_frequency_plot',
            'title': 'Fragmentomics: end motif frequency across intervals',
            'ylab': 'Frequency of fragment ends',
            'cpswitch': False,
            'stacking': None,
            'tt_decimals': 5,
        },
        'data': frequencies,
    }


def delfi_section(profiles):
    """Builds the DELFI fragmentation profile line plot.
    :param profiles <dict>: sample -> bin -> mean GC-corrected ratio
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': 'DELFI fragmentation profile',
        'description': (
            'GC-corrected short/long fragment ratio along the genome, from '
            '<code>finaletoolkit delfi</code>. Bins are in genome order and '
            f'consecutive bins are averaged into at most {MAX_PLOT_POINTS} '
            'points, which puts the profile near the megabase scale the '
            'method is normally read at; the unaveraged per-bin values are '
            'in <code>delfi/{sample}_delfi.bed</code>. Healthy samples give '
            'a comparatively flat profile, while the profile of a sample '
            'carrying tumour-derived cfDNA varies more from bin to bin.'
        ),
        'plot_type': 'linegraph',
        'pconfig': {
            'id': 'delfi_fragmentation_profile_plot',
            'title': 'Fragmentomics: DELFI fragmentation profile',
            'xlab': 'Genomic bin (genome order)',
            'ylab': 'Short/long ratio (GC-corrected)',
            'tt_decimals': 4,
        },
        'data': profiles,
    }


def tss_profile_section(section_id, name, tool, description, profiles):
    """Builds one of the three aggregate TSS profile line plots. They differ
    only in what is being plotted against distance from the TSS.
    :param section_id <str>: identifier of the section
    :param name <str>: section heading
    :param tool <str>: finaletoolkit command the profile came from
    :param description <str>: what the profile shows, appended to the source
    :param profiles <dict>: sample -> distance from TSS -> mean value
    :return section <dict>: MultiQC custom-content document
    """
    return {
        'section_name': name,
        'description': (
            f'Mean of <code>{tool}</code> over all transcription start sites, '
            'aggregated by <code>finaletoolkit agg-bw</code> and averaged '
            f'into at most {MAX_PLOT_POINTS} points across the window. '
            f'{description}'
        ),
        'plot_type': 'linegraph',
        'pconfig': {
            'id': f'{section_id}_plot',
            'title': f'Fragmentomics: {name.lower()}',
            'xlab': 'Distance from TSS (bp)',
            'ylab': name,
            'tt_decimals': 4,
        },
        'data': profiles,
    }


def main(args):
    if not os.path.isdir(args.output):
        os.makedirs(args.output, exist_ok=True)

    print('- Building MultiQC sections from the finaletoolkit outputs')
    distribution, fragment_metrics = read_frag_length_bins(args.frag_length_bins)
    interval_distribution, interval_fragment_metrics = read_frag_length_intervals(
        args.frag_length_intervals
    )
    mds_metrics = read_mds(args.mds)
    end_motifs = read_end_motifs(args.end_motifs)
    interval_end_motifs, interval_motif_metrics = read_interval_end_motifs(
        args.interval_end_motifs
    )
    coverage_metrics = read_coverage(args.coverage)
    delfi_profiles, delfi_metrics = read_delfi(args.delfi)
    wps_profiles = read_wig(args.wps_aggr, WPS_AGGR_SUFFIX)
    adjusted_wps_profiles = read_wig(args.adjusted_wps_aggr, ADJUSTED_WPS_SUFFIX)
    cleavage_profiles = read_wig(args.cleavage_aggr, CLEAVAGE_AGGR_SUFFIX)

    metrics = merge_metrics(
        fragment_metrics,
        interval_fragment_metrics,
        mds_metrics,
        interval_motif_metrics,
        coverage_metrics,
        delfi_metrics,
    )

    # Sections are written only when they have data behind them. An empty
    # custom-content file is not a no-op in MultiQC: it produces an empty
    # section, or a warning, in the report.
    sections = {}
    if metrics:
        sections['fragmentomics_general_stats'] = general_stats_section(metrics)
        sections['fragmentomics_metrics'] = metrics_section(metrics)
    if distribution:
        sections['fragment_length_distribution'] = fragment_length_section(distribution)
    if interval_distribution:
        sections['interval_fragment_length'] = interval_fragment_length_section(
            interval_distribution
        )
    if end_motifs:
        sections['end_motif_frequency'] = end_motif_section(end_motifs)
    if interval_end_motifs:
        sections['interval_end_motif_frequency'] = interval_end_motif_section(
            interval_end_motifs
        )
    if delfi_profiles:
        sections['delfi_fragmentation_profile'] = delfi_section(delfi_profiles)
    if wps_profiles:
        sections['tss_wps_profile'] = tss_profile_section(
            'tss_wps_profile',
            'WPS around the TSS',
            'finaletoolkit wps',
            'The windowed protection score is positive where fragments span '
            'the window and negative where their ends fall in it, so a '
            'nucleosome-depleted region just upstream of an expressed gene\'s '
            'TSS shows as a dip with periodic peaks either side of it.',
            wps_profiles,
        )
    if adjusted_wps_profiles:
        sections['tss_adjusted_wps_profile'] = tss_profile_section(
            'tss_adjusted_wps_profile',
            'Adjusted WPS around the TSS',
            'finaletoolkit adjust-wps',
            'The same score after the long-range trend has been subtracted, '
            'which is what makes the nucleosome periodicity comparable '
            'between samples of differing coverage.',
            adjusted_wps_profiles,
        )
    if cleavage_profiles:
        sections['tss_cleavage_profile'] = tss_profile_section(
            'tss_cleavage_profile',
            'Cleavage profile around the TSS',
            'finaletoolkit cleavage-profile',
            'The proportion of fragment ends at each position, i.e. where '
            'the nuclease cut, which is highest in the accessible region '
            'around an active TSS.',
            cleavage_profiles,
        )

    if not sections:
        # Not an error: a genome build shipping none of the reference files
        # the fragmentomics rules need leaves nothing to report, and the
        # aggregate report should still be produced from the alignment QC.
        warn(
            'no finaletoolkit outputs could be summarized, the MultiQC '
            'report will hold only the alignment QC'
        )

    for section_id, section in sections.items():
        write_section(args.output, section_id, section)

    print(f'- Wrote {len(sections)} section(s) for {len(metrics)} sample(s)')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert finaletoolkit outputs into MultiQC custom '
                    'content so the fragmentation features appear in the '
                    'aggregate report'
    )
    parser.add_argument(
        '--output',
        help='directory to write the MultiQC custom-content files into, '
             'inside the directory MultiQC scans',
        required=True,
    )
    # Every input is optional: which fragmentomics rules ran depends on the
    # reference files the selected genome build provides, so the caller
    # passes whatever exists for this run.
    parser.add_argument(
        '--frag-length-bins',
        nargs='*',
        default=[],
        help='per-sample fragment length histograms from the '
             'frag_length_bins rule',
    )
    parser.add_argument(
        '--frag-length-intervals',
        nargs='*',
        default=[],
        help='per-sample per-interval fragment length statistics from the '
             'frag_length_intervals rule',
    )
    parser.add_argument(
        '--mds',
        nargs='*',
        default=[],
        help='per-sample motif diversity scores from the mds rule',
    )
    parser.add_argument(
        '--end-motifs',
        nargs='*',
        default=[],
        help='per-sample end motif frequencies from the end_motifs rule',
    )
    parser.add_argument(
        '--interval-end-motifs',
        nargs='*',
        default=[],
        help='per-sample per-interval end motif frequencies from the '
             'interval_end_motifs rule',
    )
    parser.add_argument(
        '--coverage',
        nargs='*',
        default=[],
        help='per-sample coverage BEDs from the coverage rule',
    )
    parser.add_argument(
        '--delfi',
        nargs='*',
        default=[],
        help='per-sample fragmentation profiles from the delfi rule',
    )
    parser.add_argument(
        '--wps-aggr',
        nargs='*',
        default=[],
        help='per-sample aggregate WPS profiles from the agg_wps rule',
    )
    parser.add_argument(
        '--adjusted-wps-aggr',
        nargs='*',
        default=[],
        help='per-sample aggregate adjusted WPS profiles from the '
             'agg_adjust_wps rule',
    )
    parser.add_argument(
        '--cleavage-aggr',
        nargs='*',
        default=[],
        help='per-sample aggregate cleavage profiles from the '
             'agg_cleavage_profile rule',
    )
    main(parser.parse_args())
