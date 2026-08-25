#!/usr/bin/env python
"""Subset a staged BAM to the contigs of the selected reference genome.

Input BAMs are frequently aligned against a full assembly (primary contigs
plus hundreds of alt/decoy/patch scaffolds) while this pipeline's references
are "clean" primary-only builds. This script compares the BAM header against
the reference sequence dictionary and writes a new BAM containing only the
contigs the two agree on.

Comparison uses the @SQ SN (name) and LN (length) fields only. The UR field
of a .dict records the filesystem path of the FastA it was built from and the
M5 field records a checksum; neither is present in a typical aligner-produced
BAM header, and UR would differ between the reference and the input even for
the same assembly. Comparing them would produce false mismatches, so both are
ignored.

Contig classification:

  * matched          - same SN and same LN in both. Kept.
  * length mismatch  - same SN, different LN. This means the BAM was aligned
                       to a different assembly than the one selected (e.g. an
                       hg19 BAM against an hg38 reference), so it is fatal.
                       Subsetting cannot repair it.
  * missing from BAM - in the reference but not the BAM. Warned about and
                       skipped; the output simply lacks that contig.
  * extra in BAM     - in the BAM but not the reference (the alt/decoy
                       scaffolds). Dropped, which is the point of this script.

Kept contigs are emitted in BAM-header order rather than reference-dictionary
order so that a coordinate-sorted input stays coordinate-sorted: records are
appended in the same contig order the new header declares.
"""
import argparse
import os
import subprocess
import sys

import pysam


RED     = '\033[31m'
GREEN   = '\033[32m'
YELLOW  = '\033[33m'
RESET   = '\033[0m'


def parse_sequence_dictionary(dict_file):
    """Reads the @SQ records of a Picard-style sequence dictionary.
    Only SN and LN are extracted; UR and M5 are deliberately ignored (see the
    module docstring).
    :param dict_file <str>: path to a .dict file
    :return contigs <dict>: mapping of contig name to contig length
    """
    contigs = {}
    with open(dict_file) as fh:
        for line in fh:
            if not line.startswith('@SQ'):
                continue
            fields = {}
            for field in line.rstrip('\n').split('\t')[1:]:
                if ':' in field:
                    tag, value = field.split(':', 1)
                    fields[tag] = value
            if 'SN' in fields and 'LN' in fields:
                contigs[fields['SN']] = int(fields['LN'])

    if not contigs:
        raise ValueError(
            f'No @SQ records found in sequence dictionary: {dict_file}'
        )
    return contigs


def bam_contigs(bam_file):
    """Reads contig names and lengths from a BAM header, in header order.
    :param bam_file <str>: path to a BAM file
    :return contigs list[tuple(<str>, <int>)]: (name, length) in header order
    """
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        return list(zip(bam.references, bam.lengths))


def compare(reference, observed):
    """Classifies every contig as matched, length-mismatched, missing or extra.
    :param reference <dict>: reference contig name to length
    :param observed list[tuple(<str>, <int>)]: BAM contigs in header order
    :return classification <dict>: lists keyed by category
    """
    observed_lengths = dict(observed)
    matched, mismatched = [], []

    # Iterate the BAM header, not the reference, so kept contigs come out in
    # header (coordinate-sort) order.
    for name, length in observed:
        if name not in reference:
            continue
        if reference[name] == length:
            matched.append(name)
        else:
            mismatched.append((name, reference[name], length))

    return {
        'matched': matched,
        'mismatched': mismatched,
        'missing': [c for c in reference if c not in observed_lengths],
        'extra': [n for n, _ in observed if n not in reference],
    }


def write_report(report_file, args, reference, observed, classification,
                 dropped_mates=None):
    """Writes a human-readable record of the comparison.
    :param report_file <str>: output path
    :param args <argparse.Namespace>: parsed arguments
    :param reference <dict>: reference contig name to length
    :param observed list[tuple]: BAM contigs
    :param classification <dict>: output of compare()
    :param dropped_mates <int>: reads excluded for a dangling mate, if known
    """
    matched = classification['matched']
    with open(report_file, 'w') as fh:
        fh.write('# Reference contig validation\n')
        fh.write(f'bam\t{os.path.abspath(args.bam)}\n')
        fh.write(f'sequence_dictionary\t{os.path.abspath(args.dict)}\n')
        fh.write('# Comparison uses @SQ SN and LN only; UR and M5 ignored.\n')
        fh.write(f'reference_contigs\t{len(reference)}\n')
        fh.write(f'bam_contigs\t{len(observed)}\n')
        fh.write(f'matched\t{len(matched)}\n')
        fh.write(f'length_mismatches\t{len(classification["mismatched"])}\n')
        fh.write(f'missing_from_bam\t{len(classification["missing"])}\n')
        fh.write(f'extra_in_bam_dropped\t{len(classification["extra"])}\n')
        if dropped_mates is not None:
            fh.write(f'reads_dropped_dangling_mate\t{dropped_mates}\n')

        fh.write('\n# kept contigs (name, length)\n')
        for name in matched:
            fh.write(f'kept\t{name}\t{reference[name]}\n')

        if classification['mismatched']:
            fh.write('\n# length mismatches (name, reference, bam)\n')
            for name, ref_len, bam_len in classification['mismatched']:
                fh.write(f'length_mismatch\t{name}\t{ref_len}\t{bam_len}\n')

        if classification['missing']:
            fh.write('\n# reference contigs absent from the bam\n')
            for name in classification['missing']:
                fh.write(f'missing\t{name}\t{reference[name]}\n')

        if classification['extra']:
            fh.write('\n# bam contigs absent from the reference (dropped)\n')
            for name in classification['extra']:
                fh.write(f'dropped\t{name}\n')


def build_header(bam_file, keep):
    """Builds a SAM header text retaining only the kept @SQ records.
    Non-@SQ header lines (@HD, @RG, @PG, @CO) are preserved verbatim so read
    groups and provenance survive the subset.
    :param bam_file <str>: path to the input BAM
    :param keep list[<str>]: contig names to retain, in output order
    :return header <str>: SAM header text
    """
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        header = bam.header.to_dict()

    lines = []
    if 'HD' in header:
        hd = '\t'.join(f'{k}:{v}' for k, v in header['HD'].items())
        lines.append(f'@HD\t{hd}')

    # Emit @SQ in the caller's order so the header matches the record order.
    by_name = {sq['SN']: sq for sq in header.get('SQ', [])}
    for name in keep:
        if name not in by_name:
            continue
        sq = '\t'.join(f'{k}:{v}' for k, v in by_name[name].items())
        lines.append(f'@SQ\t{sq}')

    for tag in ('RG', 'PG'):
        for record in header.get(tag, []):
            fields = '\t'.join(f'{k}:{v}' for k, v in record.items())
            lines.append(f'@{tag}\t{fields}')

    for comment in header.get('CO', []):
        lines.append(f'@CO\t{comment}')

    return '\n'.join(lines) + '\n'


def mate_consistency_filter(keep):
    """Builds a samtools filter expression that drops reads whose mate lies on
    a dropped contig.

    Such reads are why a subset is not merely a header edit. When the mate's
    contig disappears, samtools cannot resolve RNEXT and silently rewrites it
    to '*' while leaving PNEXT and the mate-mapped flag untouched, producing a
    self-inconsistent record. Excluding the read entirely is the only way to
    keep the output internally consistent. In practice these are never proper
    pairs, so the downstream fragmentomics rules (which all filter on proper
    pairs) would have ignored them regardless.

    :param keep list[<str>]: contigs retained in the output
    :return expression <str>: samtools -e filter expression
    """
    names = ' || '.join(f'mrname=="{name}"' for name in keep)
    # flag.munmap covers single-end reads and pairs with an unmapped mate,
    # where RNEXT carries no contig to preserve.
    return f'({names} || flag.munmap)'


def count_reads(bam_file, threads, regions, expression=None):
    """Counts reads in the given regions, optionally under a filter.
    :param bam_file <str>: BAM to count
    :param threads <int>: threads for samtools
    :param regions list[<str>]: regions to restrict the count to
    :param expression <str>: optional samtools -e filter expression
    :return count <int>: matching read count
    """
    command = ['samtools', 'view', '-c', '-@', str(threads)]
    if expression:
        command += ['-e', expression]
    command += [bam_file] + regions
    result = subprocess.run(command, stdout=subprocess.PIPE, check=True)
    return int(result.stdout.strip())


def subset_bam(bam_file, output_file, keep, threads, header_text):
    """Writes a BAM containing only reads on the kept contigs.

    The records are streamed as headerless SAM text and concatenated onto the
    reduced header before being converted back to BAM. This matters for
    correctness: BAM records store the reference as an integer index into the
    @SQ list, so dropping @SQ entries would silently corrupt those indices.
    SAM text refers to contigs by name instead, so the final `samtools view`
    rebuilds the name-to-index mapping against the reduced header.

    :param bam_file <str>: input BAM
    :param output_file <str>: output BAM
    :param keep list[<str>]: contigs to retain, in header order
    :param threads <int>: threads for samtools
    :param header_text <str>: reduced SAM header
    :return dropped_mates <int>: reads excluded for having a dangling mate
    """
    header_file = output_file + '.header.sam'
    with open(header_file, 'w') as fh:
        fh.write(header_text)

    keep_expression = mate_consistency_filter(keep)
    dropped_mates = count_reads(
        bam_file, threads, keep, f'!{keep_expression}'
    )

    try:
        # Stream: reduced header, then headerless records for the kept
        # regions, into a single `samtools view -b`. Requesting regions
        # requires an index on the input, which stage_bams always produces.
        convert = subprocess.Popen(
            ['samtools', 'view', '-b', '-@', str(threads), '-o', output_file, '-'],
            stdin=subprocess.PIPE,
        )
        # The header must reach the stream before any record does, so it is
        # written and flushed before the record producer is started.
        with open(header_file, 'rb') as fh:
            convert.stdin.write(fh.read())
        convert.stdin.flush()

        view = subprocess.Popen(
            ['samtools', 'view', '--no-header', '-@', str(threads),
             '-e', keep_expression, bam_file] + keep,
            stdout=convert.stdin,
        )
        view_rc = view.wait()
        # Close only after the producer exits, or it sees EPIPE.
        convert.stdin.close()
        convert_rc = convert.wait()
    finally:
        if os.path.exists(header_file):
            os.remove(header_file)

    if view_rc != 0:
        raise RuntimeError(f'samtools view failed on {bam_file} (rc={view_rc})')
    if convert_rc != 0:
        raise RuntimeError(f'samtools view -b failed (rc={convert_rc})')

    pysam.index(output_file)
    # Prevent "index is older than data file" warnings downstream.
    os.utime(output_file + '.bai', None)
    return dropped_mates


def main(args):
    reference = parse_sequence_dictionary(args.dict)
    observed = bam_contigs(args.bam)
    classification = compare(reference, observed)
    matched = classification['matched']

    print(f'- Validating {args.bam} against {args.dict}')
    print(f'\t reference contigs: {len(reference)}, bam contigs: {len(observed)}')
    print(f'\t matched: {len(matched)}, dropped: {len(classification["extra"])}')

    # Written before the fatal checks below so the report survives as a
    # diagnostic when validation fails.
    if args.report:
        write_report(args.report, args, reference, observed, classification)

    if classification['mismatched']:
        detail = '\n'.join(
            f'    {name}: reference={ref_len} bam={bam_len}'
            for name, ref_len, bam_len in classification['mismatched'][:10]
        )
        raise SystemExit(
            f'{RED}Fatal: contig length mismatch between the input BAM and the '
            f'selected reference genome.{RESET}\n'
            f'{len(classification["mismatched"])} shared contig name(s) have '
            f'different lengths, which means this BAM was aligned to a '
            f'different assembly:\n{detail}\n'
            'Subsetting cannot correct this. Re-run with the matching --genome, '
            'or realign the input.'
        )

    if not matched:
        raise SystemExit(
            f'{RED}Fatal: no contigs are shared between {args.bam} and '
            f'{args.dict}.{RESET}\nThe BAM does not appear to be aligned to the '
            'selected reference genome.'
        )

    if classification['missing']:
        print(
            f'{YELLOW}Warning: {len(classification["missing"])} reference '
            f'contig(s) are absent from the BAM header and will be absent from '
            f'the output: {", ".join(sorted(classification["missing"])[:10])}'
            f'{RESET}',
            file=sys.stderr,
        )

    print(f'- Writing {len(matched)} contig(s) to {GREEN}{args.output}{RESET}')
    header_text = build_header(args.bam, matched)
    dropped_mates = subset_bam(
        args.bam, args.output, matched, args.threads, header_text
    )

    if dropped_mates:
        print(
            f'{YELLOW}Warning: dropped {dropped_mates} read(s) whose mate was '
            f'on an excluded contig; keeping them would leave RNEXT '
            f'unresolvable.{RESET}',
            file=sys.stderr,
        )

    # Rewrite the report now that the subset counts are known.
    if args.report:
        write_report(
            args.report, args, reference, observed, classification,
            dropped_mates=dropped_mates,
        )

    print(f'- Wrote {GREEN}{args.output}{RESET}')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Validate a BAM header against a reference sequence '
                    'dictionary and subset the BAM to the shared contigs'
    )
    parser.add_argument('--bam', required=True, help='input (staged) BAM')
    parser.add_argument(
        '--dict',
        required=True,
        help='reference sequence dictionary (.dict) for the selected genome',
    )
    parser.add_argument('--output', required=True, help='output subset BAM')
    parser.add_argument(
        '--report',
        help='optional path for the validation report',
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='threads for samtools',
    )
    main(parser.parse_args())
