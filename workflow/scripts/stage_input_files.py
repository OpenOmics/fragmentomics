#!/usr/bin/env python
import argparse
import pysam
import os
import shutil
import subprocess
import time
import re
from collections import Counter


RED     = '\033[31m'
GREEN   = '\033[32m'
BLUE    = '\033[34m'
RESET   = '\033[0m'


def get_mode(file):
    fn = os.path.basename(file).lower()
    mode = None
    if fn.endswith('.bam'):
        mode = "rb"
    elif fn.endswith('.cram'):
        mode = "rc"
    elif fn.endswith('.sam'):
        mode = "r"
    return mode


def get_clean_fn(fn):
    fn = os.path.basename(fn)
    fn = re.sub('\\.sorted', '',  fn, flags=re.IGNORECASE)
    fn = re.sub('\\.bam', '', fn, flags=re.IGNORECASE)
    fn = re.sub('\\.cram', '', fn, flags=re.IGNORECASE)
    fn = re.sub('\\.sam', '', fn, flags=re.IGNORECASE)
    fn = fn + '.sorted.bam'
    return fn


def build_filters(min_mapq, min_baseq):
    """Builds the samtools view arguments for the pipeline's read filtering
    thresholds (the frontend's --mapscore and --baseqscore). Both are inclusive
    lower bounds, and a threshold of 0 leaves that filter out entirely.
    @param min_mapq <int>:
        Minimum mapping quality (MAPQ) to keep
    @param min_baseq <int>:
        Minimum mean base quality (Phred) to keep
    @return args list[<str>]:
        samtools view arguments, empty when no filtering was requested
    """
    args = []
    if min_mapq > 0:
        # -q skips alignments with MAPQ *smaller* than the value, i.e. it keeps
        # MAPQ >= min_mapq.
        args += ['-q', str(min_mapq)]
    if min_baseq > 0:
        # Inclusive: a mean of exactly min_baseq is kept. Note that htslib
        # stores a missing quality string ('*') as 0xff per base, so avg(qual)
        # is 255 for such a record and it passes any threshold rather than
        # being dropped. Verified against samtools 1.13.
        args += ['-e', f'avg(qual) >= {min_baseq}']
    return args


def filter_and_sort(align_file, output_fp, threads, filters):
    """Coordinate sorts an alignment file into place, applying read filters on
    the way if any were requested. Filtering is piped straight into the sort so
    no intermediate alignment file is written to disk.
    @param align_file <str>:
        Input bam/cram/sam file
    @param output_fp <str>:
        Output path for the sorted bam
    @param threads <int>:
        Threads to hand to samtools
    @param filters list[<str>]:
        samtools view arguments from build_filters()
    """
    if not filters:
        pysam.sort('-@', str(threads), '-o', output_fp, align_file)
        return

    view = ['samtools', 'view', '-b', '-@', str(threads)] + filters + [align_file]
    sort = ['samtools', 'sort', '-@', str(threads), '-o', output_fp, '-']
    viewing = subprocess.Popen(view, stdout=subprocess.PIPE)
    sorting = subprocess.Popen(sort, stdin=viewing.stdout)
    # Close this process's handle on the pipe so view is signalled if sort dies
    viewing.stdout.close()
    sort_rc = sorting.wait()
    view_rc = viewing.wait()
    if view_rc != 0 or sort_rc != 0:
        raise RuntimeError(
            f'Failed to filter and sort {align_file}: '
            f'`{" ".join(view)}` exited {view_rc}, '
            f'`{" ".join(sort)}` exited {sort_rc}'
        )
    return


def validate_directory(path):
    path = os.path.abspath(path)
    if not os.path.exists(path):
        os.mkdir(path, mode=0o777)
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid directory.")
    return path


def main(args):
    start_time = time.time()
    memory = str(int(args.memory*0.7)) + 'G' # since args.memory is total job memory, use fraction of it
    filters = build_filters(args.min_mapq, args.min_baseq)
    if filters and not shutil.which('samtools'):
        raise RuntimeError(
            'Read filtering was requested (--min-mapq/--min-baseq) but samtools '
            'was not found on PATH. Either run this inside the pipeline\'s '
            'container or pass 0 for both thresholds to skip filtering.'
        )
    print(f"- Conversion started, converting files: {', '.join(args.files)}")
    if filters:
        print(f"\t > Read filters: MAPQ >= {args.min_mapq}, "
              f"mean base quality >= {args.min_baseq}")
    else:
        print("\t > Read filters: none requested")
    fn_outs = []
    for _file in args.files:
        fn_outs.append(get_clean_fn(_file))
    duplicates = list(set([item for item, count in Counter(fn_outs).items() if count > 1]))
    if duplicates:
        raise ValueError('Duplicate file basenames, will cause collbering of data.\n' + 
                         f'Please re-name your files to have distinct basenames, duplicated: {", ".join(duplicates)}')
    idxs = []
    for align_file in args.files:
        this_mode = get_mode(align_file)
        output_fn = get_clean_fn(os.path.basename(align_file))
        output_fp = os.path.join(args.output, output_fn)
        if not this_mode:
            raise ValueError(f'Unable to determine if cram/bam/sam: {align_file}')
        with pysam.AlignmentFile(align_file, this_mode) as alignment:
            print(f'\t > Converting {RED}{align_file}{RESET} to {GREEN}{output_fp}{RESET}')
            # pysam.sort("-@", str(args.threads), "-m", memory, "-o", output_fp, align_file)
            filter_and_sort(align_file, output_fp, args.threads, filters)
            pysam.index(output_fp)
            idxs.append(output_fp + '.bai')
            alignment.close()
    end_time = time.time()
    
    for idx in idxs:
        # prevent index older than bam errors
        os.utime(idx, None)

    elapsed_time = end_time - start_time
    print(f'- Convertsion(s) completed, elapsed time: {elapsed_time:.4f} seconds')
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'A script to convert a glob of bam/sam/cram files into a sorted bam')
    parser.add_argument(
        '--files',
        nargs='+', 
        help='bam/sam/cram files to convert and sort',
        required=True
    )
    parser.add_argument(
        "--output", 
        type=validate_directory,
        help="output directory",
        required=True
    )
    parser.add_argument(
        '--threads', 
        type=int,
        default=4,
        help='threads for multiprocessing alignment file conversion'
    )
    parser.add_argument(
        '--memory',
        default=4,
        type=int,
        help='maxmimum memory utilization for alignment file conversion in integers of gigabytes (Gb)'
    )
    parser.add_argument(
        '--min-mapq',
        default=0,
        type=int,
        help='minimum mapping quality (MAPQ) to keep, i.e. keep reads with MAPQ >= this value, 0 disables'
    )
    parser.add_argument(
        '--min-baseq',
        default=0,
        type=int,
        help='minimum mean base quality (Phred) to keep, i.e. keep reads with mean base quality >= this value, 0 disables'
    )
    main(parser.parse_args())