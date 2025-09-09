#!/usr/bin/env python
import argparse
import pysam
import os
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


def validate_directory(path):
    path = os.path.abspath(path)
    if not os.path.exists(path):
        os.mkdir(path, mode=0o777)
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid directory.")
    return path


def main(args):
    start_time = time.time()
    memory = str(args.memory) + 'G'
    print(f"- Conversion started, converting files: {', '.join(args.files)}")
    fn_outs = []
    for _file in args.files:
        fn_outs.append(get_clean_fn(_file))
    duplicates = list(set([item for item, count in Counter(fn_outs).items() if count > 1]))
    import ipdb; ipdb.set_trace()
    if duplicates:
        raise ValueError('Duplicate file basenames, will cause collbering of data.\n' + 
                         f'Please re-name your files to have distinct basenames, duplicated: {", ".join(duplicates)}')
    for align_file in args.files:
        this_mode = get_mode(align_file)
        output_fn = get_clean_fn(os.path.basename(align_file))
        output_fp = os.path.join(args.output, output_fn)
        if not this_mode:
            raise ValueError(f'Unable to determine if cram/bam/sam: {align_file}')
        with pysam.AlignmentFile(align_file, this_mode) as alignment:
            print(f'\t > Converting {RED}{align_file}{RESET} to {GREEN}{output_fp}{RESET}')
            pysam.sort("-@", str(args.threads), "-m", memory, "-o", output_fp, align_file)
            pysam.index(output_fp)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'- Convertsion completed, elapsed_time {elapsed_time:.4f} seconds')

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
    main(parser.parse_args())