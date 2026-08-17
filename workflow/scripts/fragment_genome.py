#!/usr/bin/env python3
"""
Fragment sequences from a FASTA file into discrete fragments and output as BED.

Uses only the Python standard library.
The final fragment of each contig is truncated to end exactly at the contig length.

Usage:
    python fragment_fasta.py <input.fasta> <fragment_size> [-o output.bed]
"""

import argparse
import sys


def parse_fasta_lengths(fasta_path):
    """Stream a FASTA file and yield (sequence_name, sequence_length) tuples."""
    name = None
    length = 0
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                # First whitespace-delimited token after '>' is the sequence name
                name = line[1:].split()[0] if len(line) > 1 else ""
                length = 0
            else:
                length += len(line)
        if name is not None:
            yield name, length


def fragment_sequence(name, length, fragment_size):
    """
    Yield (chrom, start, end) BED intervals for discrete fragments.

    Fragments tile the contig from position 0 in steps of `fragment_size`.
    The final fragment ends exactly at `length` (may be shorter than
    `fragment_size` if length is not a multiple of fragment_size).
    """
    start = 0
    while start < length:
        end = start + fragment_size
        end = min(end, length)  # final fragment ends at contig length
        yield name, start, end
        start = end


def main():
    parser = argparse.ArgumentParser(
        description="Fragment FASTA sequences into a BED file of discrete fragments "
        "(stdlib only). Final fragment of each contig ends at contig length."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument(
        "fragment_size", type=int, help="Fragment size in bp (positive int)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output BED file (default: stdout)",
    )
    args = parser.parse_args()

    if args.fragment_size <= 0:
        parser.error("fragment_size must be a positive integer")

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        for name, length in parse_fasta_lengths(args.fasta):
            if length == 0:
                continue
            for chrom, start, end in fragment_sequence(
                name, length, args.fragment_size
            ):
                out.write(f"{chrom}\t{start}\t{end}\n")
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
