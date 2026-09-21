#!/usr/bin/env python
"""Bounded-memory joins for per-sample cfDNAanalyzer feature matrices.

Each input normally has one header and one sample row. Identical schemas take
a byte-streaming fast path. Mismatched schemas are unioned in first-appearance
order and aligned one sample at a time, so the gather never retains every wide
matrix in memory. The exceptionally wide regional EMR matrix is handled by the
companion C++ helper and assembled here from its header and aligned rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


CHUNK = 4 * 1024 * 1024
INDEX_COLUMNS = ["sample", "label"]
YELLOW = "\033[33m"
RESET = "\033[0m"


@dataclass(frozen=True)
class CsvInfo:
    path: str
    expected_sample: str
    header_bytes: int
    header_sha256: str
    columns: int
    data_bytes: int
    rows: int


def validate_csv_records(path, expected_sample):
    """Validate the CSV shape and return its column and data-row counts."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or header[:2] != INDEX_COLUMNS:
            raise ValueError("Header does not begin with sample,label: {}".format(path))
        row = next(reader, None)
        if row is None:
            return len(header), 0
        if next(reader, None) is not None:
            raise ValueError("Expected at most one data row: {}".format(path))
        if len(row) != len(header):
            raise ValueError(
                "Data/header field-count mismatch ({} versus {}): {}".format(
                    len(row), len(header), path
                )
            )
        if row[0] != expected_sample:
            raise ValueError(
                "Sample field {!r} does not match {!r}: {}".format(
                    row[0], expected_sample, path
                )
            )
        return len(header), 1


def scan_csv(path, expected_sample):
    """Validate a per-sample matrix without loading its wide row."""
    path = Path(path)
    if not path.is_file():
        raise ValueError("Missing per-sample matrix: {}".format(path))
    if path.stat().st_size == 0:
        return CsvInfo(str(path), expected_sample, 0, "", 0, 0, 0)

    header_hash = hashlib.sha256()
    header_bytes = 0
    data_bytes = 0
    header_prefix = bytearray()
    in_header = True

    with path.open("rb") as handle:
        while True:
            block = handle.read(CHUNK)
            if not block:
                break
            if in_header:
                newline = block.find(b"\n")
                if newline < 0:
                    header_part, data_part = block, b""
                else:
                    header_part = block[: newline + 1]
                    data_part = block[newline + 1 :]
                    in_header = False
                header_hash.update(header_part)
                header_bytes += len(header_part)
                if len(header_prefix) < 128:
                    header_prefix.extend(header_part[: 128 - len(header_prefix)])
            else:
                data_part = block

            if data_part:
                data_bytes += len(data_part)

    if in_header:
        raise ValueError("Header is not newline-terminated: {}".format(path))
    prefix = bytes(header_prefix)
    if prefix.startswith(b"\xef\xbb\xbf"):
        prefix = prefix[3:]
    if not (
        prefix.startswith(b"sample,label,")
        or prefix.startswith(b"sample,label\r\n")
        or prefix.startswith(b"sample,label\n")
    ):
        raise ValueError("Header does not begin with sample,label: {}".format(path))

    columns, rows = validate_csv_records(path, expected_sample)

    return CsvInfo(
        str(path), expected_sample, header_bytes, header_hash.hexdigest(),
        columns, data_bytes, rows,
    )


def copy_exact(source, destination, byte_count):
    remaining = byte_count
    while remaining:
        block = source.read(min(CHUNK, remaining))
        if not block:
            raise OSError("Unexpected EOF with {} bytes left".format(remaining))
        destination.write(block)
        remaining -= len(block)


def ends_in_newline(path):
    with Path(path).open("rb") as handle:
        handle.seek(-1, os.SEEK_END)
        return handle.read(1) == b"\n"


def merge_result(method, columns, populated, sources):
    return {
        "method": method,
        "columns": columns,
        "rows": len(populated),
        "missing_samples": [item.expected_sample for item in sources if not item.rows],
    }


def build_matrix_candidate(name, sources, candidate):
    """Merge one matrix to ``candidate`` and return validation metadata."""
    candidate = Path(candidate)
    candidate.parent.mkdir(parents=True, exist_ok=True)
    header_sources = [item for item in sources if item.header_bytes]
    populated = [item for item in header_sources if item.rows]
    if not header_sources:
        candidate.write_bytes(b"sample,label\n")
        return merge_result("all_empty", 2, populated, sources)

    reference = header_sources[0]
    identical = all(
        item.header_bytes == reference.header_bytes
        and item.header_sha256 == reference.header_sha256
        for item in header_sources[1:]
    )
    if identical:
        with candidate.open("wb") as destination:
            with Path(reference.path).open("rb") as source:
                copy_exact(source, destination, reference.header_bytes)
            for item in populated:
                with Path(item.path).open("rb") as source:
                    source.seek(item.header_bytes)
                    copy_exact(source, destination, item.data_bytes)
                if not ends_in_newline(item.path):
                    destination.write(b"\n")
        return merge_result("byte_stream", reference.columns, populated, sources)

    print(
        "{}Warning: {} schemas differ; aligning their column-name union one "
        "sample at a time.{}".format(YELLOW, name, RESET),
        file=sys.stderr,
    )
    union = []
    seen = set()
    for item in header_sources:
        try:
            columns = pd.read_csv(item.path, nrows=0).columns
        except pd.errors.EmptyDataError:
            continue
        for column in columns:
            if column not in seen:
                seen.add(column)
                union.append(column)
    if len(union) < 2 or union[:2] != INDEX_COLUMNS:
        raise ValueError("Column union for {} does not begin with sample,label".format(name))

    with candidate.open("w", newline="") as destination:
        csv.writer(destination, lineterminator="\n").writerow(union)
    for number, item in enumerate(populated, 1):
        frame = pd.read_csv(item.path, dtype={"sample": str, "label": str})
        if len(frame) != 1 or str(frame.iloc[0]["sample"]) != item.expected_sample:
            raise ValueError("Unexpected row while aligning {}".format(item.path))
        frame.reindex(columns=union).to_csv(
            candidate, mode="a", header=False, index=False
        )
        print(
            "aligned={} sample={}/{} id={}".format(
                name, number, len(populated), item.expected_sample
            ),
            flush=True,
        )
    return merge_result("name_aligned_one_sample_at_a_time", len(union), populated, sources)


def publish_file(candidate, output):
    """Copy from job scratch and atomically replace one declared output."""
    candidate, output = Path(candidate), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(".{}.partial.{}".format(output.name, os.getpid()))
    try:
        with candidate.open("rb") as source, partial.open("wb") as destination:
            shutil.copyfileobj(source, destination, CHUNK)
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(str(partial), str(output))
    finally:
        if partial.exists():
            partial.unlink()


def merge_matrix(args):
    if len(args.sources) != len(args.samples):
        raise ValueError("--sources and --samples must have the same length")
    infos = [
        scan_csv(source, sample)
        for source, sample in zip(args.sources, args.samples)
    ]
    candidate = Path(args.scratch) / Path(args.output).name
    result = build_matrix_candidate(args.name, infos, candidate)
    publish_file(candidate, args.output)
    print(
        "merged={} method={} rows={} columns={} bytes={} output={}".format(
            args.name, result["method"], result["rows"], result["columns"],
            Path(args.output).stat().st_size, args.output,
        )
    )


def merge_site_lists(args):
    if len(args.sample_dirs) != len(args.samples):
        raise ValueError("--sample-dirs and --samples must have the same length")
    names = sorted({
        path.name
        for sample_dir in args.sample_dirs
        for path in (Path(sample_dir) / args.site_lists).glob("*.txt")
    })
    if not names:
        raise ValueError("No {} tables were found".format(args.site_lists))

    output_dir = Path(args.output) / args.site_lists
    scratch_dir = Path(args.scratch) / args.site_lists
    output_dir.mkdir(parents=True, exist_ok=True)
    for number, name in enumerate(names, 1):
        pairs = [
            (Path(sample_dir) / args.site_lists / name, sample)
            for sample_dir, sample in zip(args.sample_dirs, args.samples)
            if (Path(sample_dir) / args.site_lists / name).is_file()
        ]
        infos = [scan_csv(path, sample) for path, sample in pairs]
        candidate = scratch_dir / name
        build_matrix_candidate("{}/{}".format(args.site_lists, name), infos, candidate)
        publish_file(candidate, output_dir / name)
        print("site_list={}/{} name={}".format(number, len(names), name), flush=True)
    print("merged_site_lists={} tables={} output={}".format(
        args.site_lists, len(names), output_dir
    ))


def write_source_list(args):
    """Write the non-empty regional-EMR inputs consumed by the C++ helper."""
    if len(args.sources) != len(args.samples):
        raise ValueError("--sources and --samples must have the same length")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for source, sample in zip(args.sources, args.samples):
            source = Path(source)
            if source.is_file() and source.stat().st_size:
                if "\t" in str(source) or "\n" in str(source):
                    raise ValueError("Unsupported source path: {}".format(source))
                handle.write("{}\t{}\n".format(sample, source))


def assemble(args):
    candidate = Path(args.scratch) / Path(args.output).name
    candidate.parent.mkdir(parents=True, exist_ok=True)
    with candidate.open("wb") as destination:
        for path in [args.header] + args.rows:
            with Path(path).open("rb") as source:
                shutil.copyfileobj(source, destination, CHUNK)
        destination.flush()
        os.fsync(destination.fileno())
    publish_file(candidate, args.output)
    print("assembled={} rows={} bytes={} output={}".format(
        args.name, len(args.rows), Path(args.output).stat().st_size, args.output
    ))


def parser():
    root = argparse.ArgumentParser(description=__doc__)
    subcommands = root.add_subparsers(dest="command")

    matrix = subcommands.add_parser("matrix")
    matrix.add_argument("--name", required=True)
    matrix.add_argument("--sources", nargs="+", required=True)
    matrix.add_argument("--samples", nargs="+", required=True)
    matrix.add_argument("--scratch", required=True)
    matrix.add_argument("--output", required=True)
    matrix.set_defaults(func=merge_matrix)

    sites = subcommands.add_parser("site-lists")
    sites.add_argument("--sample-dirs", nargs="+", required=True)
    sites.add_argument("--samples", nargs="+", required=True)
    sites.add_argument("--site-lists", default="NP_site_list")
    sites.add_argument("--scratch", required=True)
    sites.add_argument("--output", required=True)
    sites.set_defaults(func=merge_site_lists)

    sources = subcommands.add_parser("source-list")
    sources.add_argument("--sources", nargs="+", required=True)
    sources.add_argument("--samples", nargs="+", required=True)
    sources.add_argument("--output", required=True)
    sources.set_defaults(func=write_source_list)

    assembly = subcommands.add_parser("assemble")
    assembly.add_argument("--name", required=True)
    assembly.add_argument("--header", required=True)
    assembly.add_argument("--rows", nargs="+", required=True)
    assembly.add_argument("--scratch", required=True)
    assembly.add_argument("--output", required=True)
    assembly.set_defaults(func=assemble)
    return root


def main():
    root = parser()
    args = root.parse_args()
    if not hasattr(args, "func"):
        root.error("a subcommand is required")
    try:
        args.func(args)
    except (OSError, ValueError, KeyError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
