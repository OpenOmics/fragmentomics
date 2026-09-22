#!/usr/bin/env python3
"""Unit tests for the bounded-memory cfDNAanalyzer gather helpers."""

import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/merge_cda_features.py"
# cda_dense_merge is compiled into the cfDNAanalyzer image (see
# docker/cfDNAanalyzer/cda_dense_merge.cpp), so it is not importable here and
# these tests do not build it: they state its contract against the Python model
# below, which needs no toolchain. Set CDA_DENSE_MERGE to a built binary and the
# same cases run against that instead, which is how the contract gets checked
# against the implementation the jobs actually run.
HELPER_BINARY = os.environ.get("CDA_DENSE_MERGE")
SPEC = importlib.util.spec_from_file_location("merge_cda_features", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ScratchCase(unittest.TestCase):
    def setUp(self):
        scratch = os.environ.get("TMPDIR")
        self.temp = tempfile.TemporaryDirectory(dir=scratch)
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()


class OrdinaryMergeTests(ScratchCase):
    def source(self, sample, text):
        path = self.root / "{}.csv".format(sample)
        path.write_text(text)
        return MODULE.scan_csv(path, sample)

    def test_identical_headers_stream(self):
        sources = [
            self.source("s1", "sample,label,a,b\ns1,1,2,3\n"),
            self.source("s2", "sample,label,a,b\ns2,0,4,5\n"),
        ]
        output = self.root / "merged.csv"
        result = MODULE.build_matrix_candidate("x", sources, output)
        self.assertEqual(result["method"], "byte_stream")
        self.assertEqual(
            output.read_text(), "sample,label,a,b\ns1,1,2,3\ns2,0,4,5\n"
        )

    def test_name_aligned_union_preserves_first_appearance(self):
        sources = [
            self.source("s1", "sample,label,b\ns1,1,3\n"),
            self.source("s2", "sample,label,a\ns2,0,4\n"),
        ]
        output = self.root / "merged.csv"
        result = MODULE.build_matrix_candidate("x", sources, output)
        self.assertEqual(result["method"], "name_aligned_one_sample_at_a_time")
        with output.open(newline="") as handle:
            rows = list(csv.reader(handle))
        self.assertEqual(rows[0], ["sample", "label", "b", "a"])
        self.assertEqual(rows[1], ["s1", "1", "3", ""])
        self.assertEqual(rows[2], ["s2", "0", "", "4"])

    def test_empty_and_header_only_sources_are_omitted(self):
        sources = [
            self.source("s1", ""),
            self.source("s2", "sample,label,a\n"),
        ]
        output = self.root / "merged.csv"
        result = MODULE.build_matrix_candidate("PFE", sources, output)
        self.assertEqual(result["rows"], 0)
        self.assertEqual(result["missing_samples"], ["s1", "s2"])
        self.assertEqual(output.read_text(), "sample,label\n")

    def test_header_only_columns_do_not_extend_populated_schema(self):
        sources = [
            self.source("s1", "sample,label,fallback\n"),
            self.source("s2", "sample,label,measured\ns2,0,4\n"),
        ]
        output = self.root / "merged.csv"
        result = MODULE.build_matrix_candidate("x", sources, output)
        self.assertEqual(result["method"], "byte_stream")
        self.assertEqual(
            output.read_text(), "sample,label,measured\ns2,0,4\n"
        )

    def test_rejects_wrong_sample(self):
        path = self.root / "wrong.csv"
        path.write_text("sample,label,a\nother,1,2\n")
        with self.assertRaises(ValueError):
            MODULE.scan_csv(path, "expected")

    def test_name_aligned_union_quotes_column_names(self):
        sources = [
            self.source("s1", 'sample,label,"a,a"\ns1,1,2\n'),
            self.source("s2", "sample,label,b\ns2,0,3\n"),
        ]
        output = self.root / "quoted.csv"
        MODULE.build_matrix_candidate("x", sources, output)
        with output.open(newline="") as handle:
            rows = list(csv.reader(handle))
        self.assertEqual(rows[0], ["sample", "label", "a,a", "b"])
        self.assertEqual(rows[1], ["s1", "1", "2", ""])
        self.assertEqual(rows[2], ["s2", "0", "", "3"])


class HelperError(RuntimeError):
    """What cda_dense_merge exiting non-zero looks like to these tests."""


def read_source_list(path):
    """The `sample<TAB>path` rows the schema mode is handed."""
    entries = []
    for line in Path(path).read_text().splitlines():
        sample, tab, source = line.partition("\t")
        if not tab or not sample or not source:
            raise HelperError("invalid source-list row: " + line)
        entries.append((sample, source))
    if not entries:
        raise HelperError("source list is empty")
    return entries


def feature_names(path):
    """
    One source's feature columns under the helper's naming rules: the header
    must start `sample,label` and its features must already be sorted, and a
    name repeated within one header takes a `.N` suffix so that every column of
    a source stays addressable by name.
    """
    with open(path, newline="") as handle:
        header = next(csv.reader(handle), [])
    if header[:2] != ["sample", "label"]:
        raise HelperError("header does not begin with sample,label in " + str(path))
    names = []
    previous_raw = previous_emitted = None
    duplicates = 0
    for raw in header[2:]:
        if previous_raw is not None and raw < previous_raw:
            raise HelperError(
                "feature header is not sorted in {}: {} follows {}".format(
                    path, raw, previous_raw
                )
            )
        if raw == previous_raw:
            duplicates += 1
            current = "{}.{}".format(raw, duplicates)
        else:
            duplicates = 0
            current = raw
        # Resolves a header whose raw columns are already a,a,a.1: the suffix
        # the duplicate above earns is a name the next column may also carry.
        if previous_emitted is not None and current <= previous_emitted:
            current += ".1"
        if previous_emitted is not None and current <= previous_emitted:
            raise HelperError(
                "duplicate-name normalization is not monotonic in " + str(path)
            )
        names.append(current)
        previous_raw, previous_emitted = raw, current
    return names


def write_csv_row(path, fields):
    with open(path, "w", newline="") as handle:
        csv.writer(handle, lineterminator="\n").writerow(fields)


def model_schema(listing, schema_path, header_path):
    """
    The schema mode: the union of every source's feature names, ascending.

    The helper heap-merges the sorted headers so that it holds one name per
    source at a time; at these sizes a sorted set is the same answer, since
    every source's names are ascending by the time feature_names() returns.
    """
    union = set()
    for _, source in read_source_list(listing):
        union.update(feature_names(source))
    ordered = sorted(union)
    with open(schema_path, "w") as handle:
        handle.writelines(name + "\n" for name in ordered)
    write_csv_row(header_path, ["sample", "label"] + ordered)


def model_row(source, sample, schema_path, output_path):
    """
    The row mode: one sample's data row rewritten against the union schema,
    empty where the sample does not carry that feature.
    """
    features = feature_names(source)
    taken = 0
    present = []
    for feature in Path(schema_path).read_text().splitlines():
        if taken < len(features) and features[taken] < feature:
            raise HelperError(
                "source feature is absent from schema: " + features[taken]
            )
        present.append(taken < len(features) and features[taken] == feature)
        if present[-1]:
            taken += 1
    if taken < len(features):
        raise HelperError("schema ended before source header")

    with open(source, newline="") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2:
        # A feature cfDNAanalyzer dropped the sample for leaves a header and no
        # data row. An empty output is how that sample stays out of the matrix.
        Path(output_path).write_bytes(b"")
        return
    data = rows[1]
    if len(data) < 2:
        raise HelperError("missing label in " + str(source))
    if data[0] != sample:
        raise HelperError("data-row sample differs in " + str(source))
    values = data[2:]
    if len(values) < len(features):
        raise HelperError("data row ended before its header in " + str(source))
    if len(values) > len(features):
        raise HelperError(
            "data row has more fields than its header in " + str(source)
        )
    fields = data[:2]
    cursor = 0
    for bit in present:
        fields.append(values[cursor] if bit else "")
        cursor += int(bit)
    write_csv_row(output_path, fields)


def run_binary(arguments):
    result = subprocess.run(
        [HELPER_BINARY] + arguments, text=True, capture_output=True
    )
    if result.returncode != 0:
        raise HelperError(result.stderr.strip())


def dense_schema(listing, schema_path, header_path):
    if HELPER_BINARY is None:
        model_schema(listing, schema_path, header_path)
        return
    run_binary([
        "schema", "--list", str(listing),
        "--schema", str(schema_path), "--header", str(header_path),
    ])


def dense_row(source, sample, schema_path, output_path):
    if HELPER_BINARY is None:
        model_row(source, sample, schema_path, output_path)
        return
    run_binary([
        "row", "--source", str(source), "--sample", sample,
        "--schema", str(schema_path), "--output", str(output_path),
    ])


class DenseHelperTests(ScratchCase):
    """
    The contract the cda_emr_schema and cda_emr_row rules depend on, stated
    against the Python model above and, when CDA_DENSE_MERGE names a binary,
    against that binary as well.
    """

    def run_schema(self, sources):
        listing = self.root / "sources.tsv"
        listing.write_text("".join(
            "{}\t{}\n".format(sample, path) for sample, path in sources
        ))
        schema = self.root / "schema.txt"
        header = self.root / "header.csv"
        dense_schema(listing, schema, header)
        return schema, header

    def test_schema_and_rows(self):
        first = self.root / "s1.csv"
        second = self.root / "s2.csv"
        first.write_text("sample,label,a,b\ns1,1,10,20\n")
        second.write_text("sample,label,b,c\ns2,0,30,40\n")
        schema, header = self.run_schema([("s1", first), ("s2", second)])
        self.assertEqual(schema.read_text(), "a\nb\nc\n")
        self.assertEqual(header.read_text(), "sample,label,a,b,c\n")
        row = self.root / "s2.row.csv"
        dense_row(second, "s2", schema, row)
        self.assertEqual(row.read_text(), "s2,0,,30,40\n")

    def test_header_only_source_yields_no_row(self):
        source = self.root / "empty.csv"
        source.write_text("sample,label,a\n")
        schema, _ = self.run_schema([("s1", source)])
        row = self.root / "empty.row.csv"
        dense_row(source, "s1", schema, row)
        self.assertEqual(row.read_bytes(), b"")

    def test_source_list_excludes_header_only_schema_columns(self):
        header_only = self.root / "s1.csv"
        populated = self.root / "s2.csv"
        header_only.write_text("sample,label,fallback\n")
        populated.write_text("sample,label,measured\ns2,0,4\n")
        listing = self.root / "sources.tsv"
        MODULE.write_source_list(SimpleNamespace(
            sources=[str(header_only), str(populated)],
            samples=["s1", "s2"],
            output=str(listing),
        ))
        self.assertEqual(listing.read_text(), "s2\t{}\n".format(populated))
        schema = self.root / "schema.txt"
        header = self.root / "header.csv"
        dense_schema(listing, schema, header)
        self.assertEqual(schema.read_text(), "measured\n")
        self.assertEqual(header.read_text(), "sample,label,measured\n")

    def test_source_list_excludes_zero_byte_and_all_header_only_sources(self):
        zero_byte = self.root / "zero.csv"
        header_only = self.root / "header.csv"
        zero_byte.write_bytes(b"")
        header_only.write_text("sample,label,fallback\n")
        listing = self.root / "sources.tsv"
        MODULE.write_source_list(SimpleNamespace(
            sources=[str(zero_byte), str(header_only)],
            samples=["s1", "s2"],
            output=str(listing),
        ))
        self.assertEqual(listing.read_bytes(), b"")

    def test_source_list_ignores_blank_lines_after_header(self):
        source = self.root / "blank-lines.csv"
        source.write_bytes(b"sample,label,fallback\r\n\r\n\n")
        listing = self.root / "sources.tsv"
        MODULE.write_source_list(SimpleNamespace(
            sources=[str(source)], samples=["s1"], output=str(listing)
        ))
        self.assertEqual(listing.read_bytes(), b"")

    def test_source_list_rejects_malformed_header(self):
        source = self.root / "malformed.csv"
        source.write_text("sample,wrong,a\ns1,1,2\n")
        with self.assertRaises(ValueError) as raised:
            MODULE.has_data_record(source)
        self.assertIn("does not begin with sample,label", str(raised.exception))

    def test_source_list_rejects_unterminated_header(self):
        source = self.root / "unterminated.csv"
        source.write_text("sample,label,a")
        with self.assertRaises(ValueError) as raised:
            MODULE.has_data_record(source)
        self.assertIn("not newline-terminated", str(raised.exception))

    def test_source_list_scans_large_header_in_bounded_blocks(self):
        source = self.root / "large.csv"
        with source.open("wb") as handle:
            handle.write(b"sample,label,")
            handle.write(b"a" * (MODULE.CHUNK + 17))
            handle.write(b"\ns1,1,2\n")
        self.assertTrue(MODULE.has_data_record(source))

    def test_duplicate_columns_are_preserved(self):
        source = self.root / "duplicate.csv"
        source.write_text("sample,label,a,a,b\ns1,1,10,11,20\n")
        schema, header = self.run_schema([("s1", source)])
        self.assertEqual(schema.read_text(), "a\na.1\nb\n")
        self.assertEqual(header.read_text(), "sample,label,a,a.1,b\n")
        row = self.root / "duplicate.row.csv"
        dense_row(source, "s1", schema, row)
        self.assertEqual(row.read_text(), "s1,1,10,11,20\n")

    def test_quoted_column_names_round_trip(self):
        source = self.root / "quoted.csv"
        source.write_text('sample,label,"a,a",b\ns1,1,10,20\n')
        schema, header = self.run_schema([("s1", source)])
        self.assertEqual(schema.read_text(), "a,a\nb\n")
        self.assertEqual(header.read_text(), 'sample,label,"a,a",b\n')

    def test_unsorted_header_is_rejected(self):
        source = self.root / "unsorted.csv"
        source.write_text("sample,label,b,a\ns1,1,10,11\n")
        with self.assertRaises(HelperError) as raised:
            self.run_schema([("s1", source)])
        self.assertIn("not sorted", str(raised.exception))

    def test_feature_missing_from_schema_is_rejected(self):
        source = self.root / "extra.csv"
        source.write_text("sample,label,a,b\ns1,1,10,20\n")
        schema = self.root / "short-schema.txt"
        schema.write_text("b\n")
        with self.assertRaises(HelperError) as raised:
            dense_row(source, "s1", schema, self.root / "extra.row.csv")
        self.assertIn("absent from schema", str(raised.exception))

    def test_wrong_sample_is_rejected(self):
        source = self.root / "wrong.csv"
        source.write_text("sample,label,a\nother,1,10\n")
        schema, _ = self.run_schema([("s1", source)])
        with self.assertRaises(HelperError) as raised:
            dense_row(source, "s1", schema, self.root / "wrong.row.csv")
        self.assertIn("sample differs", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
