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


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/merge_cda_features.py"
HELPER_SOURCE = ROOT / "workflow/scripts/cfdnaanalyzer_dense_csv_merge.cpp"
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
        self.assertEqual(output.read_text(), "sample,label,a\n")

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


class DenseHelperTests(ScratchCase):
    def setUp(self):
        super().setUp()
        self.helper = self.root / "cda_dense_merge"
        subprocess.run([
            "/usr/bin/g++", "-O2", "-std=c++17", "-Wall", "-Wextra",
            "-pedantic", str(HELPER_SOURCE), "-o", str(self.helper),
        ], check=True)

    def run_schema(self, sources):
        listing = self.root / "sources.tsv"
        listing.write_text("".join(
            "{}\t{}\n".format(sample, path) for sample, path in sources
        ))
        schema = self.root / "schema.txt"
        header = self.root / "header.csv"
        subprocess.run([
            str(self.helper), "schema", "--list", str(listing),
            "--schema", str(schema), "--header", str(header),
        ], check=True)
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
        subprocess.run([
            str(self.helper), "row", "--source", str(second),
            "--sample", "s2", "--schema", str(schema), "--output", str(row),
        ], check=True)
        self.assertEqual(row.read_text(), "s2,0,,30,40\n")

    def test_header_only_source_yields_no_row(self):
        source = self.root / "empty.csv"
        source.write_text("sample,label,a\n")
        schema, _ = self.run_schema([("s1", source)])
        row = self.root / "empty.row.csv"
        subprocess.run([
            str(self.helper), "row", "--source", str(source),
            "--sample", "s1", "--schema", str(schema), "--output", str(row),
        ], check=True)
        self.assertEqual(row.read_bytes(), b"")

    def test_duplicate_columns_are_preserved(self):
        source = self.root / "duplicate.csv"
        source.write_text("sample,label,a,a,b\ns1,1,10,11,20\n")
        schema, header = self.run_schema([("s1", source)])
        self.assertEqual(schema.read_text(), "a\na.1\nb\n")
        self.assertEqual(header.read_text(), "sample,label,a,a.1,b\n")
        row = self.root / "duplicate.row.csv"
        subprocess.run([
            str(self.helper), "row", "--source", str(source),
            "--sample", "s1", "--schema", str(schema), "--output", str(row),
        ], check=True)
        self.assertEqual(row.read_text(), "s1,1,10,11,20\n")

    def test_unsorted_header_is_rejected(self):
        source = self.root / "unsorted.csv"
        source.write_text("sample,label,b,a\ns1,1,10,11\n")
        listing = self.root / "sources.tsv"
        listing.write_text("s1\t{}\n".format(source))
        result = subprocess.run([
            str(self.helper), "schema", "--list", str(listing),
            "--schema", str(self.root / "bad-schema"),
            "--header", str(self.root / "bad-header"),
        ], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not sorted", result.stderr)


if __name__ == "__main__":
    unittest.main()
