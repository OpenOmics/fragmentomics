# <code>finaletoolkit <b>mds</b></code>

## 1. About

> **Reads k-mer frequencies from a file and calculates a motif diversity score (MDS) using normalized Shannon entropy.**

Collapses the 256-entry end-motif frequency vector from [`end-motifs`](end-motifs.md) into a single number per sample. MDS is the **normalized Shannon entropy** of the motif distribution, as described by Jiang *et al.* (2020):

- A **high** MDS means fragment ends are spread evenly across the possible 4-mers — diverse, less sequence-specific cleavage.
- A **low** MDS means ends concentrate on a few motifs — a more sequence-specific nuclease profile.

Because it is one number per sample, MDS is the most directly comparable fragmentomics feature in the pipeline's output and a natural candidate for cohort-level comparison.

Unlike every other analysis, this step does **not** read the BAM — it reads the end-motif TSV, so it is fast and cheap.

## 2. As the pipeline runs it

```text
mds_score=$(finaletoolkit mds end_motifs/{sample}_endmotif.tsv)
```

The score is captured from standard output and written to the sample's TSV with a header line.

### 2.1 Arguments

  `file_path`
> **The end-motif frequency table.**
> *value:* `end_motifs/{sample}_endmotif.tsv`
>
> Two columns — motif and frequency — as produced by [`end-motifs`](end-motifs.md). FinaleToolkit reads from standard input when this is omitted; the pipeline passes the file explicitly.

!!! note "`mds` has no `-o` option"

    Unlike the other subcommands, `mds` writes its result to **standard output** rather than to a file. Its only options are `-s, --sep` (field separator) and `--header` (number of header rows to skip), neither of which the pipeline overrides — the defaults match the tab-delimited, headerless output of `end-motifs`. The workflow therefore captures stdout in a shell variable and formats the TSV itself.

## 3. Output

`mds/{sample}_mds.tsv` — a two-column TSV with a header, the sample name and its score:

```text
Sample	MDS_score
Sample_A	0.9214
```

One row per file, so comparing MDS across a cohort is a matter of concatenating the per-sample TSVs.

## 4. Controlling it

MDS itself has nothing to tune. What changes the score is the upstream end-motif tally:

| What | How |
|------|-----|
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` (applied by [`end-motifs`](end-motifs.md)) |
| k-mer length | fixed at `k=4` in the workflow — MDS is normalized against 256 categories |
| Threads / memory / walltime | the `mds` entry of `config/cluster.json` |

```json
"mds": {
    "threads": 16,
    "mem": "2G",
    "time": "08:00:00",
    "partition": "norm"
}
```

The allocation is generous relative to the work; the step reads a 256-row TSV and computes an entropy.

## 5. Requires

`ref2bit` for the selected genome build — not because `mds` needs it, but because [`end-motifs`](end-motifs.md) does, and this step depends on its output.

## 6. Reference

Jiang P, *et al.* Plasma DNA End-Motif Profiling as a Fragmentomic Marker in Cancer, Pregnancy, and Transplantation. *Cancer Discovery*, 2020.
