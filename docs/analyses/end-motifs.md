# <code>finaletoolkit <b>end-motifs</b></code>

## 1. About

> **Measures frequency of k-mer 5' end motifs.**

Nucleases do not cut cfDNA at random positions — the enzymes fragmenting DNA in plasma have sequence preferences, so the bases immediately at a fragment's 5' end are informative about which nuclease produced it. This step reads the reference sequence at every fragment end and tallies the frequency of each 4-mer.

The resulting 256-entry frequency vector is a cfDNA feature in its own right, and it is the direct input to [`mds`](mds.md), which collapses it to a single motif diversity score.

## 2. As the pipeline runs it

```text
finaletoolkit end-motifs bams/{sample}.sorted.bam <ref2bit> \
    -q 30 \
    -k 4 \
    -min <--fragment-minimum> \
    -max <--fragment-maximum> \
    -o end_motifs/{sample}_endmotif.tsv \
    -w <threads> \
    -v
```

### 2.1 Arguments

  `input_file`
> **The analysis BAM.**
> *value:* `bams/{sample}.sorted.bam`

---
  `refseq_file`
> **2bit reference sequence.**
> *value:* the `ref2bit` file for the selected `--genome` build
>
> The motif is read from the *reference*, not from the read's own bases, so sequencing errors at the fragment end do not perturb the tally. This is why the 2bit file must be the same assembly the BAM was aligned to — see [Reference contig filter](../pipeline/contig-filter.md) for how that is enforced on the BAM input path.

---
  `-k`
> **k-mer length.**
> *value:* `4`
>
> 4-mers give 256 possible motifs, which is the convention in the cfDNA end-motif literature and what the motif diversity score is normalized against. Changing `k` would change the number of categories and therefore the scale of the MDS.

---
  `-q, --quality-threshold`
> **Minimum mapping quality.**
> *value:* `30`

---
  `-min, --min-length` / `-max, --max-length`
> **Fragment length window.**
> *value:* `--fragment-minimum` / `--fragment-maximum` (defaults 50 / 500)

!!! note "Both strands are counted"

    FinaleToolkit's `-B, --no-both-strands` flag restricts the tally to one strand, with `-n` selecting the negative strand. The pipeline passes **neither**, so the default applies and the 5' ends of *both* strands contribute. Each fragment therefore contributes two motifs — one from each end — which doubles the counts available and avoids any strand-specific bias in the frequency vector.

    Watch for the flag-name difference if you compare this to [`interval-end-motifs`](interval-end-motifs.md): there the same `-B` is spelled `--single-strand`.

## 3. Output

`end_motifs/{sample}_endmotif.tsv` — a two-column TSV of each 4-mer and its frequency, 256 rows:

```text
AAAA	0.0083
AAAC	0.0041
AAAG	0.0052
...
```

This file is consumed directly by [`mds`](mds.md).

## 4. Controlling it

| What | How |
|------|-----|
| Fragment length window | `--fragment-minimum`, `--fragment-maximum` |
| Reference sequence | the `ref2bit` entry of `config/genome.json` for the selected build |
| Threads / memory / walltime | the `end_motifs` entry of `config/cluster.json` |

```json
"end_motifs": {
    "threads": 24,
    "mem": "4G",
    "time": "08:00:00",
    "partition": "norm"
}
```

The `-k 4` k-mer length is fixed in the workflow rather than exposed as a `run` option, because [`mds`](mds.md) downstream assumes a 256-motif vector.

## 5. Requires

`ref2bit` for the selected genome build. Both bundled builds provide it.

## 6. Related

- [`interval-end-motifs`](interval-end-motifs.md) — the same measurement per genomic window instead of genome-wide
- [`mds`](mds.md) — collapses this frequency vector to a single diversity score
