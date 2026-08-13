# <code>fragmentomics <b>run</b></code>

## 1. About 
The `fragmentomics` executable is composed of several inter-related sub commands. Please see `fragmentomics -h` for all available options.

This part of the documentation describes options and concepts for <code>fragmentomics <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running fragmentomics pipeline. 

Setting up the fragmentomics pipeline is fast and easy! In its most basic form, <code>fragmentomics <b>run</b></code> only has *three required arguments*.

## 2. Synopsis

```text
$ fragmentomics run [--help] [--overwrite-pipeline-template] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--fragment-minimum FRAGMENT_MINIMUM] \
      [--fragment-maximum FRAGMENT_MAXIMUM] \
      [--left-tss-flank LEFT_TSS_FLANK] \
      [--right-tss-flank RIGHT_TSS_FLANK] \
      [--split-interval SPLIT_INTERVAL] \
      [--bin-size BIN_SIZE] \
      [--mapscore MAPSCORE] [--baseqscore BASEQSCORE] \
      -g {{hg38,hg19}} \
      --input INPUT [INPUT ...] \
      --output OUTPUT
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of Illumina FastQ or BAM files (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, and a reference genome build via the `--genome` argument.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input Paired-end Illumina FastQ or BAM file(s).**  
> *type: file(s)*  
> 
> One or more paired-end Illumina FastQ files **or** one or more BAM files can be provided. Only a single input type is allowed per run (i.e all inputs must either be all FastQ or BAM files). Input FastQ and BAM files **cannot** be mixed in a single run (i.e `--input *.bam *.fastq.gz`). The input type is auto-detected from the files you provide. It is worth noting that single-end FastQ files are not supported. From the command-line, each input file should be seperated by a space. Globbing is supported! This makes selecting input files easy.
> 
> ***Example:*** `--input .tests/*.bam`  
> ***Example:*** `--input .tests/*.R1.fastq.gz .tests/*.R2.fastq.gz`
>
> **How each input type is processed.** The pipeline converges both input types
> on a single coordinate-sorted, indexed BAM per sample
> (`bams/{sample}.sorted.bam`) that feeds every downstream fragmentomics
> analysis:
> - **FastQ input** runs the *full* pipeline. Paired-end reads are quality
>   filtered, aligned to the selected `--genome` build with `bwa-mem2`, and
>   reduced to properly-paired, primary, mapped reads with duplicates marked
>   (but kept) before the analysis steps run. See
>   [FastQ alignment](../pipeline/fastq-alignment.md).
> - **BAM input** runs the *full pipeline minus alignment*. Provided alignments
>   are staged (quality filtered, coordinate-sorted and indexed), stripped of
>   singletons and of the orphaned mates that per-read filtering leaves behind,
>   checked against the selected genome build's sequence dictionary, and subset
>   to the contigs the two share. See
>   [BAM normalization](../pipeline/bam-normalization.md) and
>   [Reference contig filter](../pipeline/contig-filter.md).

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output pipeline_output`

---  
  `--genome {hg38,hg19}`
> **Reference genome build.**  
> *type: string*   
> *default: hg38*  
> 
> Selects the bundled set of reference files the pipeline uses to characterize cfDNA fragmentation features. Choosing a genome build determines which chromosome sizes, 2bit reference sequence, genomic interval and TSS files, and any blacklist or gap files are used throughout the analysis. For FastQ input, it additionally selects the `bwa-mem2` reference index that reads are aligned against. For BAM input, it selects the sequence dictionary the input alignments are validated against. This argument is **required** for both FastQ and BAM input (BAM inputs still need the build to select the analysis references, and are assumed to already be aligned to it). Vaild options include: `hg38` or `hg19`.
>
> Every file each build provides, and which analysis uses it, is documented in [Reference files](../pipeline/references.md).
>  
> ***Example:*** `--genome hg38` 

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. What each analysis does with these values is documented per analysis in the [Analyses](../analyses/index.md) section.

  `--fragment-minimum FRAGMENT_MINIMUM`  
> **Minimum fragment length.**  
> *type: int*  
> *default: 50*
> 
> Minimum fragment length, in base pairs. Fragments shorter than this length are excluded from the fragmentation analyses: [coverage](../analyses/coverage.md), [fragment-length bins](../analyses/frag-length-bins.md) and [intervals](../analyses/frag-length-intervals.md), [end motifs](../analyses/end-motifs.md) and [interval end motifs](../analyses/interval-end-motifs.md), and [cleavage profiles](../analyses/cleavage-profile.md).
>
> Two analyses do not use this value: [`wps`](../analyses/wps.md) uses a fixed 120–180 bp window (the definition of L-WPS), and [`delfi`](../analyses/delfi.md) partitions short from long internally.
> 
> ***Example:*** `--min 50`

---  
  `--fragment-maximum FRAGMENT_MAXIMUM`  
> **Maximum fragment length.**  
> *type: int*  
> *default: 500*
> 
> Maximum fragment length, in base pairs. Fragments longer than this length are excluded from the same analyses listed under `--fragment-minimum`, and with the same two exceptions. Raise this if you need the dinucleosome and trinucleosome range represented in the fragment-length distribution.
> 
> ***Example:*** `--max 500`

---  
  `--left-tss-flank LEFT_TSS_FLANK`  
> **Left TSS flank size.**  
> *type: int*  
> *default: 2000*
> 
> Size, in base pairs, of the region to include to the left of each transcription start site (TSS) when computing [cleavage profiles](../analyses/cleavage-profile.md). The TSS reference file holds point positions, so a region has to be built around each one before a profile can be computed. Left and right flanks are set separately, so the window need not be symmetric.
> 
> ***Example:*** `--left-tss-flank 2000`

---  
  `--right-tss-flank RIGHT_TSS_FLANK`  
> **Right TSS flank size.**  
> *type: int*  
> *default: 2000*
> 
> Size, in base pairs, of the region to include to the right of each transcription start site (TSS) when computing [cleavage profiles](../analyses/cleavage-profile.md). Widening this extends the profile further into the gene body, where the downstream nucleosome array is most regular. Note that `cleavage_profile` is the most expensive step in the pipeline, and its cost scales with the total flank width.
> 
> ***Example:*** `--right-tss-flank 2000`

---  
  `--split-interval SPLIT_INTERVAL`  
> **WPS interval size.**  
> *type: int*  
> *default: 5000*
> 
> Interval size, in base pairs, used when computing and adjusting [window protection scores](../analyses/wps.md) around transcription start sites. Each TSS is expanded to an interval of this size, centered on the site, and WPS is computed across every base of it. The same value is passed to [`adjust-wps`](../analyses/adjust-wps.md) so the two stages agree.
>
> Setting this to `0` also disables the [`frag-length-intervals`](../analyses/frag-length-intervals.md) analysis.
> 
> ***Example:*** `--split-interval 5000`

---  
  `--bin-size BIN_SIZE`  
> **Fragment-length bin size.**  
> *type: int*  
> *default: 1*
> 
> Bin size, in base pairs, for the [fragment-length distribution](../analyses/frag-length-bins.md) histogram. At the default of `1`, every distinct fragment length gets its own row, which is what preserves the ~10 bp periodicity in the distribution. This value is also embedded in the fragment-length bin output filenames, so runs at different bin sizes do not overwrite each other.
> 
> ***Example:*** `--bin-size 1`

---  
  `--mapscore MAPSCORE`, `--m MAPSCORE`  
> **Minimum mapping quality.**  
> *type: int*  
> *default: 20*
> 
> Minimum mapping quality (MAPQ) a read must have to be kept: a read survives when its `MAPQ >= MAPSCORE`. The filter is applied with `samtools` on both input paths — on the [FastQ path](../pipeline/fastq-alignment.md) immediately after alignment, and on the [BAM path](../pipeline/bam-normalization.md) while the input alignments are staged. Pass `0` to keep every alignment regardless of mapping quality.
> 
> The same value is passed to every FinaleToolkit analysis as its `-q` threshold, so one setting governs both what is written into the analysis BAM and what the analyses read out of it. See [shared conventions](../analyses/index.md#3-shared-conventions).
> 
> ***Example:*** `--mapscore 20`

---  
  `--baseqscore BASEQSCORE`, `--b BASEQSCORE`  
> **Minimum mean base quality.**  
> *type: int*  
> *default: 20*
> 
> Minimum mean base quality (Phred) a read must have to be kept: a read survives when its mean base quality is `>= BASEQSCORE`. On the [FastQ path](../pipeline/fastq-alignment.md) this is applied by `fastp` before the reads are aligned, so a failing pair is never aligned at all; on the [BAM path](../pipeline/bam-normalization.md) it is applied to the staged alignments by `samtools` using the filter expression `avg(qual) >= BASEQSCORE`. Pass `0` to keep every read regardless of base quality.
> 
> Both paths therefore mean the same thing by this option: the *mean* quality over a read's bases, not a per-base cutoff.
> 
> On the FastQ path the mean is measured *after* [adapter trimming](../pipeline/fastq-alignment.md#34-adapter-trimming-and-read-filtering), which `fastp` performs in the same pass, so low-quality adapter read-through at the 3′ end cannot push an otherwise good pair below the threshold. Adapter trimming is not configurable by this option and is always on for FastQ input.
> 
> ***Example:*** `--baseqscore 20`

### 2.3 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running fragmentomics in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a slurm execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*  
> *default: pipeline_fragmentomics*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to `pipeline_fragmentomics`.
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `/path/to/output/directory/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `fragmentomics cache` subcommand can be used to create a local SIF cache. Please see `fragmentomics cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running fragmentomics with this option when ever possible.
> 
> ***Example:*** `--sif-cache /data/OpenOmics/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for local processes.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for local process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`


---  
  `--tmp-dir TMP_DIR`   
> **Path for writing temporary files.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /data/scratch/$USER/`

---  
  `--overwrite-pipeline-template`   
> **Overwrite pipeline template in output directory.**  
> *type: boolean flag*
> 
> Overwrite pipeline template files in an existing output directory. When this option is provided, the pipeline replaces the existing `workflow/`, `resources/`, and `config/` directories in the output directory with fresh copies from the current pipeline installation. This is mainly useful for developers who want to re-run or test an existing output directory after updating the pipeline template or workflow. 
> 
> ***Example:*** `--overwrite-pipeline-template`

### 2.4 Miscellaneous options  

Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example

### 3.1 Running the pipeline with input BAM files

```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline with
# Input BAM files
./fragmentomics run --input .tests/*.bam \
    --sif-cache /data/OpenOmics/SIFs \
    --output /data/$USER/output \
    --mode slurm \
    --dry-run

# Step 2B.) Run the fragmentomics pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./fragmentomics run --input .tests/*.bam \
    --sif-cache /data/OpenOmics/SIFs \
    --output /data/$USER/output \
    --mode slurm
```

### 3.2 Running the pipeline with input FastQ files

```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline with
# Input paired-end FastQ files, please
# note single-end data is not supported!
./fragmentomics run \
    --input .tests/*.fastq.gz \
    --output /data/$USER/output \
    --genome hg38 \
    --sif-cache /data/OpenOmics/SIFs \
    --mode slurm \
    --dry-run

# Step 2B.) Run the fragmentomics pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./fragmentomics run \
    --input .tests/*.fastq.gz \
    --output /data/$USER/output \
    --genome hg38 \
    --sif-cache /data/OpenOmics/SIFs \
    --mode slurm
```

## 4. Results

The `--output` directory holds one subdirectory per analysis, plus project-level QC. The full layout is documented in [Pipeline overview §4](../pipeline/overview.md#4-output-directory-layout); the files to look at first are:

  `multiqc/multiqc_report.html`
> **Aggregate QC across all samples.** Alignment rates, duplicate rates, read quality, insert-size distributions and per-contig read counts, plus a **Fragmentomics** section summarizing the analysis results themselves — fragment lengths, end motifs, MDS, DELFI and the TSS profiles — for every sample side by side. Written for either input type. See [Quality control](../pipeline/quality-control.md).

---
  `coverage/coverage_summary.xlsx`
> **Merged coverage workbook.** A per-sample summary sheet and the full interval × sample coverage matrix.

---
  `frag_length_bins/{sample}_frag_bin1.png`
> **Fragment-length histogram.** The fastest per-sample sanity check — look for the ~167 bp mononucleosome peak.

---
  `qc/{sample}.contig_validation.txt`
> **Contig validation record** *(BAM input only)*. Which contigs were kept, which were dropped, and how many reads were affected. See [Reference contig filter](../pipeline/contig-filter.md).

Each analysis output is documented on its own page under [Analyses](../analyses/index.md).