# Python standard library
from os.path import join
import os
from textwrap import dedent


# configuration
split_interval                  = int(config['options']['split_interval'])   
max_fragment_len                = int(config['options']['fragment_maximum'])
min_fragment_len                = int(config['options']['fragment_minimum'])
right_flank                     = int(config['options']['right_tss_flank'])
left_flank                      = int(config['options']['left_tss_flank'])
bin_size                        = int(config['options']['bin_size'])
sample_stems                    = config['samples']
# Read filtering thresholds from the frontend's --mapscore/--baseqscore. Both
# are inclusive lower bounds (a read is kept when its score is >= the value)
# and 0 disables the filter. Mapping quality is enforced with samtools on both
# input paths; mean base quality is enforced by fastp on the FastQ path (before
# alignment, on the adapter-trimmed read) and by samtools on the BAM path.
min_mapping_quality             = int(config['options']['mapscore'])
min_base_quality                = int(config['options']['baseqscore'])

# genome linked artifacts
genome_files                    = config["references"][genome]
chrom_sizes                     = genome_files["chrom_sizes"]
ref2bit                         = genome_files["ref2bit"]
intervals                       = genome_files["intervals"]
tss                             = genome_files['tss']
tss_interval                    = genome_files['tss_interval']
gap                             = genome_files.get("gap", None)
blacklist                       = genome_files.get("blacklist", None)
# Picard-style sequence dictionary, used by the bam input path to check that
# staged alignments were made against the selected genome build.
sequence_dict                   = genome_files.get("dict", None)
# FastQ alignment references (only used by the illumina_fastq input path)
reference_fa                    = genome_files.get("reference_fa", None)
bwamem2_index                   = genome_files.get("bwamem2_index", None)

# Resolved input type recorded by the frontend (src/run.py). One input type is
# enforced per run, so the workflow either aligns FastQs (illumina_fastq) or
# stages ready-made BAMs (bam); both converge on bams/{sample}.sorted.bam.
input_type                      = config["project"]["input_type"]

# directories
data_dir                         = config["project"]["datapath"]
all_input_files                  = config['options']['input']
output_dir                       = config['options']['output']
inputs_dir                       = join(output_dir, 'inputs')
bin_dir                          = join(output_dir, 'workflow', 'scripts')
tmpdir                           = config['options']['tmp_dir']
bam_dir                          = join(output_dir, 'bams')
# Landing area for the bam input path. stage_bams sorts user BAMs here, then
# filter_reference_contigs subsets them into bam_dir. Kept as a sibling of
# bam_dir rather than a subdirectory because Snakemake wildcards match across
# path separators, so bams/staged/x would also satisfy bams/{sid}.
staged_bam_dir                   = join(output_dir, 'staged_bams')
bed_dir                          = join(output_dir, 'beds')
coverage_dir                     = join(output_dir, 'coverage')
fragment_length_dir              = join(output_dir, 'frag_length_bins')
fragment_length_int_dir          = join(output_dir, 'frag_length_intervals')
end_motifs_dir                   = join(output_dir, 'end_motifs')
interval_end_motifs_dir          = join(output_dir, 'interval_end_motifs')
mds_dir                          = join(output_dir, 'mds')
delfi_dir                        = join(output_dir, 'delfi')
wps_dir                          = join(output_dir, 'wps')
adjust_wps_dir                   = join(output_dir, 'adjust_wps')
cleavage_profile_dir             = join(output_dir, 'cleavage_profile')
qc_dir                           = join(output_dir, 'qc')
multiqc_dir                      = join(output_dir, 'multiqc')

# project-level (non scatter-per-sample) outputs
coverage_xlsx                    = join(coverage_dir, 'coverage_summary.xlsx')
multiqc_report                   = join(multiqc_dir, 'multiqc_report.html')

# default resources
default_threads                  = cluster['__default__']['threads']


# Input staging. Exactly one input type runs per pipeline invocation (the
# frontend, src/run.py, auto-detects and enforces a single type). Both branches
# below produce the same target, bams/{sample}.sorted.bam, which every
# downstream finaletoolkit rule consumes:
#   • illumina_fastq -> align_fastq: bwa-mem2 alignment + proper-pair filtering
#     + duplicate marking (full pipeline work).
#   • bam            -> stage_bams:  sort/index ready-made alignments (full
#     pipeline work minus alignment).

if input_type == "illumina_fastq":
    rule align_fastq:
        """
        Adapter-trim and quality-filter paired-end Illumina FastQ reads, align
        them to the reference genome with bwa-mem2, then produce an
        analysis-ready coordinate-sorted BAM: proper-pair/primary/
        mapping-quality filter -> fixmate -> orphan re-check -> sort -> markdup.

        Duplicates are marked but deliberately *not* removed. markdup is the
        last stage of the pipe, downstream of every samtools view filter, so
        the 0x400 flags it sets cannot be acted on by those filters; the
        analysis BAM keeps every properly-paired primary read and simply
        records which ones are duplicates. Each downstream rule is then free to
        include or exclude them.

        Filtering runs *before* fixmate so that orphaned mates can be cleaned
        up. Dropping reads on flags or mapping quality can remove one mate of a
        pair while keeping the other, and the survivor would still carry 0x2
        plus mate coordinates pointing at a read that is no longer in the file.
        Because the filter runs while the stream is still name-collated,
        fixmate sees the survivor as a singleton and clears its paired flags
        (0x1/0x2) and mate fields; the second, flag-only `view -f 2` then drops
        those de-paired records. fixmate -r additionally removes unmapped and
        secondary leftovers. Without this, downstream steps would infer
        fragment coordinates from a mate that does not exist.

        Only one sort is performed. samtools fixmate needs name-collated input,
        which is exactly what bwa-mem2 emits for paired FastQ input (both mates
        of a pair are adjacent), so the query-name sort that would otherwise
        precede fixmate is redundant, and the filter and orphan re-check both
        run inside that same name-collated stretch of the pipe. The single
        coordinate sort that follows is the one markdup and every downstream
        rule require.

        Read filtering is driven by the frontend's --baseqscore and --mapscore
        thresholds, both inclusive lower bounds. fastp applies the base-quality
        threshold to the raw reads before they are aligned, discarding any pair
        whose mean base quality falls below it, and the samtools view stage
        applies the mapping-quality threshold to the alignments. Passing 0 for
        either disables that filter.

        fastp also trims adapter read-through before aligning. For paired-end
        input it locates the adapter by overlap analysis - the two mates of a
        short fragment overlap, which reveals where the insert ends and the
        adapter begins - and --detect_adapter_for_pe additionally auto-detects
        the adapter sequence itself for pairs that do not overlap. Trimming runs
        before the quality filter, so --baseqscore is evaluated on the trimmed
        read rather than on adapter bases that are about to be removed, and
        fastp's default --length_required of 15 drops any pair whose read trims
        below 15 bp. --unqualified-percent-limit is pinned to 100 so mean base
        quality remains the only *quality* criterion fastp applies, which keeps
        the meaning of --baseqscore identical on both input paths; the BAM path
        enforces the same threshold with samtools, though it cannot trim.

        FastQC runs on both sides of the alignment - on the raw mates before,
        and on the analysis BAM after - and markdup writes its duplicate
        report, as does fastp for the reads it filtered. All of them land in
        qc/, which is the directory the multiqc rule scans, so they are picked
        up by the aggregate report.

        bwa-mem2 ships one binary per SIMD instruction set rather than a single
        portable executable, so the binary to run is resolved at runtime by
        python_cpu_arch.py. Detection runs in the shell block rather than at
        DAG-build time because the submitting host and the compute node need
        not share a CPU generation, and a binary built for absent instructions
        dies with SIGILL.

        @Input:
            Paired-end FastQ mates (scatter-per-sample), staged by the
            frontend into the output directory's inputs/ folder.
        @Output:
            Coordinate-sorted, indexed analysis BAM with duplicates flagged,
            the markdup duplicate report, and FastQC reports for the raw mates
            and for the analysis BAM.
        """
        input:
            r1                  = join(inputs_dir, "{sid}.R1.fastq.gz"),
            r2                  = join(inputs_dir, "{sid}.R2.fastq.gz"),
        output:
            bam                 = join(bam_dir, "{sid}.sorted.bam"),
            bai                 = join(bam_dir, "{sid}.sorted.bam.bai"),
            markdup             = join(qc_dir, "{sid}.markdup.stats.txt"),
            # FastQC derives these names from its input filenames, so they are
            # not free-form: {sid}.R1.fastq.gz -> {sid}.R1_fastqc.zip and
            # {sid}.sorted.bam -> {sid}.sorted_fastqc.zip. The matching .html
            # reports are written alongside them.
            fastqc_r1           = join(qc_dir, "{sid}.R1_fastqc.zip"),
            fastqc_r2           = join(qc_dir, "{sid}.R2_fastqc.zip"),
            fastqc_bam          = join(qc_dir, "{sid}.sorted_fastqc.zip"),
            fastp_json          = join(qc_dir, "{sid}.fastp.json"),
            fastp_html          = join(qc_dir, "{sid}.fastp.html"),
        container:
            config['images']['bwamem2']
        resources:
            partition = allocated("partition", "align_fastq", cluster),
            mem       = allocated("mem",  "align_fastq", cluster),
            time      = allocated("time", "align_fastq", cluster),
            gres      = allocated("gres", "align_fastq", cluster),
        threads:
            int(allocated("threads", "align_fastq", cluster))
        params:
            rname               = "align_fastq",
            sid                 = "{sid}",
            index               = bwamem2_index,
            arch_script         = join(bin_dir, 'python_cpu_arch.py'),
            # 0x2 (proper pair) kept; 2828 excludes unmapped, mate-unmapped,
            # secondary, qcfail and supplementary reads. This is the usual 3852
            # minus 0x400 (duplicate): duplicate marking runs after this filter
            # precisely so that marked duplicates are not dropped here.
            # keep_flag is used twice: once with drop_flag/map_quality for the
            # main filter, and once on its own after fixmate to drop mates that
            # fixmate de-paired because their partner did not survive.
            keep_flag           = "2",
            drop_flag           = "2828",
            # --mapscore, applied by samtools view: -q keeps alignments whose
            # MAPQ is >= this value. --baseqscore, applied by fastp: a pair is
            # kept when its mean base quality, measured after adapter trimming,
            # is >= this value. 0 disables either filter (samtools -q 0 keeps
            # everything, and 0 is fastp's own "no requirement" value for
            # --average_qual).
            map_quality         = min_mapping_quality,
            base_quality        = min_base_quality,
            # fastp ignores anything above 16 worker threads.
            fastp_threads       = min(int(allocated("threads", "align_fastq", cluster)), 16),
            qc_dir              = qc_dir,
            tmpdir              = tmpdir,
        shell:
            dedent("""
            # Resolve the bwa-mem2 build matching this node's instruction set
            # before doing any work, so an unsupported CPU fails immediately.
            bwa_bin=$(python {params.arch_script}) || exit 1

            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf "${{tmp}}"' EXIT
            mkdir -p \"{params.qc_dir}\"

            # Pre-alignment read QC on the raw mates. Two threads because
            # FastQC parallelizes across input files, not within one.
            fastqc \\
                --threads 2 \\
                --dir \"${{tmp}}\" \\
                --outdir {params.qc_dir} \\
                {input.r1} {input.r2}

            # Adapter trimming plus the base-quality filter (--baseqscore):
            # trim adapter read-through, then drop pairs whose mean base quality
            # is below the threshold. Trimming happens first, so --baseqscore is
            # evaluated on the trimmed read rather than on adapter bases that
            # are about to be removed. The trimmed mates are written to the
            # node's temporary directory and consumed by the aligner below, so
            # no preprocessed FastQ is kept after the step.
            fastp \\
                --in1 {input.r1} \\
                --in2 {input.r2} \\
                --out1 ${{tmp}}/{params.sid}.R1.trimmed.fastq.gz \\
                --out2 ${{tmp}}/{params.sid}.R2.trimmed.fastq.gz \\
                --detect_adapter_for_pe \\
                --average_qual {params.base_quality} \\
                --unqualified_percent_limit 100 \\
                --json {output.fastp_json} \\
                --html {output.fastp_html} \\
                --thread {params.fastp_threads}

            # Alignment and BAM finishing, fused into one pipe so no
            # intermediate BAM is written. Stage order matters:
            #   view   filter on flags + mapping quality while the stream is
            #          still name-collated, so orphans stay adjacent
            #   fixmate fill mate coordinates/ISIZE, add the ms tag markdup
            #          needs, and de-pair any mate whose partner was filtered
            #   view   drop those de-paired orphans (flag-only, no -q)
            #   sort   the one coordinate sort, required by markdup
            #   markdup flag duplicates last so nothing can drop them
            \"${{bwa_bin}}\" mem \\
                -t {threads} \\
                -R "@RG\\tID:{params.sid}\\tSM:{params.sid}\\tPL:ILLUMINA\\tLB:{params.sid}" \\
                {params.index} \\
                ${{tmp}}/{params.sid}.R1.trimmed.fastq.gz \\
                ${{tmp}}/{params.sid}.R2.trimmed.fastq.gz \\
            | samtools view -b -h -f {params.keep_flag} -F {params.drop_flag} \\
                -q {params.map_quality} -@ {threads} - \\
            | samtools fixmate -m -r -@ {threads} - - \\
            | samtools view -b -f {params.keep_flag} -@ {threads} - \\
            | samtools sort -@ {threads} -T ${{tmp}}/psort -O bam - \\
            | samtools markdup -@ {threads} -T ${{tmp}}/mkdup \\
                -f {output.markdup} - {output.bam}

            samtools index -@ {threads} {output.bam}

            # Post-alignment QC on the analysis BAM.
            fastqc \\
                --format bam \\
                --dir \"${{tmp}}\" \\
                --outdir {params.qc_dir} \\
                {output.bam}
            """)

else:
    # Where stage_bams writes. When the selected genome ships a sequence
    # dictionary, staging lands in staged_bam_dir and filter_reference_contigs
    # produces the canonical bams/{sid}.sorted.bam. Without a dictionary there
    # is nothing to validate against, so staging writes bam_dir directly and
    # the filter rule is not defined at all.
    stage_dir = staged_bam_dir if sequence_dict else bam_dir

    rule stage_bams:
        """
        Sort, filter and index ready-made alignments into the staging area.

        The frontend's --mapscore and --baseqscore thresholds are enforced here
        with samtools, since this path has no reads to preprocess: mapping
        quality with `view -q`, and mean base quality with the filter
        expression `avg(qual) >= <threshold>`. Both are inclusive lower bounds
        and either is skipped entirely when its threshold is 0.
        @Input:
            User-provided BAM/CRAM/SAM alignments (gather).
        @Output:
            Coordinate-sorted, indexed, filtered BAMs.
        """
        input:
            all_input_files
        output:
            bams = expand(join(stage_dir, "{sample}.sorted.bam"), sample=sample_stems),
            # stage_input_files.py indexes each staged BAM; declaring the
            # indices lets filter_reference_contigs request regions from them.
            bais = expand(join(stage_dir, "{sample}.sorted.bam.bai"), sample=sample_stems),
        container:
            config['images']['finaletoolkit']
        resources:
            partition = allocated("partition", "stage_bams", cluster),
            mem       = allocated("mem",  "stage_bams", cluster),
            time      = allocated("time", "stage_bams", cluster),
            gres      = allocated("gres", "stage_bams", cluster),
        threads:
            int(allocated("threads", "stage_bams", cluster))
        params:
            rname                   = "stage_bams",
            bam_dir                  = stage_dir,
            python_script            = join(bin_dir, 'stage_input_files.py'),
            memory                   = str(allocated("mem",  "stage_bams", cluster)).replace('G' ,''),
            map_quality              = min_mapping_quality,
            base_quality             = min_base_quality,
        shell:
            dedent("""
            python {params.python_script} \\
                --files {input} \\
                --output {params.bam_dir} \\
                --threads {threads} \\
                --memory {params.memory} \\
                --min-mapq {params.map_quality} \\
                --min-baseq {params.base_quality}
            """)


if input_type != "illumina_fastq" and sequence_dict:

    rule filter_reference_contigs:
        """
        Verify that a staged BAM was aligned to the selected reference genome
        and subset it to the contigs the two share.

        Only the @SQ SN (name) and LN (length) fields of the reference
        sequence dictionary are compared. A .dict records the path of the
        FastA it was built from in UR and a checksum in M5; neither appears in
        a typical aligner-produced BAM header, and UR legitimately differs
        between the reference and the input for the same assembly, so
        comparing either would raise false mismatches.

        Contigs present in the BAM but not the reference (alt/decoy/patch
        scaffolds, often several hundred of them) are dropped. A shared contig
        name with a differing length means the BAM was aligned to a different
        build entirely, which subsetting cannot repair, so the rule fails.
        @Input:
            Staged, coordinate-sorted BAM (scatter-per-sample).
        @Output:
            BAM containing only the reference's contigs, plus a report of the
            comparison.
        """
        input:
            bam                 = join(staged_bam_dir, "{sid}.sorted.bam"),
            bai                 = join(staged_bam_dir, "{sid}.sorted.bam.bai"),
        output:
            bam                 = join(bam_dir, "{sid}.sorted.bam"),
            bai                 = join(bam_dir, "{sid}.sorted.bam.bai"),
            report              = join(qc_dir, "{sid}.contig_validation.txt"),
        container:
            config['images']['finaletoolkit']
        resources:
            partition = allocated("partition", "filter_reference_contigs", cluster),
            mem       = allocated("mem",  "filter_reference_contigs", cluster),
            time      = allocated("time", "filter_reference_contigs", cluster),
            gres      = allocated("gres", "filter_reference_contigs", cluster),
        threads:
            int(allocated("threads", "filter_reference_contigs", cluster))
        params:
            rname                   = "filter_reference_contigs",
            python_script           = join(bin_dir, 'filter_reference_contigs.py'),
            # Declared as a param rather than an input, matching every other
            # reference artifact in this workflow. The dictionary lives on the
            # shared reference filesystem and is not produced by any rule, so
            # listing it as an input only makes DAG construction fail wherever
            # that filesystem is absent (e.g. CI dry-runs).
            seq_dict                = sequence_dict,
        shell:
            dedent("""
            python {params.python_script} \\
                --bam {input.bam} \\
                --dict {params.seq_dict} \\
                --output {output.bam} \\
                --report {output.report} \\
                --threads {threads}
            """)


rule bam_to_bed:
    """
    Convert the canonical analysis BAM into a compressed, tabix-indexed BED6
    of aligned intervals, one record per alignment.

    Reads the same canonical BAM every finaletoolkit rule consumes, so the BED
    is produced identically on both input paths and inherits the filtering
    already applied upstream (proper-pair/primary flags and the
    --mapscore/--baseqscore thresholds). Duplicates are flagged but not
    dropped in that BAM, and this rule does not drop them either: the BED is a
    faithful interval representation of the analysis BAM, and consumers that
    want duplicates excluded are free to filter on the read name or re-derive
    them.

    Compression is bgzip at level 9, the highest the deflate format defines, so
    the output is as small as gzip can make it; the cost is paid once here
    rather than on every read of the result. bgzip rather than gzip because
    BGZF is a gzip-compatible container - `zcat`/`gunzip`/`pandas.read_csv`
    read the result exactly as they would a plain .gz - while additionally
    being block-compressed, which is what makes the tabix index below possible.
    A flat gzip stream at the same level is no smaller and cannot be indexed.
    The rule's threads go to bgzip, which is the slow half of the pipe at level
    9 and the only half that scales; bedtools is single-threaded.

    bedtools writes to stdout and the stream is compressed inside the pipe, so
    no uncompressed BED ever lands on disk - worth avoiding, since the
    plain-text form of a whole-genome cfDNA BAM is several times the size of
    the BAM itself.

    tabix then builds the companion .tbi so consumers can pull a single locus
    out of a whole-genome BED without decompressing it, the same random access
    the .bai gives for the BAM. It is valid here because the records come out
    of bedtools in input order and the input is coordinate-sorted, which is
    exactly tabix's requirement; no intermediate sort is needed. -f overwrites
    a stale index left behind by an interrupted run rather than failing on it.

    @Input:
        Coordinate-sorted, indexed analysis BAM (scatter-per-sample).
    @Output:
        BGZF-compressed (level 9) BED of aligned intervals, plus its tabix
        index.
    """
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
        bai                     = join(bam_dir, "{sid}.sorted.bam.bai"),
    output:
        bed                     = join(bed_dir, "{sid}.bed.gz"),
        tbi                     = join(bed_dir, "{sid}.bed.gz.tbi"),
    container:
        config['images']['bwamem2']
    resources:
        partition = allocated("partition", "bam_to_bed", cluster),
        mem       = allocated("mem",  "bam_to_bed", cluster),
        time      = allocated("time", "bam_to_bed", cluster),
        gres      = allocated("gres", "bam_to_bed", cluster),
    threads:
        int(allocated("threads", "bam_to_bed", cluster))
    params:
        rname                   = "bam_to_bed",
    shell:
        dedent("""
        bedtools bamtobed \\
            -i {input.bam} \\
        | bgzip --compress-level 9 --threads {threads} -c > {output.bed}

        tabix -f -p bed {output.bed}
        """)


rule coverage:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = join(coverage_dir, "{sid}_coverage.bed")
    resources:
        partition = allocated("partition", "coverage", cluster),
        mem       = allocated("mem",  "coverage", cluster),
        time      = allocated("time", "coverage", cluster),
        gres      = allocated("gres", "coverage", cluster),
    threads:
        int(allocated("threads", "coverage", cluster))
    container: 
        config['images']['finaletoolkit']
    params:
        rname                   = "coverage",
        map_quality             = min_mapping_quality,
        intervals               = tss_interval,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len
    shell:
        dedent("""
        finaletoolkit \\
            coverage {input.bam} {params.intervals} \\
            -n \\
            --scale-factor 1000000 \\
            -q {params.map_quality} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            --intersect-policy any \\
            -o {output.bed} \\
            -t {threads} \\
            -v
        """)


rule frag_length_bins:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = join(fragment_length_dir, "{sid}_frag_bin" + str(bin_size) + ".tsv"),
        png                     = join(fragment_length_dir, "{sid}_frag_bin" + str(bin_size) + ".png")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "frag_length_bins", cluster),
        mem       = allocated("mem",  "frag_length_bins", cluster),
        time      = allocated("time", "frag_length_bins", cluster),
        gres      = allocated("gres", "frag_length_bins", cluster),
    threads:
        int(allocated("threads", "frag_length_bins", cluster))
    params:
        rname                   = "frag_length_bins",
        map_quality             = min_mapping_quality,
        bin_size                = str(bin_size),
        min_len                 = str(min_fragment_len),
        max_len                 = str(max_fragment_len),
    shell:
        dedent("""
        finaletoolkit \\
            frag-length-bins {input.bam} \\
            -q {params.map_quality} \\
            --bin-size {params.bin_size} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            -p midpoint \\
            -o {output.tsv} \\
            --histogram {output.png} \\
            -v
        """)


rule frag_length_intervals:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = join(fragment_length_int_dir, "{sid}_frag_interval.bed")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "frag_length_intervals", cluster),
        mem       = allocated("mem",  "frag_length_intervals", cluster),
        time      = allocated("time", "frag_length_intervals", cluster),
        gres      = allocated("gres", "frag_length_intervals", cluster),
    threads:
        int(allocated("threads", "frag_length_intervals", cluster))
    params:
        rname                   = "frag_length_intervals",
        map_quality             = min_mapping_quality,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        intervals               = intervals
    shell:
        dedent("""
        finaletoolkit \\
            frag-length-intervals {input.bam} {params.intervals} \\
            -q {params.map_quality} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            --intersect-policy any \\
            -o {output.bed} \\
            -t {threads} \\
            -v
        """)


rule end_motifs:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = join(end_motifs_dir, "{sid}_endmotif.tsv"),
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "end_motifs", cluster),
        mem       = allocated("mem",  "end_motifs", cluster),
        time      = allocated("time", "end_motifs", cluster),
        gres      = allocated("gres", "end_motifs", cluster),
    threads:
        int(allocated("threads", "end_motifs", cluster))
    params:
        rname                   = "end_motifs",
        map_quality             = min_mapping_quality,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        ref2bit                 = ref2bit,
    shell:
        dedent("""
        finaletoolkit \\
            end-motifs {input.bam} {params.ref2bit} \\
            -q {params.map_quality} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            -o {output.tsv} \\
            -t {threads} \\
            -v
        """)


rule interval_end_motifs:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = join(interval_end_motifs_dir, "{sid}_endmotif_interval.tsv"),
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "interval_end_motifs", cluster),
        mem       = allocated("mem",  "interval_end_motifs", cluster),
        time      = allocated("time", "interval_end_motifs", cluster),
        gres      = allocated("gres", "interval_end_motifs", cluster),
    threads:
        int(allocated("threads", "interval_end_motifs", cluster))
    params:
        rname                   = "interval_end_motifs",
        map_quality             = min_mapping_quality,
        ref2bit                 = ref2bit,
        intervals               = intervals,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        tmpdir                  = tmpdir
    shell:
        dedent("""
        if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{params.tmpdir}\")
        trap 'ls -al ${{tmp}}; rm -rf "${{tmp}}"' EXIT

        finaletoolkit \\
            interval-end-motifs {input.bam} {params.ref2bit} {params.intervals} \\
            -q {params.map_quality} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            -o {output.tsv} \\
            -t {threads} \\
            -v
        """)


rule mds:
    input:
        endmotif                = join(end_motifs_dir, "{sid}_endmotif.tsv"),
    output:
        tsv                     = join(mds_dir, "{sid}_mds.tsv"),
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "mds", cluster),
        mem       = allocated("mem",  "mds", cluster),
        time      = allocated("time", "mds", cluster),
        gres      = allocated("gres", "mds", cluster),
    threads:
        int(allocated("threads", "mds", cluster))
    params:
        rname                   = "mds",
        sid                     = "{sid}"
    shell:
        dedent("""
        mds_score=$(finaletoolkit mds {input.endmotif})
        echo "Sample\tMDS_score" > {output.tsv}
        echo "{params.sid}\t${{mds_score}}" >> {output.tsv}
        """)


rule delfi:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = join(delfi_dir, "{sid}_delfi.bed"),
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "delfi", cluster),
        mem       = allocated("mem",  "delfi", cluster),
        time      = allocated("time", "delfi", cluster),
        gres      = allocated("gres", "delfi", cluster),
    threads:
        int(allocated("threads", "delfi", cluster))
    params:
        rname                   = "delfi",
        map_quality             = min_mapping_quality,
        chrom_sizes             = chrom_sizes,
        ref2bit                 = ref2bit,
        intervals               = intervals,
        blacklist_cmd           = f" --blacklist {blacklist}" if blacklist else "",
        gap_cmd                 = f" -g {gap}" if gap else ""
    shell:
        dedent("""
        finaletoolkit delfi {input.bam} {params.chrom_sizes} {params.ref2bit} {params.intervals} \\
            -q {params.map_quality}{params.blacklist_cmd}{params.gap_cmd} \\
            -o {output.bed} \\
            -t {threads} \\
            -v \\
            --no-merge-bins
        """)


rule wps:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        bw                      = join(wps_dir, "{sid}_wps_out_tss.bw")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "wps", cluster),
        mem       = allocated("mem",  "wps", cluster),
        time      = allocated("time", "wps", cluster),
        gres      = allocated("gres", "wps", cluster),
    threads:
        int(allocated("threads", "wps", cluster))
    params:
        rname                   = "wps",
        map_quality             = min_mapping_quality,
        intervals               = split_interval,
        tss                     = tss
    shell:
        """
        finaletoolkit \\
            wps {input.bam} {params.tss} \\
            -i {params.intervals} \\
            -W 120 \\
            --min-length 120 \\
            --max-length 180 \\
            -q {params.map_quality} \\
            -o {output.bw} \\
            -t {threads} \\
            -v
        """


rule adjust_wps:
    input:
        wps_bw                  = join(wps_dir, "{sid}_wps_out_tss.bw")
    output:
        bw                      = join(adjust_wps_dir, "{sid}_wps_out_tss_adjusted.bw")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "adjust_wps", cluster),
        mem       = allocated("mem",  "adjust_wps", cluster),
        time      = allocated("time", "adjust_wps", cluster),
        gres      = allocated("gres", "adjust_wps", cluster),
    threads:
        int(allocated("threads", "adjust_wps", cluster))
    params:
        rname                   = "adjust_wps",
        intervals               = split_interval,
        tss_interval            = tss_interval,
        chrom_sizes             = chrom_sizes
    shell:
        dedent("""
        finaletoolkit \\
            adjust-wps {input.wps_bw} {params.tss_interval} {params.chrom_sizes} \\
            -o {output.bw} \\
            -i {params.intervals} \\
            -m 200 \\
            --subtract-edges \\
            --no-savgol \\
            -v
        """)


rule cleavage_profile:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
    output:
        bw                      = join(cleavage_profile_dir, "{sid}_cleavage_profile_tss.bw")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "cleavage_profile", cluster),
        mem       = allocated("mem",  "cleavage_profile", cluster),
        time      = allocated("time", "cleavage_profile", cluster),
        gres      = allocated("gres", "cleavage_profile", cluster),
    threads:
        int(allocated("threads", "cleavage_profile", cluster))
    params:
        rname                   = "cleavage_profile",
        map_quality             = min_mapping_quality,
        l                       = left_flank,
        r                       = right_flank,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        tss                     = tss,
        chrom_sizes             = chrom_sizes
    shell:
        dedent("""
        finaletoolkit \\
            cleavage-profile {input.bam} {params.tss} {params.chrom_sizes} \\
            -o {output.bw} \\
            --pad-left {params.l} \\
            --pad-right {params.r} \\
            -q {params.map_quality} \\
            --min-length {params.min_len} \\
            --max-length {params.max_len} \\
            -t {threads} \\
            -v
        """)


rule agg_wps:
    input:
        bw                      = join(wps_dir, "{sid}_wps_out_tss.bw")
    output:
        wig                     = join(wps_dir, "{sid}_wps_out_tss_aggr.wig")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "agg_wps", cluster),
        mem       = allocated("mem",  "agg_wps", cluster),
        time      = allocated("time", "agg_wps", cluster),
        gres      = allocated("gres", "agg_wps", cluster),
    threads:
        int(allocated("threads", "agg_wps", cluster))
    params:
        rname                   = "agg_wps",
        tss_interval            = tss_interval
    shell:
        dedent("""
        finaletoolkit \\
            agg-bw {input.bw} {params.tss_interval} \\
            -o {output.wig} \\
            --mean \\
            -v
        """)

rule agg_adjust_wps:
    input:
        bw                      = join(adjust_wps_dir, "{sid}_wps_out_tss_adjusted.bw")
    output:
        wig                     = join(adjust_wps_dir, "{sid}_wps_out_tss_adj_aggr.wig")
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "agg_adjust_wps", cluster),
        mem       = allocated("mem",  "agg_adjust_wps", cluster),
        time      = allocated("time", "agg_adjust_wps", cluster),
        gres      = allocated("gres", "agg_adjust_wps", cluster),
    threads:
        int(allocated("threads", "agg_adjust_wps", cluster))
    params:
        rname                   = "agg_adjust_wps",
        tss_interval            = tss_interval
    shell:
        dedent("""
        finaletoolkit \\
            agg-bw {input.bw} {params.tss_interval} \\
            -o {output.wig} \\
            --mean \\
            -v
        """)


rule bam_stats:
    """
    Collect alignment QC metrics from the analysis BAM with samtools. The
    three reports (stats, flagstat, idxstats) are all natively parsed by
    MultiQC, which gives the aggregate report something to summarize
    regardless of which input path produced the BAM.
    @Input:
        Coordinate-sorted, indexed analysis BAM (scatter-per-sample).
    @Output:
        samtools stats, flagstat and idxstats reports.
    """
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
        bai                     = join(bam_dir, "{sid}.sorted.bam.bai"),
    output:
        stats                   = join(qc_dir, "{sid}.samtools.stats.txt"),
        flagstat                = join(qc_dir, "{sid}.flagstat.txt"),
        idxstats                = join(qc_dir, "{sid}.idxstats.txt"),
    container:
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "bam_stats", cluster),
        mem       = allocated("mem",  "bam_stats", cluster),
        time      = allocated("time", "bam_stats", cluster),
        gres      = allocated("gres", "bam_stats", cluster),
    threads:
        int(allocated("threads", "bam_stats", cluster))
    params:
        rname                   = "bam_stats",
    shell:
        dedent("""
        samtools stats -@ {threads} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        samtools idxstats {input.bam} > {output.idxstats}
        """)


rule merge_coverage_excel:
    """
    Merge every per-sample finaletoolkit coverage BED into one Excel
    workbook: a `summary` sheet of per-sample coverage statistics and a
    `coverage` sheet holding the merged interval x sample matrix. This is a
    project-level (gather) rule, so it waits on all samples.
    @Input:
        Per-sample coverage BEDs from the coverage rule (gather).
    @Output:
        Single Excel workbook summarizing coverage across all samples.
    """
    input:
        beds                    = expand(join(coverage_dir, "{sample}_coverage.bed"), sample=sample_stems),
    output:
        xlsx                    = coverage_xlsx,
    container:
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "merge_coverage_excel", cluster),
        mem       = allocated("mem",  "merge_coverage_excel", cluster),
        time      = allocated("time", "merge_coverage_excel", cluster),
        gres      = allocated("gres", "merge_coverage_excel", cluster),
    threads:
        int(allocated("threads", "merge_coverage_excel", cluster))
    params:
        rname                   = "merge_coverage_excel",
        python_script           = join(bin_dir, 'merge_coverage.py'),
    shell:
        dedent("""
        python {params.python_script} \\
            --beds {input.beds} \\
            --output {output.xlsx}
        """)


# QC artifacts that only the FastQ input path produces: the FastQC reports
# align_fastq writes on either side of the alignment, the samtools markdup
# duplicate report, and the fastp trimming/filtering report. They are guaranteed
# to exist by the time multiqc runs
# (align_fastq also produces the BAM that bam_stats consumes), but they are
# requested explicitly so the aggregate report's dependency on them is visible
# in the DAG. Empty for a BAM-input run, where no such reports exist.
if input_type == "illumina_fastq":
    fastq_qc_reports = (
        expand(join(qc_dir, "{sample}.R1_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.R2_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.sorted_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.markdup.stats.txt"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.fastp.json"), sample=sample_stems)
    )
else:
    fastq_qc_reports = []


rule multiqc:
    """
    Aggregate the per-sample QC reports into a single interactive HTML report.
    Only the qc/ directory is scanned, so MultiQC does not walk the large
    bigwig/bed outputs of the fragmentomics rules. The coverage workbook is
    taken as an input so the report is generated once the project-level
    coverage gather has finished.
    @Input:
        Per-sample samtools stats/flagstat/idxstats reports (gather), plus the
        FastQC and markdup reports on a FastQ run, and the merged coverage
        workbook.
    @Output:
        MultiQC HTML report.
    """
    input:
        stats                   = expand(join(qc_dir, "{sample}.samtools.stats.txt"), sample=sample_stems),
        flagstat                = expand(join(qc_dir, "{sample}.flagstat.txt"), sample=sample_stems),
        idxstats                = expand(join(qc_dir, "{sample}.idxstats.txt"), sample=sample_stems),
        fastq_qc                = fastq_qc_reports,
        xlsx                    = coverage_xlsx,
    output:
        report                  = multiqc_report,
        data                    = directory(join(multiqc_dir, "multiqc_report_data")),
    container:
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "multiqc", cluster),
        mem       = allocated("mem",  "multiqc", cluster),
        time      = allocated("time", "multiqc", cluster),
        gres      = allocated("gres", "multiqc", cluster),
    threads:
        int(allocated("threads", "multiqc", cluster))
    params:
        rname                   = "multiqc",
        qc_dir                  = qc_dir,
        multiqc_dir             = multiqc_dir,
        # MultiQC appends .html itself, so the report basename is passed
        # without its extension.
        report_name             = "multiqc_report",
    shell:
        dedent("""
        multiqc {params.qc_dir} \\
            --outdir {params.multiqc_dir} \\
            --filename {params.report_name} \\
            --force
        """)


rule agg_cleavage_profile:
    input:
        bw                      = join(cleavage_profile_dir, "{sid}_cleavage_profile_tss.bw"),
    output:
        wig                     = join(cleavage_profile_dir, "{sid}_cleavage_profile_aggr.wig"),
    container: 
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "agg_cleavage_profile", cluster),
        mem       = allocated("mem",  "agg_cleavage_profile", cluster),
        time      = allocated("time", "agg_cleavage_profile", cluster),
        gres      = allocated("gres", "agg_cleavage_profile", cluster),
    threads:
        int(allocated("threads", "agg_cleavage_profile", cluster))
    params:
        rname                   = "agg_cleavage_profile",
        tss_interval            = tss_interval     
    shell:
        dedent("""
        finaletoolkit \\
            agg-bw {input.bw} {params.tss_interval} \\
            -o {output.wig} \\
            --mean \\
            -v
        """)
