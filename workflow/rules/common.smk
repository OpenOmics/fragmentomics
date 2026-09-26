# Python standard library
from os.path import join
import os
from textwrap import dedent


def local_tmp_dir(configured):
    """Resolve a --tmp-dir that names a SLURM variable for a run that has no
    cluster job to name it after. See the tmpdir assignment below for why this
    is needed. The configured path is returned unchanged when it references no
    SLURM variable, since then there is nothing to resolve.
    @param configured <str>:
        The --tmp-dir recorded in config.json, e.g. /lscratch/$SLURM_JOBID/
    @return <str>:
        A temporary directory that this machine can actually write to
    """
    slurm_vars = (
        '${SLURM_JOBID}', '$SLURM_JOBID', '${SLURM_JOB_ID}', '$SLURM_JOB_ID'
    )
    if not any(var in configured for var in slurm_vars):
        return configured

    # A local run started from an interactive allocation (sinteractive) does
    # have a job id, and reusing it keeps node-local disk in play.
    job_id = os.environ.get('SLURM_JOBID', os.environ.get('SLURM_JOB_ID', ''))
    if job_id:
        resolved = configured
        for var in slurm_vars:
            resolved = resolved.replace(var, job_id)
        # /lscratch/$SLURM_JOBID only exists if the allocation asked for the
        # lscratch gres, so the substitution is usable only when the directory
        # it names is really there.
        if os.path.isdir(resolved):
            return resolved

    return os.environ.get('TMPDIR', '/tmp')


# configuration
split_interval                  = int(config['options']['split_interval'])   
max_fragment_len                = int(config['options']['fragment_maximum'])
min_fragment_len                = int(config['options']['fragment_minimum'])
right_flank                     = int(config['options']['right_tss_flank'])
left_flank                      = int(config['options']['left_tss_flank'])
bin_size                        = int(config['options']['bin_size'])
# Width of the fixed-size genomic windows make_intervals tiles the reference
# into, as a whole number of bases (--interval normalizes the user's SI-prefixed
# base unit), plus the short label derived from it that names the BED.
interval_width                  = int(config['options']['interval'])
interval_name                   = interval_label(interval_width)
sample_stems                    = config['samples']
# Read filtering thresholds from the frontend's --mapscore/--baseqscore. Both
# are inclusive lower bounds (a read is kept when its score is >= the value)
# and 0 disables the filter. Mapping quality is enforced with samtools on both
# input paths; mean base quality is enforced by fastp on the FastQ path (before
# alignment, on the adapter-trimmed read) and by samtools on the BAM path.
min_mapping_quality             = int(config['options']['mapscore'])
min_base_quality                = int(config['options']['baseqscore'])

# genome linked artifacts
#
# Every one of these is read with .get() because a build is not obliged to define
# all of them: the bundled builds do, but --genome also takes a config file
# describing references of the user's own, and such a build may only have some of
# them. Which analyses a run performs is decided by which of these keys the
# selected build supplies (see the gating at the bottom of this file and in the
# Snakefile), so a key that is absent leaves the rules that read it defined but
# never requested, rather than failing the workflow as it is being built.
genome_files                    = config["references"][genome]
chrom_sizes                     = genome_files.get("chrom_sizes", None)
ref2bit                         = genome_files.get("ref2bit", None)
# Note: the genomic interval BED is not a bundled reference. It is generated per
# run by make_intervals, so its path is derived from the output directory below
# rather than read from config/genome.json.
tss                             = genome_files.get('tss', None)
tss_interval                    = genome_files.get('tss_interval', None)
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

# Executor the frontend recorded for this run: 'slurm' submits every rule as
# its own cluster job, 'local' runs them all on the current machine.
run_mode                        = config['options']['mode']

# directories
data_dir                         = config["project"]["datapath"]
all_input_files                  = config['options']['input']
output_dir                       = config['options']['output']
inputs_dir                       = join(output_dir, 'inputs')
bin_dir                          = join(output_dir, 'workflow', 'scripts')
# Scratch space handed to every rule that stages intermediate files, as its
# params.tmpdir. --tmp-dir defaults to /lscratch/$SLURM_JOBID, which is right
# only under the slurm executor: each rule runs in its own allocation, so the
# job id has to stay unexpanded here and be filled in by the shell of whichever
# job ends up running the rule. A local run has no allocation to name, and
# because Snakemake runs shell blocks under `set -u` the unexpanded reference
# aborts the rule outright with "SLURM_JOBID: unbound variable" before it does
# any work. So resolve the path up front for anything but a cluster job.
tmpdir                           = config['options']['tmp_dir'] if run_mode == 'slurm' \
                                   else local_tmp_dir(config['options']['tmp_dir'])
bam_dir                          = join(output_dir, 'bams')
# Landing area for the bam input path. stage_bams sorts user BAMs here, then
# filter_reference_contigs subsets them into bam_dir. Kept as a sibling of
# bam_dir rather than a subdirectory because Snakemake wildcards match across
# path separators, so bams/staged/x would also satisfy bams/{sid}.
staged_bam_dir                   = join(output_dir, 'staged_bams')
bed_dir                          = join(output_dir, 'beds')
intervals_dir                    = join(output_dir, 'intervals')
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
# Fixed-size windows tiling the reference, built once per run by make_intervals
# and shared by every interval-based analysis. The width is in the filename so a
# run at a different --interval writes a new BED rather than silently reusing an
# existing one built at another size.
intervals                        = join(
                                     intervals_dir,
                                     f'{genome}_{interval_name}_intervals.bed'
                                   )

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
        mapping-quality filter -> name sort -> fixmate -> orphan re-check ->
        coordinate sort -> markdup.

        Duplicates are marked but deliberately *not* removed. markdup is the
        last stage of the pipe, downstream of every samtools view filter, so
        the 0x400 flags it sets cannot be acted on by those filters; the
        analysis BAM keeps every properly-paired primary read and simply
        records which ones are duplicates. Each downstream rule is then free to
        include or exclude them.

        Two sorts are performed, and each is required by the stage it feeds.
        samtools fixmate needs query-name-ordered input - it pairs records by
        walking adjacent reads of the same name - so `sort -n` precedes it
        rather than relying on the ordering bwa-mem2 happens to emit. The
        coordinate sort that follows is the one markdup and every downstream
        rule require. Both spill to the node's temporary directory, so this
        rule's lscratch allocation covers two passes over the read data.

        Filtering runs *before* the name sort so that orphaned mates can be
        cleaned up, and so the sort only handles reads that survive. Dropping
        reads on flags or mapping quality can remove one mate of a pair while
        keeping the other, and the survivor would still carry 0x2 plus mate
        coordinates pointing at a read that is no longer in the file. With the
        stream name-ordered, fixmate sees the survivor as a singleton and clears
        its paired flags (0x1/0x2) and mate fields; the second, flag-only
        `view -f 2` then drops those de-paired records. fixmate -r additionally
        removes unmapped and secondary leftovers. Without this, downstream steps
        would infer fragment coordinates from a mate that does not exist.

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
            #   view   filter on flags + mapping quality first, so the name
            #          sort below only handles the reads that survive
            #   sort -n query-name sort, the input ordering fixmate requires
            #   fixmate fill mate coordinates/ISIZE, add the ms tag markdup
            #          needs, and de-pair any mate whose partner was filtered
            #   view   drop those de-paired orphans (flag-only, no -q)
            #   sort   coordinate sort, required by markdup
            #   markdup flag duplicates last so nothing can drop them
            \"${{bwa_bin}}\" mem \\
                -t {threads} \\
                -R "@RG\\tID:{params.sid}\\tSM:{params.sid}\\tPL:ILLUMINA\\tLB:{params.sid}" \\
                {params.index} \\
                ${{tmp}}/{params.sid}.R1.trimmed.fastq.gz \\
                ${{tmp}}/{params.sid}.R2.trimmed.fastq.gz \\
            | samtools view -b -h -f {params.keep_flag} -F {params.drop_flag} \\
                -q {params.map_quality} -@ {threads} - \\
            | samtools sort -n -@ {threads} -T ${{tmp}}/nsort -O bam - \\
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
    # The bam input path is a three-step chain, and every step but the last
    # writes into staged_bam_dir:
    #   stage_bams            sort + quality filter user alignments
    #   remove_orphan_reads   repair the pairs that filtering broke
    #   filter_reference_contigs  subset to the reference's contigs
    # Whichever step runs last produces the canonical bams/{sid}.sorted.bam.
    # When the selected genome ships no sequence dictionary there is nothing to
    # validate against, so filter_reference_contigs is not defined at all and
    # remove_orphan_reads writes bam_dir directly; otherwise it hands a
    # temporary intermediate to the filter rule, which Snakemake deletes once
    # that rule has consumed it (three copies of a whole-genome BAM on disk at
    # once is a lot of space to hold for no reason).
    stage_dir = staged_bam_dir
    if sequence_dict:
        paired_bam  = join(staged_bam_dir, "{sid}.paired.bam")
        paired_out  = temp(paired_bam)
        paired_bai  = temp(paired_bam + '.bai')
    else:
        paired_bam  = join(bam_dir, "{sid}.sorted.bam")
        paired_out  = paired_bam
        paired_bai  = paired_bam + '.bai'

    rule stage_bams:
        """
        Sort, filter and index ready-made alignments into the staging area.

        The frontend's --mapscore and --baseqscore thresholds are enforced here
        with samtools, since this path has no reads to preprocess: mapping
        quality with `view -q`, and mean base quality with the filter
        expression `avg(qual) >= <threshold>`. Both are inclusive lower bounds
        and either is skipped entirely when its threshold is 0.

        Both filters act on individual records, so either can keep one mate of
        a pair and drop the other. Repairing that is remove_orphan_reads' job,
        which is why this rule writes into the staging area rather than
        producing the analysis BAM directly.
        @Input:
            User-provided BAM alignments (gather).
        @Output:
            Coordinate-sorted, indexed, filtered BAMs in the staging area.
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

    rule remove_orphan_reads:
        """
        Drop singletons and orphaned mates from a staged BAM, so every read
        that reaches the analysis BAM has its mate beside it.

        stage_bams filters on mapping quality and mean base quality, and both
        act on individual records: a pair whose R1 passes and whose R2 fails
        loses only R2, and the surviving R1 keeps its 0x2 (properly paired) flag
        plus mate coordinates pointing at a record that is no longer in the
        file. Every fragmentomics feature here is derived from a pair's
        coordinates, so such a read is at best dead weight and at worst a
        fragment inferred from a mate that does not exist. Reads whose mate went
        unmapped are the same problem arriving from the input BAM itself.

        The pipe is the same shape as the FastQ path's, for the same reason
        (see align_fastq):
          view    keep primary paired records with both mates mapped, which is
                  also what makes the pairing invariant well defined - a
                  secondary or supplementary record shares its name with a
                  primary one, so "every name appears twice" cannot hold while
                  they are present
          sort -n query-name order, the input ordering fixmate requires
          fixmate recompute mate coordinates/ISIZE and de-pair any read whose
                  mate is missing, clearing its 0x1/0x2 flags and mate fields
          view    drop those de-paired records, flags only
          sort    back to coordinate order for every downstream rule
        fixmate is run without -m: that flag exists to add the mate-score tag
        samtools markdup needs, and this path does not mark duplicates - the
        input BAM's own duplicate flags are preserved untouched.

        Single-end input is passed through unchanged rather than emptied. A BAM
        with no paired records at all has no orphans to remove but would not
        survive `view -f 1`, so the read count of the result is checked and the
        staged BAM is copied through with a warning if the pipe removed
        everything. Both counts are read from the BAM indices rather than from a
        pass over the reads, so the check is free in the normal case.

        Verified against samtools 1.13 on a staged BAM carrying 710 orphans
        (555,122 records): 554,412 records survive, flagstat reports 0
        singletons, and every surviving read name appears exactly twice.
        @Input:
            Staged, coordinate-sorted BAM (scatter-per-sample).
        @Output:
            Coordinate-sorted, indexed BAM whose every read is one of a
            complete, mapped pair.
        """
        input:
            bam                 = join(staged_bam_dir, "{sid}.sorted.bam"),
            bai                 = join(staged_bam_dir, "{sid}.sorted.bam.bai"),
        output:
            bam                 = paired_out,
            bai                 = paired_bai,
        container:
            config['images']['finaletoolkit']
        resources:
            partition = allocated("partition", "remove_orphan_reads", cluster),
            mem       = allocated("mem",  "remove_orphan_reads", cluster),
            time      = allocated("time", "remove_orphan_reads", cluster),
            gres      = allocated("gres", "remove_orphan_reads", cluster),
        threads:
            int(allocated("threads", "remove_orphan_reads", cluster))
        params:
            rname               = "remove_orphan_reads",
            # 0x1 (paired) required; 2316 excludes unmapped (0x4), mate-unmapped
            # (0x8), secondary (0x100) and supplementary (0x800) records. The
            # mask deliberately omits 0x2 (proper pair), so discordant pairs are
            # kept as long as both mates are present, and 0x400 (duplicate), so
            # the input BAM's duplicate marking survives this step. keep_flag is
            # used twice: once with drop_flag for the main filter, and once on
            # its own after fixmate to drop the mates fixmate de-paired.
            keep_flag           = "1",
            drop_flag           = "2316",
            tmpdir              = tmpdir,
        shell:
            dedent("""
            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf "${{tmp}}"' EXIT

            samtools view -b -f {params.keep_flag} -F {params.drop_flag} \\
                -@ {threads} {input.bam} \\
            | samtools sort -n -@ {threads} -T ${{tmp}}/nsort -O bam - \\
            | samtools fixmate -@ {threads} - - \\
            | samtools view -b -f {params.keep_flag} -@ {threads} - \\
            | samtools sort -@ {threads} -T ${{tmp}}/psort -O bam -o {output.bam} -

            samtools index -@ {threads} {output.bam}

            # Nothing left implies there was nothing paired to begin with:
            # single-end input, which has no orphans but cannot survive -f 1.
            # Both totals come from the .bai rather than from the reads.
            kept=$(samtools idxstats {output.bam} \\
                | awk '{{total += $3 + $4}} END {{print total + 0}}')
            if [ "${{kept}}" -eq 0 ]; then
                staged=$(samtools idxstats {input.bam} \\
                    | awk '{{total += $3 + $4}} END {{print total + 0}}')
                if [ "${{staged}}" -gt 0 ]; then
                    echo "WARNING: no paired reads in {input.bam}, so there are" \\
                         "no orphans to remove. Copying it through unchanged;" \\
                         "note that this pipeline's fragment-level analyses" \\
                         "expect paired-end data." >&2
                    cp {input.bam} {output.bam}
                    samtools index -@ {threads} {output.bam}
                fi
            fi
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

        Reads whose mate lies on a dropped contig are excluded along with the
        contig, so subsetting cannot undo the pairing guarantee
        remove_orphan_reads established upstream of it.
        @Input:
            Staged, coordinate-sorted, pairing-cleaned BAM
            (scatter-per-sample).
        @Output:
            BAM containing only the reference's contigs, plus a report of the
            comparison.
        """
        input:
            bam                 = join(staged_bam_dir, "{sid}.paired.bam"),
            bai                 = join(staged_bam_dir, "{sid}.paired.bam.bai"),
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


if input_type != "illumina_fastq":

    # Read QC for the bam input path. The FastQ path gets both of these from
    # align_fastq as a side effect of the work it is already doing - FastQC on
    # the raw mates and on the BAM it just built, fastp on the reads it trims
    # before alignment - so on a BAM run the aggregate report would otherwise
    # have no FastQC or fastp section at all. Neither rule changes the
    # analysis BAM; both only report on it.

    rule fastqc_bam:
        """
        Read QC on the analysis BAM, the BAM path's counterpart to the
        post-alignment FastQC that align_fastq runs.

        A BAM run has no raw reads to inspect, but the analysis BAM itself is
        a format FastQC reads natively (--format bam), so the post-alignment
        half of the FastQ path's read QC can be produced for it: per-base and
        per-sequence quality, GC content, read-length distribution, adapter
        and overrepresented-sequence content.

        Report names are derived by FastQC from its input filename rather than
        chosen here: {sid}.sorted.bam -> {sid}.sorted_fastqc.zip, with the
        HTML copy alongside. That is the same name align_fastq's
        post-alignment report carries, so MultiQC labels the sample
        identically on either input path and its FastQC section means the same
        thing in both.

        One caveat applies here exactly as it does to align_fastq's
        post-alignment report: FastQC estimates duplication and
        overrepresented sequences from the first 100,000 reads it sees, and in
        a coordinate-sorted BAM those all come from the start of the first
        contig instead of from across the library. Read the quality and
        content modules from this report, and take the duplicate rate from
        `samtools flagstat` (see bam_stats), which counts flags over the whole
        file.
        @Input:
            Coordinate-sorted analysis BAM (scatter-per-sample).
        @Output:
            FastQC report for the analysis BAM.
        """
        input:
            bam                 = join(bam_dir, "{sid}.sorted.bam"),
        output:
            zip                 = join(qc_dir, "{sid}.sorted_fastqc.zip"),
            html                = join(qc_dir, "{sid}.sorted_fastqc.html"),
        container:
            config['images']['bwamem2']
        resources:
            partition = allocated("partition", "fastqc_bam", cluster),
            mem       = allocated("mem",  "fastqc_bam", cluster),
            time      = allocated("time", "fastqc_bam", cluster),
            gres      = allocated("gres", "fastqc_bam", cluster),
        threads:
            int(allocated("threads", "fastqc_bam", cluster))
        params:
            rname               = "fastqc_bam",
            qc_dir              = qc_dir,
            tmpdir              = tmpdir,
        shell:
            dedent("""
            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf "${{tmp}}"' EXIT
            mkdir -p {params.qc_dir}

            # No --threads here: FastQC parallelizes across input files, not
            # within one, and this rule hands it a single BAM. --dir keeps its
            # temporary files on the node rather than in the output directory.
            fastqc \\
                --format bam \\
                --dir \"${{tmp}}\" \\
                --outdir {params.qc_dir} \\
                {input.bam}
            """)


    rule fastp_bam:
        """
        fastp report for the analysis BAM, produced by converting it back to
        FastQ on the fly.

        On the FastQ path fastp is a processing step - it trims adapter
        read-through and drops pairs below the --baseqscore mean-quality
        threshold before alignment - and its report is a by-product of that
        work. A BAM run cannot do any of it: the reads are already aligned,
        and --mapscore/--baseqscore were already applied to them with samtools
        during staging. What it can still have is the report, which is what
        this rule produces: the analysis BAM is streamed back to FastQ and
        fastp is pointed at the result with every transformation switched off,
        so it acts purely as a reporter.

        Nothing is written back out. fastp given neither --out1 nor --out2
        discards the reads it read and writes only its reports, and the FastQ
        pair it read lives in the node's temporary directory, so the step
        leaves behind two report files and nothing else.

        Trimming and both filters are disabled deliberately, so that every
        number in the report describes the reads that are actually in the
        analysis BAM. Left enabled, fastp would report adapter bases removed
        and reads dropped that no downstream step will lose, which reads as
        data loss this run never had - and worse, the same MultiQC section
        would mean the opposite thing on a FastQ run, where those removals are
        real. The cost is that the before/after halves of the report are
        identical by construction; what is worth reading is everything that
        describes the reads themselves - the per-base quality and content
        curves, Q20/Q30 rates, GC content, duplication rate and the
        insert-size distribution fastp estimates from the overlap between
        mates. Adapter content is not lost either, it is in the FastQC report
        for the same BAM (fastqc_bam).

        Conversion is collate -> fastq rather than sort -> fastq. fastp reads
        the two mate files in lockstep and needs the nth record of one to be
        the mate of the nth record of the other, which is what collating by
        name gives it, and collating is cheaper than a full coordinate sort.
        The BAM stream between the two stages is uncompressed and never
        touches disk, -n leaves the read names untouched so both mates keep
        the same name as they would in a raw FastQ pair, and secondary and
        supplementary alignments are excluded by `samtools fastq`'s own
        default filter.

        Paired or single-end mode is read off the sample's flagstat report
        rather than by counting the BAM a second time: a staged BAM is
        whatever the user provided, and pointing fastp's paired-end mode at an
        unpaired BAM would fail the step. Reads whose mate did not survive
        filtering (samtools calls these singletons) are left out of the
        report, since they cannot be reported as pairs; there are normally
        very few of them.
        @Input:
            Coordinate-sorted analysis BAM and the flagstat report bam_stats
            wrote for it (scatter-per-sample).
        @Output:
            fastp JSON report, which is the one MultiQC parses, plus its HTML
            copy.
        """
        input:
            bam                 = join(bam_dir, "{sid}.sorted.bam"),
            flagstat            = join(qc_dir, "{sid}.flagstat.txt"),
        output:
            json                = join(qc_dir, "{sid}.fastp.json"),
            html                = join(qc_dir, "{sid}.fastp.html"),
        container:
            config['images']['bwamem2']
        resources:
            partition = allocated("partition", "fastp_bam", cluster),
            mem       = allocated("mem",  "fastp_bam", cluster),
            time      = allocated("time", "fastp_bam", cluster),
            gres      = allocated("gres", "fastp_bam", cluster),
        threads:
            int(allocated("threads", "fastp_bam", cluster))
        params:
            rname               = "fastp_bam",
            sid                 = "{sid}",
            qc_dir              = qc_dir,
            tmpdir              = tmpdir,
            # Switches off everything fastp can do to a read, leaving only
            # the measuring: adapter trimming, the mean/per-base quality
            # filter, the post-trim length filter, and polyG tail trimming
            # (which fastp turns on by itself for NextSeq/NovaSeq data).
            qc_only             = "--disable_adapter_trimming "
                                  "--disable_quality_filtering "
                                  "--disable_length_filtering "
                                  "--disable_trim_poly_g",
            # fastp ignores anything above 16 worker threads.
            fastp_threads       = min(int(allocated("threads", "fastp_bam", cluster)), 16),
        shell:
            dedent("""
            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf "${{tmp}}"' EXIT
            mkdir -p {params.qc_dir}

            # Read layout of the staged BAM, taken from the report bam_stats
            # already wrote for it. Empty (no such line) is treated as
            # single-end below.
            paired=$(awk '/paired in sequencing/ {{print $1; exit}}' {input.flagstat})

            # BAM -> FastQ. collate brings both mates of a pair together,
            # which is what samtools fastq needs to write matched mate files;
            # -u keeps the stream between the two uncompressed. Reads with
            # neither mate flag go to the unpaired file, reads whose mate was
            # filtered away to the singleton file.
            samtools collate -u -O -@ {threads} {input.bam} ${{tmp}}/collate \\
                | samtools fastq -n -@ {threads} \\
                    -1 ${{tmp}}/{params.sid}.R1.fastq.gz \\
                    -2 ${{tmp}}/{params.sid}.R2.fastq.gz \\
                    -0 ${{tmp}}/{params.sid}.unpaired.fastq.gz \\
                    -s ${{tmp}}/{params.sid}.singleton.fastq.gz \\
                    -

            if [ "${{paired:-0}}" -gt 0 ]; then
                fastp \\
                    --in1 ${{tmp}}/{params.sid}.R1.fastq.gz \\
                    --in2 ${{tmp}}/{params.sid}.R2.fastq.gz \\
                    {params.qc_only} \\
                    --json {output.json} \\
                    --html {output.html} \\
                    --thread {params.fastp_threads}
            else
                fastp \\
                    --in1 ${{tmp}}/{params.sid}.unpaired.fastq.gz \\
                    {params.qc_only} \\
                    --json {output.json} \\
                    --html {output.html} \\
                    --thread {params.fastp_threads}
            fi
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


rule make_intervals:
    """
    Tile the reference genome into fixed-size windows, as the interval BED the
    interval-based analyses summarize over.

    These windows used to be built by hand and their path stored in
    config/genome.json, which meant the size was fixed per genome build and
    changing it was an out-of-band edit to a shared reference. They are now
    derived from the build's reference FastA at the width the run asked for
    (--interval), so the window size is a property of the run and is recorded in
    both config.json and the output filename.

    Gathers over nothing and scatters over nothing: one BED per run, shared by
    every sample and every interval-based analysis downstream. The FastA is
    streamed a line at a time to sum contig lengths, so the memory cost is flat
    regardless of genome size, but a ~3 GB assembly still takes a few minutes to
    read through.

    @Input:
        Reference FastA for the genome build (project-level).
    @Output:
        BED of contiguous, non-overlapping windows tiling every contig in the
        FastA, named for the width they were built at.
    """
    output:
        bed                     = intervals,
    container:
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "make_intervals", cluster),
        mem       = allocated("mem",  "make_intervals", cluster),
        time      = allocated("time", "make_intervals", cluster),
        gres      = allocated("gres", "make_intervals", cluster),
    threads:
        int(allocated("threads", "make_intervals", cluster))
    params:
        rname                   = "make_intervals",
        python_script           = join(bin_dir, 'fragment_genome.py'),
        # Declared as a param rather than an input, matching every other
        # reference artifact in this workflow. The FastA lives on the shared
        # reference filesystem and is not produced by any rule, so listing it as
        # an input only makes DAG construction fail wherever that filesystem is
        # absent (e.g. CI dry-runs).
        reference_fa            = reference_fa,
        interval_width          = interval_width,
    shell:
        dedent("""
        python {params.python_script} \\
            {params.reference_fa} \\
            {params.interval_width} \\
            -o {output.bed}
        """)


rule frag_length_intervals:
    input:
        bam                     = join(bam_dir, "{sid}.sorted.bam"),
        # Generated by make_intervals rather than read from a bundled reference,
        # so unlike the other reference artifacts in this workflow it is a real
        # input: it is what puts make_intervals in the DAG ahead of this rule.
        intervals               = intervals,
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
    shell:
        dedent("""
        finaletoolkit \\
            frag-length-intervals {input.bam} {input.intervals} \\
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
        intervals               = intervals,
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
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        tmpdir                  = tmpdir
    shell:
        dedent("""
        if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
        tmp=$(mktemp -d -p \"{params.tmpdir}\")
        trap 'ls -al ${{tmp}}; rm -rf "${{tmp}}"' EXIT

        finaletoolkit \\
            interval-end-motifs {input.bam} {params.ref2bit} {input.intervals} \\
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
        intervals               = intervals,
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
        blacklist_cmd           = f" --blacklist {blacklist}" if blacklist else "",
        gap_cmd                 = f" -g {gap}" if gap else ""
    shell:
        dedent("""
        finaletoolkit delfi {input.bam} {params.chrom_sizes} {params.ref2bit} {input.intervals} \\
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


# Read QC reports, i.e. everything in the aggregate report that describes the
# reads rather than the alignments. Both input paths produce a FastQC report on
# the analysis BAM and a fastp report, but from different rules and with
# different meanings, so the list is assembled per path:
#   • illumina_fastq -> FastQC on the raw mates and on the BAM, the samtools
#     markdup duplicate report, and fastp's report on the trimming and
#     filtering it performed before alignment (align_fastq).
#   • bam            -> FastQC on the analysis BAM (fastqc_bam) and a
#     descriptive-only fastp report of the same BAM converted back to FastQ
#     (fastp_bam). There are no raw reads to report on, and duplicates are
#     taken from flagstat rather than marked here, so neither a pre-alignment
#     FastQC report nor a markdup report exists.
# Each is guaranteed to exist by the time multiqc runs, but they are requested
# explicitly so the aggregate report's dependency on them is visible in the
# DAG.
if input_type == "illumina_fastq":
    read_qc_reports = (
        expand(join(qc_dir, "{sample}.R1_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.R2_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.sorted_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.markdup.stats.txt"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.fastp.json"), sample=sample_stems)
    )
else:
    read_qc_reports = (
        expand(join(qc_dir, "{sample}.sorted_fastqc.zip"), sample=sample_stems)
        + expand(join(qc_dir, "{sample}.fastp.json"), sample=sample_stems)
    )


# finaletoolkit results to summarise in the aggregate report, gathered over all
# samples. The set is not fixed: most fragmentomics rules only exist when the
# selected genome build supplies the reference files they need, so the gating
# below mirrors the Snakefile's, target for target, and the summary covers
# whichever rules actually ran. Fragment length bins and coverage need nothing
# beyond the analysis BAM and are therefore always present.
#
# The keys double as command-line flags for scripts/finaletoolkit_multiqc.py
# (underscores become dashes), so a key rename has to be made in both places.
finaletoolkit_qc_inputs = {
    'frag_length_bins': expand(
        join(fragment_length_dir, "{sample}_frag_bin" + str(bin_size) + ".tsv"),
        sample=sample_stems
    ),
    'coverage': expand(join(coverage_dir, "{sample}_coverage.bed"), sample=sample_stems),
}

# 'reference_fa' stands in for the interval BED in the gating below: the BED is
# tiled from that FastA by make_intervals rather than shipped with the build, so
# its presence is what decides whether the interval-based analyses can run.
if split_interval and 'reference_fa' in genome_files:
    finaletoolkit_qc_inputs['frag_length_intervals'] = expand(
        join(fragment_length_int_dir, "{sample}_frag_interval.bed"), sample=sample_stems
    )

if 'ref2bit' in genome_files:
    finaletoolkit_qc_inputs['end_motifs'] = expand(
        join(end_motifs_dir, "{sample}_endmotif.tsv"), sample=sample_stems
    )
    finaletoolkit_qc_inputs['mds'] = expand(
        join(mds_dir, "{sample}_mds.tsv"), sample=sample_stems
    )
    if 'reference_fa' in genome_files:
        finaletoolkit_qc_inputs['interval_end_motifs'] = expand(
            join(interval_end_motifs_dir, "{sample}_endmotif_interval.tsv"),
            sample=sample_stems
        )
        if 'chrom_sizes' in genome_files:
            finaletoolkit_qc_inputs['delfi'] = expand(
                join(delfi_dir, "{sample}_delfi.bed"), sample=sample_stems
            )

if 'tss' in genome_files and 'tss_interval' in genome_files:
    finaletoolkit_qc_inputs['wps_aggr'] = expand(
        join(wps_dir, "{sample}_wps_out_tss_aggr.wig"), sample=sample_stems
    )
    if 'chrom_sizes' in genome_files:
        finaletoolkit_qc_inputs['adjusted_wps_aggr'] = expand(
            join(adjust_wps_dir, "{sample}_wps_out_tss_adj_aggr.wig"), sample=sample_stems
        )
        finaletoolkit_qc_inputs['cleavage_aggr'] = expand(
            join(cleavage_profile_dir, "{sample}_cleavage_profile_aggr.wig"), sample=sample_stems
        )

# The same dictionary as the argument list for the summary script.
finaletoolkit_qc_args = ' '.join(
    '--{0} {1}'.format(key.replace('_', '-'), ' '.join(paths))
    for key, paths in finaletoolkit_qc_inputs.items()
)


rule finaletoolkit_multiqc:
    """
    Summarise the finaletoolkit results as MultiQC custom content, so the
    fragmentomics measurements appear in the aggregate report alongside the
    read and alignment QC.

    MultiQC has no finaletoolkit module, and cannot be given one here, so
    nothing any fragmentomics rule produces would otherwise reach the report -
    its tables and wigs and bigwigs are simply files MultiQC does not
    recognise. This rule reads them and rewrites what is summarisable as
    custom content: the per-sample tables and profiles MultiQC renders
    natively, plus four headline numbers (median fragment length, short
    fragment fraction, MDS and mean coverage) that join the general statistics
    table next to the samtools and FastQC columns. See
    scripts/finaletoolkit_multiqc.py for the sections it writes and how each is
    derived.

    The output is a directory rather than a fixed file list because the set of
    sections varies with the genome build's reference files (see
    finaletoolkit_qc_inputs above); the script writes a section per input group
    it was given. MultiQC scans qc/ recursively, so the directory only has to
    land there to be picked up, and taking it as an input of multiqc keeps the
    ordering explicit.

    Only the files listed as inputs are read. Nothing here re-reads a BAM or
    recomputes a measurement, so the step stays cheap as the cohort grows; the
    per-interval end motif tables are the one large read, and the script takes
    those in chunks rather than loading a sample at a time.
    @Input:
        Per-sample finaletoolkit results, gathered over all samples: fragment
        length bins and coverage beds always, plus per-interval fragment
        lengths, end motifs (genome-wide and per interval), MDS, DELFI bins
        and the aggregate TSS profiles when the run produces them.
    @Output:
        Directory of MultiQC custom-content documents, one per section.
    """
    input:
        **finaletoolkit_qc_inputs
    output:
        directory(join(qc_dir, "finaletoolkit")),
    container:
        config['images']['finaletoolkit']
    resources:
        partition = allocated("partition", "finaletoolkit_multiqc", cluster),
        mem       = allocated("mem",  "finaletoolkit_multiqc", cluster),
        time      = allocated("time", "finaletoolkit_multiqc", cluster),
        gres      = allocated("gres", "finaletoolkit_multiqc", cluster),
    threads:
        int(allocated("threads", "finaletoolkit_multiqc", cluster))
    params:
        rname                   = "finaletoolkit_multiqc",
        python_script           = join(bin_dir, 'finaletoolkit_multiqc.py'),
        result_args             = finaletoolkit_qc_args,
    shell:
        dedent("""
        python {params.python_script} \\
            {params.result_args} \\
            --output {output}
        """)


rule multiqc:
    """
    Aggregate the per-sample QC reports into a single interactive HTML report.
    Only the qc/ directory is scanned, so MultiQC does not walk the large
    bigwig/bed outputs of the fragmentomics rules; what it shows of those comes
    from the custom-content sections finaletoolkit_multiqc writes into qc/. The
    coverage workbook is taken as an input so the report is generated once the
    project-level coverage gather has finished.
    @Input:
        Per-sample samtools stats/flagstat/idxstats reports and read QC reports
        (gather), the finaletoolkit custom-content directory, and the merged
        coverage workbook.
    @Output:
        MultiQC HTML report.
    """
    input:
        stats                   = expand(join(qc_dir, "{sample}.samtools.stats.txt"), sample=sample_stems),
        flagstat                = expand(join(qc_dir, "{sample}.flagstat.txt"), sample=sample_stems),
        idxstats                = expand(join(qc_dir, "{sample}.idxstats.txt"), sample=sample_stems),
        read_qc                 = read_qc_reports,
        finaletoolkit           = join(qc_dir, "finaletoolkit"),
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
