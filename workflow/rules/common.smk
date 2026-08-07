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
#   • illumina_fastq -> align_fastq: bwa-mem2 alignment + dedup + proper-pair
#     filtering (full pipeline work).
#   • bam            -> stage_bams:  sort/index ready-made alignments (full
#     pipeline work minus alignment).
if input_type == "illumina_fastq":

    rule align_fastq:
        """
        Align paired-end Illumina FastQ reads to the reference genome with
        bwa-mem2, then produce an analysis-ready coordinate-sorted BAM:
        fixmate -> sort -> markdup (mark duplicates) -> keep only properly
        paired, primary, mapped, non-duplicate reads. This yields the same
        bams/{sid}.sorted.bam that the BAM input path stages directly, so all
        downstream fragmentomics rules are input-type agnostic.
        @Input:
            Paired-end FastQ mates (scatter-per-sample), staged by the
            frontend into the output directory's inputs/ folder.
        @Output:
            Coordinate-sorted, indexed analysis BAM.
        """
        input:
            r1                  = join(inputs_dir, "{sid}.R1.fastq.gz"),
            r2                  = join(inputs_dir, "{sid}.R2.fastq.gz"),
        output:
            bam                 = join(bam_dir, "{sid}.sorted.bam"),
            bai                 = join(bam_dir, "{sid}.sorted.bam.bai"),
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
            # 0x2 (proper pair) kept; 3852 excludes unmapped, mate-unmapped,
            # secondary, qcfail, duplicate, and supplementary reads.
            keep_flag           = "2",
            drop_flag           = "3852",
            tmpdir              = tmpdir,
        shell:
            dedent("""
            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf "${{tmp}}"' EXIT

            bwa-mem2 mem \\
                -t {threads} \\
                -R "@RG\\tID:{params.sid}\\tSM:{params.sid}\\tPL:ILLUMINA\\tLB:{params.sid}" \\
                {params.index} {input.r1} {input.r2} \\
            | samtools sort -n -@ {threads} -T ${{tmp}}/nsort -O bam - \\
            | samtools fixmate -m -@ {threads} - - \\
            | samtools sort -@ {threads} -T ${{tmp}}/psort -O bam - \\
            | samtools markdup -@ {threads} -T ${{tmp}}/mkdup - - \\
            | samtools view -b -f {params.keep_flag} -F {params.drop_flag} \\
                -@ {threads} -o {output.bam} -

            samtools index -@ {threads} {output.bam}
            """)

else:

    # Where stage_bams writes. When the selected genome ships a sequence
    # dictionary, staging lands in staged_bam_dir and filter_reference_contigs
    # produces the canonical bams/{sid}.sorted.bam. Without a dictionary there
    # is nothing to validate against, so staging writes bam_dir directly and
    # the filter rule is not defined at all.
    stage_dir = staged_bam_dir if sequence_dict else bam_dir

    rule stage_bams:
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
            memory                   = str(allocated("mem",  "stage_bams", cluster)).replace('G' ,'')
        shell:
            dedent("""
            python {params.python_script} \\
                --files {input} \\
                --output {params.bam_dir} \\
                --threads {threads} \\
                --memory {params.memory}
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
            seq_dict            = sequence_dict,
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
        shell:
            dedent("""
            python {params.python_script} \\
                --bam {input.bam} \\
                --dict {input.seq_dict} \\
                --output {output.bam} \\
                --report {output.report} \\
                --threads {threads}
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
        intervals               = tss_interval,
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len
    shell:
        dedent("""
        finaletoolkit \\
            coverage {input.bam} {params.intervals} \\
            -n \\
            --scale-factor 1000000 \\
            -q 30 \\
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
        bin_size                = str(bin_size),
        min_len                 = str(min_fragment_len),
        max_len                 = str(max_fragment_len),
    shell:
        dedent("""
        finaletoolkit \\
            frag-length-bins {input.bam} \\
            -q 30 \\
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
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        intervals               = intervals
    shell:
        dedent("""
        finaletoolkit \\
            frag-length-intervals {input.bam} {params.intervals} \\
            -q 30 \\
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
        min_len                 = min_fragment_len,
        max_len                 = max_fragment_len,
        ref2bit                 = ref2bit,
    shell:
        dedent("""
        finaletoolkit \\
            end-motifs {input.bam} {params.ref2bit} \\
            -q 30 \\
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
            -q 30 \\
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
        chrom_sizes             = chrom_sizes,
        ref2bit                 = ref2bit,
        intervals               = intervals,
        blacklist_cmd           = f" --blacklist {blacklist}" if blacklist else "",
        gap_cmd                 = f" -g {gap}" if gap else ""
    shell:
        dedent("""
        finaletoolkit delfi {input.bam} {params.chrom_sizes} {params.ref2bit} {params.intervals} \\
            -q 30{params.blacklist_cmd}{params.gap_cmd} \\
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
            -q 30 \\
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
            -q 30 \\
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


rule multiqc:
    """
    Aggregate the per-sample samtools QC reports into a single interactive
    HTML report. Only the qc/ directory is scanned, so MultiQC does not walk
    the large bigwig/bed outputs of the fragmentomics rules. The coverage
    workbook is taken as an input so the report is generated once the
    project-level coverage gather has finished.
    @Input:
        Per-sample samtools stats/flagstat/idxstats reports (gather) and the
        merged coverage workbook.
    @Output:
        MultiQC HTML report.
    """
    input:
        stats                   = expand(join(qc_dir, "{sample}.samtools.stats.txt"), sample=sample_stems),
        flagstat                = expand(join(qc_dir, "{sample}.flagstat.txt"), sample=sample_stems),
        idxstats                = expand(join(qc_dir, "{sample}.idxstats.txt"), sample=sample_stems),
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
