import os
from textwrap import dedent


# configuration
split_interval                  = int(config['options']['interval'])
max_fragment_len                = int(config['options']['max'])
min_fragment_len                = int(config['options']['min'])
right_flank                     = int(config['options']['right'])
left_flank                      = int(config['options']['left'])
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


# directories
data_dir                         = config["project"]["datapath"]
all_input_files                  = config['options']['input']
output_dir                       = config['options']['output']
bin_dir                          = config['binpath']
tmpdir                           = config['options']['tmp_dir']
bam_dir                          = os.path.join(output_dir, 'bams')
coverage_dir                     = os.path.join(output_dir, 'coverage')
fragment_length_dir              = os.path.join(output_dir, 'frag_length_bins')
fragment_length_int_dir          = os.path.join(output_dir, 'frag_length_intervals')
end_motifs_dir                   = os.path.join(output_dir, 'end_motifs')
interval_end_motifs_dir          = os.path.join(output_dir, 'interval_end_motifs')
mds_dir                          = os.path.join(output_dir, 'mds')
delfi_dir                        = os.path.join(output_dir, 'delfi')
wps_dir                          = os.path.join(output_dir, 'wps')
adjust_wps_dir                   = os.path.join(output_dir, 'adjust_wps')
cleavage_profile_dir             = os.path.join(output_dir, 'cleavage_profile')


# default resources
default_threads                  = config['cluster']['__default__']['threads']


rule stage_bams:
    input:
        all_input_files
    output:
        expand(os.path.join(bam_dir, "{sample}.sorted.bam"), sample=sample_stems)
    container: 
        config['images']['finaletoolkit']
    threads: 
        config["cluster"]["stage_bams"].get("threads", default_threads)
    params:
        rname                   = "stage_bams",
        bam_dir                  = bam_dir,
        python_script            = os.path.join(bin_dir, 'stage_input_files.py'),
        memory                   = str(config["cluster"]["stage_bams"].get("mem", default_threads)).replace('G' ,'')
    shell:
        dedent("""
        python {params.python_script} \\
            --files {input} \\
            --output {params.bam_dir} \\
            --threads {threads} \\
            --memory {params.memory}
        """)


rule coverage:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = os.path.join(coverage_dir, "{sid}_coverage.bed")
    threads: 
        config["cluster"]["coverage"].get("threads", default_threads)
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
            --scale 1e6 \\
            -q 30 \\
            -min {params.min_len} \\
            -max {params.max_len} \\
            -p any \\
            -o {output.bed} \\
            -w {threads} \\
            -v
        """)


rule frag_length_bins:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = os.path.join(fragment_length_dir, "{sid}_frag_bin" + str(bin_size) + ".tsv"),
        png                     = os.path.join(fragment_length_dir, "{sid}_frag_bin" + str(bin_size) + ".png")
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["frag_length_bins"].get("threads", default_threads)
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
            -min {params.min_len} \\
            -max {params.max_len} \\
            -p midpoint \\
            -o {output.tsv} \\
            --histogram-path {output.png} \\
            -v
        """)


rule frag_length_intervals:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = os.path.join(fragment_length_int_dir, "{sid}_frag_interval.bed")
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["frag_length_intervals"].get("threads", default_threads)
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
            -min {params.min_len} \\
            -max {params.max_len} \\
            -p any \\
            -o {output.bed} \\
            -w {threads} \\
            -v
        """)


rule end_motifs:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = os.path.join(end_motifs_dir, "{sid}_endmotif.tsv"),
    container: 
        config['images']['finaletoolkit']
    threads: 
        config["cluster"]["end_motifs"].get("threads", default_threads)
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
            -k 4 \\
            -min {params.min_len} \\
            -max {params.max_len} \\
            -o {output.tsv} \\
            -w {threads} \\
            -v
        """)


rule interval_end_motifs:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        tsv                     = os.path.join(interval_end_motifs_dir, "{sid}_endmotif_interval.tsv"),
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["interval_end_motifs"].get("threads", default_threads)
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
            -k 4 \\
            -min {params.min_len} \\
            -max {params.max_len} \\
            -o {output.tsv} \\
            -w {threads} \\
            -v
        """)


rule mds:
    input:
        endmotif                = os.path.join(end_motifs_dir, "{sid}_endmotif.tsv"),
    output:
        tsv                     = os.path.join(mds_dir, "{sid}_mds.tsv"),
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["mds"].get("threads", default_threads)
    params:
        rname                   = "mds",
        sid                     = "{sid}"
    shell:
        dedent("""
        mds_score=$(finaletoolkit mds {input.endmotif})
        echo "Sample\tMDS_score" > {output.tsv}
        echo "{params.sid}\t${{mds_score}}" > {output.tsv}
        """)


rule delfi:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        bed                     = os.path.join(delfi_dir, "{sid}_delfi.bed"),
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["delfi"].get("threads", default_threads)
    params:
        rname                   = "delfi",
        chrom_sizes             = chrom_sizes,
        ref2bit                 = ref2bit,
        intervals               = intervals,
        blacklist_cmd           = f"--blacklist-file {blacklist} " if blacklist else "",
        gap_cmd                 = f"-g {gap}" if gap else ""
    shell:
        """
        finaletoolkit delfi {input.bam} {params.chrom_sizes} {params.ref2bit} {params.intervals} \\
            -q 30 {params.blacklist_cmd}{params.gap_cmd} \\
            -o {output.bed} \\
            -w {threads} \\
            -v \\
            --no-merge-bins
        """


rule wps:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        bw                      = os.path.join(wps_dir, "{sid}_wps_out_tss.bw")
    container: 
        config['images']['finaletoolkit']
    threads: 
        config["cluster"]["wps"].get("threads", default_threads)
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
            -min 120 \\
            -max 180 \\
            -q 30 \\
            -o {output.bw} \\
            -w {threads} \\
            -v
        """


rule adjust_wps:
    input:
        wps_bw                  = os.path.join(wps_dir, "{sid}_wps_out_tss.bw")
    output:
        bw                      = os.path.join(adjust_wps_dir, "{sid}_wps_out_tss_adjusted.bw")
    container: 
        config['images']['finaletoolkit']
    threads: 
        config["cluster"]["adjust_wps"].get("threads", default_threads)
    params:
        rname                   = "adjust_wps",
        intervals               = split_interval,
        tss_interval            = tss_interval,
        chrom_sizes             = chrom_sizes
    shell:
        """
        finaletoolkit \\
            adjust-wps {input.wps_bw} {params.tss_interval} {params.chrom_sizes} \\
            -o {output.bw} \\
            -i {params.intervals} \\
            -m 200 \\
            -S \\
            --subtract-edges \\
            -v
        """


rule cleavage_profile:
    input:
        bam                     = os.path.join(bam_dir, "{sid}.sorted.bam"),
    output:
        bw                      = os.path.join(cleavage_profile_dir, "{sid}_cleavage_profile_tss.bw")
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["cleavage_profile"].get("threads", default_threads)
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
            -l {params.l} \\
            -r {params.r} \\
            -q 30 \\
            -min {params.min_len} \\
            -max {params.max_len} \\
            -w {threads} \\
            -v
        """)


rule agg_wps:
    input:
        bw                      = os.path.join(wps_dir, "{sid}_wps_out_tss.bw")
    output:
        wig                     = os.path.join(wps_dir, "{sid}_wps_out_tss_aggr.wig")
    container: 
        config['images']['finaletoolkit']
    threads: 
        config["cluster"]["agg_wps"].get("threads", default_threads)
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
        bw                      = os.path.join(adjust_wps_dir, "{sid}_wps_out_tss_adjusted.bw")
    output:
        wig                     = os.path.join(adjust_wps_dir, "{sid}_wps_out_tss_adj_aggr.wig")
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["agg_adjust_wps"].get("threads", default_threads)
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


rule agg_cleavage_profile:
    input:
        bw                      = os.path.join(cleavage_profile_dir, "{sid}_cleavage_profile_tss.bw"),
    output:
        wig                     = os.path.join(cleavage_profile_dir, "{sid}_cleavage_profile_aggr.wig"),
    container: 
        config['images']['finaletoolkit']
    threads:
        config["cluster"]["agg_cleavage_profile"].get("threads", default_threads)
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