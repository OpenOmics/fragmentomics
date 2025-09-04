from scripts.common import abstract_location


rule coverage:
    input:
        
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        intervals               = GENOME_PARAMS["intervals"]
    output:
        bed                     = "coverage/{sample}_coverage.bed"
    threads: 
        config["resources"]["coverage"]["threads"]
    resources:
        rname                   = "preseq",
        mem                     = config["resources"]["coverage"]["mem"],
        time                    = config["resources"]["coverage"]["time"],
        partition               = config["resources"]["coverage"]["partition"]
    params:
        rname                   = "coverage",
        min_len                 = config["min"],
        max_len                 = config["max"]
    log:
        "logs/coverage/{sample}.log"
    shell:
        """
        mkdir -p coverage logs/coverage
        finaletoolkit coverage {input.bam} {input.intervals} -n --scale 1e6 -q 30 -min {params.min_len} -max {params.max_len} -p any -o {output.bed} -w {threads} -v > {log} 2>&1
        """


rule frag_length_bins:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam"
    output:
        tsv                     = "frag_length_bins/{sample}_frag_bin{bin_size}.tsv"
    threads: 
        config["resources"]["frag_length_bins"]["threads"]
    resources:
        rname                   = "frag_length_bins",
        mem                     = config["resources"]["frag_length_bins"]["mem"],
        time                    = config["resources"]["frag_length_bins"]["time"],
        partition               = config["resources"]["frag_length_bins"]["partition"]
    wildcard_constraints:
        bin_size                = "\d+"
    params:
        bin_size                = config["bin_size"],
        min_len                 = config["min"],
        max_len                 = config["max"],
        png                     = "frag_length_bins/{sample}_frag_bin{bin_size}.png"
    log:
        "logs/frag_length_bins/{sample}_bin{bin_size}.log"
    shell:
        """
        mkdir -p frag_length_bins logs/frag_length_bins
        finaletoolkit frag-length-bins {input.bam} -q 30 --bin-size {params.bin_size} -min {params.min_len} -max {params.max_len} -p midpoint -o {output.tsv} --histogram-path {params.png} -v > {log} 2>&1
        """


rule frag_length_intervals:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        intervals               = GENOME_PARAMS["intervals"]
    output:
        bed                     = "frag_length_intervals/{sample}_frag_interval.bed"
    threads: 
        config["resources"]["frag_length_intervals"]["threads"]
    resources:
        rname                   = "frag_length_intervals",
        mem                     = config["resources"]["frag_length_intervals"]["mem"],
        time                    = config["resources"]["frag_length_intervals"]["time"],
        partition               = config["resources"]["frag_length_intervals"]["partition"]
    params:
        min_len                 = config["min"],
        max_len                 = config["max"]
    log:
        "logs/frag_length_intervals/{sample}.log"
    shell:
        """
        mkdir -p frag_length_intervals logs/frag_length_intervals
        finaletoolkit frag-length-intervals {input.bam} {input.intervals} -q 30 -min {params.min_len} -max {params.max_len} -p any -o {output.bed} -w {threads} -v > {log} 2>&1
        """


rule end_motifs:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        ref2bit                 = GENOME_PARAMS["ref2bit"]
    output:
        tsv                     = "end_motifs/{sample}_endmotif.tsv"
    threads: 
        config["resources"]["end_motifs"]["threads"]
    resources:
        mem                     = config["resources"]["end_motifs"]["mem"],
        time                    = config["resources"]["end_motifs"]["time"],
        partition               = config["resources"]["end_motifs"]["partition"]
    params:
        min_len                 = config["min"],
        max_len                 = config["max"]
    log:
        "logs/end_motifs/{sample}.log"
    shell:
        """
        mkdir -p end_motifs logs/end_motifs
        finaletoolkit end-motifs {input.bam} {input.ref2bit} -q 30 -k 4 -min {params.min_len} -max {params.max_len} -o {output.tsv} -w {threads} -v > {log} 2>&1
        """


rule interval_end_motifs:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        ref2bit                 = GENOME_PARAMS["ref2bit"],
        intervals               = GENOME_PARAMS["intervals"]
    output:
        tsv                     = "interval_end_motifs/{sample}_endmotif_interval.tsv"
    threads: 
        config["resources"]["interval_end_motifs"]["threads"]
    resources:
        mem                     = config["resources"]["interval_end_motifs"]["mem"],
        time                    = config["resources"]["interval_end_motifs"]["time"],
        partition               = config["resources"]["interval_end_motifs"]["partition"]
    params:
        min_len                 = config["min"],
        max_len                 = config["max"]
    log:
        "logs/interval_end_motifs/{sample}.log"
    shell:
        """
        mkdir -p interval_end_motifs logs/interval_end_motifs
        finaletoolkit interval-end-motifs {input.bam} {input.ref2bit} {input.intervals} -q 30 -k 4 -min {params.min_len} -max {params.max_len} -o {output.tsv} -w {threads} -v > {log} 2>&1
        """


rule mds:
    input:
        endmotif                = "end_motifs/{sample}_endmotif.tsv"
    output:
        tsv                     = "mds/{sample}_mds.tsv"
    threads: 
        config["resources"]["mds"]["threads"]
    resources:
        mem                     = config["resources"]["mds"]["mem"],
        time                    = config["resources"]["mds"]["time"],
        partition               = config["resources"]["mds"]["partition"]
    params:
        sample                  = "{sample}"
    log:
        "logs/mds/{sample}.log"
    shell:
        """
        mkdir -p mds logs/mds
        mds_score=$(finaletoolkit mds {input.endmotif} 2> {log})
        echo "{params.sample}\t${{mds_score}}" > {output.tsv}
        """


rule delfi:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        chrom_sizes             = GENOME_PARAMS["chrom_sizes"],
        ref2bit                 = GENOME_PARAMS["ref2bit"],
        intervals               = GENOME_PARAMS["intervals"],
        blacklist               = GENOME_PARAMS["blacklist"] if GENOME_PARAMS.get("blacklist") else "",
        gap                     = GENOME_PARAMS["gap"] if GENOME_PARAMS.get("gap") else ""
    output:
        bed                     = "delfi/{sample}_delfi.bed"
    threads: 
        config["resources"]["delfi"]["threads"]
    resources:
        mem                     = config["resources"]["delfi"]["mem"],
        time                    = config["resources"]["delfi"]["time"],
        partition               = config["resources"]["delfi"]["partition"]
    params:
        blacklist_cmd           = lambda wildcards, input: f"--blacklist-file {input.blacklist}" if input.blacklist else "",
        gap_cmd                 = lambda wildcards, input: f"-g {input.gap}" if input.gap else ""
    log:
        "logs/delfi/{sample}.log"
    shell:
        """
        mkdir -p delfi logs/delfi
        finaletoolkit delfi {input.bam} {input.chrom_sizes} {input.ref2bit} {input.intervals} -q 30 -o {output.bed} -w {threads} -v --no-merge-bins {params.blacklist_cmd} {params.gap_cmd} > {log} 2>&1
        """


rule wps:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        tss                     = GENOME_PARAMS["tss"]
    output:
        bw                      = "wps/{sample}_wps_out_tss.bw"
    threads: 
        config["resources"]["wps"]["threads"]
    resources:
        mem                     = config["resources"]["wps"]["mem"],
        time                    = config["resources"]["wps"]["time"],
        partition               = config["resources"]["wps"]["partition"]
    params:
        i                       = config["i"]
    log:
        "logs/wps/{sample}.log"
    shell:
        """
        mkdir -p wps logs/wps
        finaletoolkit wps {input.bam} {input.tss} -i {params.i} -W 120 -min 120 -max 180 -q 30 -o {output.bw} -w {threads} -v > {log} 2>&1
        """


rule adjust_wps:
    input:
        wps_bw                  = "wps/{sample}_wps_out_tss.bw",
        tss_interval            = GENOME_PARAMS["tss_interval"],
        chrom_sizes             = GENOME_PARAMS["chrom_sizes"]
    output:
        bw                      = "adjust_wps/{sample}_wps_out_tss_adjusted.bw"
    threads: 
        config["resources"]["adjust_wps"]["threads"]
    resources:
        mem                     = config["resources"]["adjust_wps"]["mem"],
        time                    = config["resources"]["adjust_wps"]["time"],
        partition               = config["resources"]["adjust_wps"]["partition"]
    params:
        i                       = config["i"]
    log:
        "logs/adjust_wps/{sample}.log"
    shell:
        """
        mkdir -p adjust_wps logs/adjust_wps
        finaletoolkit adjust-wps {input.wps_bw} {input.tss_interval} {input.chrom_sizes} -o {output.bw} -i {params.i} -m 200 -S --subtract-edges -w {threads} -v > {log} 2>&1
        """


rule cleavage_profile:
    input:
        bam                     = lambda wildcards: f"{bam_dir}/{wildcards.sample}.bam",
        tss                     = GENOME_PARAMS["tss"],
        chrom_sizes             = GENOME_PARAMS["chrom_sizes"]
    output:
        bw                      = "cleavage_profile/{sample}_cleavage_profile_tss.bw"
    threads: 
        config["resources"]["cleavage_profile"]["threads"]
    resources:
        mem                     = config["resources"]["cleavage_profile"]["mem"],
        time                    = config["resources"]["cleavage_profile"]["time"],
        partition               = config["resources"]["cleavage_profile"]["partition"]
    params:
        l                       = config["l"],
        r                       = config["r"],
        min_len                 = config["min"],
        max_len                 = config["max"]
    log:
        "logs/cleavage_profile/{sample}.log"
    shell:
        """
        mkdir -p cleavage_profile logs/cleavage_profile
        finaletoolkit cleavage-profile {input.bam} {input.tss} -o {output.bw} -c {input.chrom_sizes} -l {params.l} -r {params.r} -q 30 -min {params.min_len} -max {params.max_len} -w {threads} -v > {log} 2>&1
        """

rule agg_wps:
    input:
        bw                      = "wps/{sample}_wps_out_tss.bw",
        tss_interval            = GENOME_PARAMS["tss_interval"]
    output:
        wig                     = "wps/{sample}_wps_out_tss_aggr.wig"
    threads: 
        config["resources"]["agg_wps"]["threads"]
    resources:
        mem                     = config["resources"]["agg_wps"]["mem"],
        time                    = config["resources"]["agg_wps"]["time"],
        partition               = config["resources"]["agg_wps"]["partition"]
    log:
        "logs/agg_wps/{sample}.log"
    shell:
        """
        mkdir -p logs/agg_wps
        finaletoolkit agg-bw {input.bw} {input.tss_interval} -o {output.wig} --mean -v > {log} 2>&1
        """

rule agg_adjust_wps:
    input:
        bw                      = "adjust_wps/{sample}_wps_out_tss_adjusted.bw",
        tss_interval            = GENOME_PARAMS["tss_interval"]
    output:
        wig                     = "adjust_wps/{sample}_wps_out_tss_adj_aggr.wig"
    threads: 
        config["resources"]["agg_adjust_wps"]["threads"]
    resources:
        mem                     = config["resources"]["agg_adjust_wps"]["mem"],
        time                    = config["resources"]["agg_adjust_wps"]["time"],
        partition               = config["resources"]["agg_adjust_wps"]["partition"]
    log:
        "logs/agg_adjust_wps/{sample}.log"
    shell:
        """
        mkdir -p logs/agg_adjust_wps
        finaletoolkit agg-bw {input.bw} {input.tss_interval} -o {output.wig} --mean -v > {log} 2>&1
        """

rule agg_cleavage_profile:
    input:
        bw                      = "cleavage_profile/{sample}_cleavage_profile_tss.bw",
        tss_interval            = GENOME_PARAMS["tss_interval"]
    output:
        wig                     = "cleavage_profile/{sample}_cleavage_profile_aggr.wig"
    threads: 
        config["resources"]["agg_cleavage_profile"]["threads"]
    resources:
        mem                     = config["resources"]["agg_cleavage_profile"]["mem"],
        time                    = config["resources"]["agg_cleavage_profile"]["time"],
        partition               = config["resources"]["agg_cleavage_profile"]["partition"]
    log:
        "logs/agg_cleavage_profile/{sample}.log"
    shell:
        """
        mkdir -p logs/agg_cleavage_profile
        finaletoolkit agg-bw {input.bw} {input.tss_interval} -o {output.wig} --mean -v > {log} 2>&1
        """