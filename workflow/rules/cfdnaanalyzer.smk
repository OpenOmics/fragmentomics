# =============================================================
# cfDNAanalyzer (CDA) feature extraction
# =============================================================
#
# cfDNAanalyzer (https://github.com/LiymLab/cfDNAanalyzer) is a second
# fragmentomics suite the pipeline can run over the same analysis BAMs the
# finaletoolkit rules in common.smk consume. It bundles eleven feature
# extractors, from copy number calling with ichorCNA to nucleosome profiling
# with Griffin, and the frontend's --cda-features selects which of them run.
#
# Only its feature *extraction* stage is used. cfDNAanalyzer also ships feature
# processing/selection and machine learning modules that fit classifiers over
# the extracted matrices, cross-validate them and emit per-sample class
# predictions and probabilities. None of that runs here. Passing --noDA skips
# it: that flag guards a block which closes on the final line of cfDNAanalyzer's
# driver, so every inference step it can perform is inside it. cda_extract
# additionally asserts afterwards that neither inference output directory was
# created, so the guarantee is checked per job rather than assumed. What this
# file produces are feature matrices and nothing predicted from them.
#
# The label column every matrix carries is a cfDNAanalyzer output format
# requirement rather than a real annotation, and nothing reads it: its driver
# refuses to assemble the matrices without a label file, and its matrix builder
# enumerates samples from that file, so a placeholder label of 0 is synthesized
# per sample purely to reach the extraction output. It is kept in the merged
# matrices so they stay drop-in comparable with a native cfDNAanalyzer run.
#
# Extraction is scattered one job per sample and gathered afterwards, which is
# the only way to get any parallelism out of cfDNAanalyzer: its driver walks the
# BAM list serially inside a single process. Every step that gathering has to
# undo the split for is per-sample-independent, so the merged matrices are the
# same ones a single all-samples invocation would have written:
#   • the matrix builder derives each matrix's columns from the sample's own
#     output files, and every feature's columns come from an input that is
#     shared across samples (the region BED, the bundled site lists, the
#     256 end motifs, the fixed genomic bins), so the columns agree
#   • PFE standardization z-scores each row (sample) against itself
#   • the NP site list tables are one row per sample, so they concatenate
#
# See docs/analyses/cfdnaanalyzer.md for the feature descriptions and the
# caveats worth knowing before selecting a feature.

# Python standard library
from os.path import join
from textwrap import dedent
# Local imports
from scripts.common import (
    allocated,
    cda_analyses,
    cda_genome_build,
    cda_matrices,
    cda_regions,
    provided,
    CDA_REGION_FEATURES
)


# configuration
# Features to extract, already filtered down to the ones this genome build has
# the reference files for. Everything in this file is gated on the list being
# non-empty, so a run that did not ask for cfDNAanalyzer defines none of it.
cda_features                     = cda_analyses(config, genome_files)
# The build name cfDNAanalyzer is handed, which for a build of the user's own is
# the assembly it declares itself coordinate compatible with rather than its own
# name. Never None where it is used: cda_analyses() returns nothing at all when
# there is no build cfDNAanalyzer can run against, so every rule below is gated
# off in exactly that case.
cda_genome                       = cda_genome_build(config, genome_files)
# Whether any of the requested features is measured over a region BED, which is
# what decides if the region BED has to be prepared at all.
cda_by_region                    = any(
                                     feature in CDA_REGION_FEATURES
                                     for feature in cda_features
                                   )
# Region BED the region-specific features are scored over: --cda-regions when
# given, otherwise the genome build's TSS intervals.
cda_region_source                = cda_regions(config, genome_files)
# Genome bin size for the CNA feature, in kilobases. Unrelated to --bin-size,
# which bins the fragment length histogram.
cda_cna_bin_size                 = int(config['options']['cda_cna_bin_size'])
# Flanks the TSSC feature averages coverage over, reusing the same TSS flanks
# the cleavage profile analysis is run with so both describe the same window.
cda_tss_left_flank               = left_flank
cda_tss_right_flank              = right_flank
# Every matrix the selected features write, which is what makes up this
# analysis's final targets.
cda_matrix_names                 = cda_matrices(cda_features)

# directories
cda_dir                          = join(output_dir, 'cfdnaanalyzer')
cda_sample_dir                   = join(cda_dir, 'samples')
cda_features_dir                 = join(cda_dir, 'features')
cda_merge_dir                    = join(cda_dir, '.merge')

# Regional EMR is qualitatively wider than every other matrix (15.6 million
# columns in the 72-sample production run).  It gets a streaming, split gather;
# the remaining matrices are safe to align one at a time in Python.
cda_regional_emr                 = 'EMR_region_motif_frequency'
cda_ordinary_matrix_names        = [
                                     matrix for matrix in cda_matrix_names
                                     if matrix != cda_regional_emr
                                   ]
cda_emr_merge_dir                = join(cda_merge_dir, cda_regional_emr)
cda_emr_schema                   = join(cda_emr_merge_dir, 'schema.txt')
cda_emr_header                   = join(cda_emr_merge_dir, 'header.csv')
cda_emr_rows_dir                 = join(cda_emr_merge_dir, 'rows')

# The BED3 the region-specific features are handed, normalized from
# cda_region_source by cda_regions_bed below.
cda_region_bed                   = join(cda_dir, 'regions.bed')

# Per-site-list nucleosome profile tables. The NP feature scores every
# transcription factor site list bundled with Griffin, and how many of those
# there are is a property of the image (or of a site list directory pointed at
# inside it), not something the workflow can enumerate up front, so the tables
# are tracked as a directory rather than as one target per site list.
cda_np_site_lists                = join(cda_features_dir, 'NP_site_list')


if cda_features and cda_by_region:
    rule cda_regions_bed:
        """
        Data-processing step to normalize the regions cfDNAanalyzer measures its
        region-specific features over into the sorted BED3 it expects. The
        source is --cda-regions when the user gave one and the genome build's
        TSS interval BED otherwise; either way only the first three columns are
        kept, so any BED flavor can be handed in.

        The BED is a run-level input rather than a bundled reference (it can be
        either a user file or a reference file), so it is read from params for
        the same reason every other reference path in this pipeline is: naming
        it as an input would make DAG construction fail wherever that
        filesystem is absent, i.e. in a CI dry-run.
        @Input:
            None (the source BED is a params path)
        @Output:
            Sorted BED3 of the regions to score (project-level)
        """
        output:
            bed                     = cda_region_bed,
        container:
            config['images']['cfdnaanalyzer']
        resources:
            partition = allocated("partition", "cda_regions_bed", cluster),
            mem       = allocated("mem",  "cda_regions_bed", cluster),
            time      = allocated("time", "cda_regions_bed", cluster),
            gres      = allocated("gres", "cda_regions_bed", cluster),
        threads:
            int(allocated("threads", "cda_regions_bed", cluster))
        params:
            rname                   = "cda_regions_bed",
            source                  = cda_region_source,
        shell:
            dedent("""
            # Track lines and comments are dropped rather than passed through:
            # every consumer of this file reads it as plain three-column
            # intervals, and a header line would be scored as a region.
            grep -v -e '^#' -e '^track' -e '^browser' {params.source} \\
                | cut -f1,2,3 \\
                | sort -k1,1V -k2,2n \\
                > {output.bed}

            echo "Regions to score: $(wc -l < {output.bed})"
            """)


if cda_features:
    rule cda_extract:
        """
        Extract the selected cfDNAanalyzer features for one sample. This is the
        scatter half of the analysis: cfDNAanalyzer is invoked with a one-line
        BAM list so each sample is a separate job, and its per-sample feature
        matrices are collected by the split gather rules afterwards.

        Three things about how cfDNAanalyzer is packaged shape the shell below,
        and none of them are optional:

          1. Sample naming. cfDNAanalyzer names a sample after its BAM with the
             .bam suffix stripped, so the analysis BAM is linked in as
             <sample>.bam. Handing it bams/<sample>.sorted.bam directly would
             label every result "<sample>.sorted".

          2. The installation tree is written to, so CDA is not invoked in
             place. The image puts a CDA symlink on PATH, and the driver it
             points at derives everything it uses -- bundled references,
             helpers, binaries -- from its own location, resolved with
             `readlink -f "$0"`. Because that call follows the symlink all the
             way back to the real installation, invoking CDA directly would
             point the driver at a read-only tree, and the NP feature writes
             generated Griffin snakemake configs and workflows *inside* it.
             The rule therefore resolves CDA to find the installation, then
             builds a shadow copy of it in scratch space: the driver is copied
             (linking it would resolve back out again), everything beside it is
             symlinked, and only Griffin's snakemakes directory is copied, since
             that is the part written to.

          3. The NP feature runs snakemake itself. A nested snakemake keeps its
             state in ./.snakemake and is invoked with --unlock, so it has to
             run somewhere other than the pipeline's working directory; leaving
             the working directory as-is would point that --unlock at the outer
             workflow's own lock. The whole invocation is therefore run from
             scratch space.

        Note that cfDNAanalyzer re-filters the BAM it is given (MAPQ 30, and
        unmapped/secondary/QC-fail/duplicate reads dropped) before extracting
        anything. That is hard-coded in its driver and cannot be configured, so
        its features are measured on a slightly stricter subset of the analysis
        BAM than the finaletoolkit analyses are, and --mapscore does not apply
        to them. Its own filter is a floor, not a ceiling: a --mapscore above 30
        still holds, because the analysis BAM was already filtered to it.
        @Input:
            Analysis BAM and its index (scatter-per-sample)
        @Output:
            One CSV per feature matrix for this sample, plus the per-site-list
            nucleosome profile tables when NP was selected
        """
        input:
            bam                     = join(bam_dir, "{sid}.sorted.bam"),
            bai                     = join(bam_dir, "{sid}.sorted.bam.bai"),
            flagstat                = join(qc_dir, "{sid}.flagstat.txt"),
            regions                 = provided([cda_region_bed], cda_by_region),
        output:
            csvs                    = expand(
                                        join(cda_sample_dir, "{{sid}}", "{matrix}.csv"),
                                        matrix=cda_matrix_names
                                      ),
            sites                   = provided(
                                        [directory(join(cda_sample_dir, "{sid}", "NP_site_list"))],
                                        'NP' in cda_features
                                      ),
        container:
            config['images']['cfdnaanalyzer']
        resources:
            partition = allocated("partition", "cda_extract", cluster),
            mem       = allocated("mem",  "cda_extract", cluster),
            time      = allocated("time", "cda_extract", cluster),
            gres      = allocated("gres", "cda_extract", cluster),
        threads:
            int(allocated("threads", "cda_extract", cluster))
        params:
            rname                   = "cda_extract",
            tmpdir                  = tmpdir,
            sample_dir              = join(cda_sample_dir, "{sid}"),
            features                = ','.join(cda_features),
            matrices                = ' '.join(cda_matrix_names),
            genome                  = cda_genome,
            cna_bin_size            = cda_cna_bin_size,
            left_flank              = cda_tss_left_flank,
            right_flank             = cda_tss_right_flank,
            fasta_option            = "-f {0}".format(reference_fa) if reference_fa else "",
            regions_option          = "-b {0}".format(cda_region_bed) if cda_by_region else "",
        shell:
            dedent("""
            if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
            tmp=$(mktemp -d -p \"{params.tmpdir}\")
            trap 'rm -rf \"${{tmp}}\"' EXIT

            # cfDNAanalyzer takes a file of BAM paths, one per line, and names
            # each sample after the basename of its BAM without the .bam
            # suffix. Link the analysis BAM in under the sample name so that
            # is what its output is labelled with.
            mkdir -p \"${{tmp}}/bams\"
            ln -s \"$(readlink -f {input.bam})\" \"${{tmp}}/bams/{wildcards.sid}.bam\"
            ln -s \"$(readlink -f {input.bai})\" \"${{tmp}}/bams/{wildcards.sid}.bam.bai\"
            echo \"${{tmp}}/bams/{wildcards.sid}.bam\" > \"${{tmp}}/bams.txt\"

            # Placeholder labels. --noDA skips every analysis that would read
            # them, but the driver still exits without a label file and its
            # matrix builder enumerates samples from it, so one is synthesized
            # to reach the feature matrices.
            printf 'sample,label\\n{wildcards.sid},0\\n' > \"${{tmp}}/labels.csv\"

            # Read layout, taken from the report bam_stats already wrote for
            # this BAM, the same way fastp_bam does. This is not cosmetic: -s
            # decides whether cfDNAanalyzer treats the fragments as pairs, and
            # the TSSC feature changes how it extends reads based on it. Seven
            # of the eleven features (EM, FP, NP, OCF, EMR, FPR, PFE) are
            # paired-end only and cfDNAanalyzer refuses to run them under
            # 'single', so a single-end run can only select CNA, NOF, WPS and
            # TSSC. Empty (no such line in the report) is treated as single-end.
            paired=$(awk '/paired in sequencing/ {{print $1; exit}}' {input.flagstat})
            if [ \"${{paired:-0}}\" -gt 0 ]; then
                sequencing=\"pair\"
            else
                sequencing=\"single\"
            fi
            echo \"Read layout: ${{sequencing}}\"

            # Locate the installation through the CDA entrypoint the image puts
            # on PATH, rather than hardcoding a path here. CDA is a symlink to
            # the driver, and the driver derives every bundled reference, helper
            # and binary it uses from its own resolved location, so resolving
            # the link gives both the command to run and the tree it reads.
            cda_bin=$(readlink -f \"$(command -v CDA)\")
            cda_root=$(dirname \"${{cda_bin}}\")
            echo \"cfDNAanalyzer installation: ${{cda_root}}\"

            # Shadow installation, so the NP feature has somewhere writable to
            # generate its Griffin workflows.
            #
            # CDA cannot simply be invoked in place. The driver locates itself
            # with `readlink -f \"$0\"`, which resolves the symlink all the way
            # back to the real installation, so calling CDA directly would set
            # its script_dir to that read-only tree and NP would fail writing
            # its generated configs and snakefiles into Griffin/snakemakes.
            # The driver is therefore *copied* into scratch (a symlink would
            # resolve straight back out again), keeping the CDA name so the
            # logs still show what was run; everything beside it is symlinked,
            # apart from Griffin's snakemakes directory, which is the part that
            # gets written to and so has to be a real copy.
            mkdir -p \"${{tmp}}/cda/Griffin\"
            cp \"${{cda_bin}}\" \"${{tmp}}/cda/CDA\"
            driver=$(basename \"${{cda_bin}}\")
            for entry in \"${{cda_root}}\"/*; do
                base=$(basename \"${{entry}}\")
                if [ \"${{base}}\" = \"${{driver}}\" ] || [ \"${{base}}\" = \"Griffin\" ]; then
                    continue
                fi
                ln -s \"${{entry}}\" \"${{tmp}}/cda/${{base}}\"
            done
            for entry in \"${{cda_root}}\"/Griffin/*; do
                base=$(basename \"${{entry}}\")
                if [ \"${{base}}\" = \"snakemakes\" ]; then
                    continue
                fi
                ln -s \"${{entry}}\" \"${{tmp}}/cda/Griffin/${{base}}\"
            done
            cp -r \"${{cda_root}}/Griffin/snakemakes\" \"${{tmp}}/cda/Griffin/snakemakes\"
            # cp preserves the source's mode, and the whole point of this copy
            # is that NP writes into it, so make sure it is writable regardless
            # of how the installation in the image was permissioned.
            chmod -R u+w \"${{tmp}}/cda/Griffin/snakemakes\"

            # Make libR.so findable before anything runs. The NOF feature drives
            # DANPOS3, which reaches R through rpy2, and rpy2 embeds R in the
            # Python process instead of going through R's launcher script. That
            # launcher is what normally puts $R_HOME/lib on the loader's search
            # path -- it sources $R_HOME/etc/ldpaths, which prepends it to
            # LD_LIBRARY_PATH -- and the directory is on no default path
            # otherwise: no /etc/ld.so.conf.d entry, so libR.so is not in the
            # ldconfig cache, and the R binary carries no RUNPATH. Embedded, R
            # therefore starts and then fails to dyn.load its own base packages,
            # every one of which links libR.so, and NOF dies before writing a
            # single wig while the ten other features carry on.
            #
            # The image sets this too, but exporting it here as well means NOF
            # works against a published tag that predates that, and costs
            # nothing once it does not. $R_HOME/lib holds only libR.so, so
            # nothing else is shadowed by this.
            export LD_LIBRARY_PATH=\"$(R RHOME)/lib${{LD_LIBRARY_PATH:+:${{LD_LIBRARY_PATH}}}}\"

            # Run from scratch space: the NP feature's nested snakemake writes
            # its state to ./.snakemake and calls --unlock on it, which must
            # not be the outer workflow's working directory.
            cd \"${{tmp}}\"
            bash \"${{tmp}}/cda/CDA\" \\
                -I \"${{tmp}}/bams.txt\" \\
                -o \"${{tmp}}/cda/out\" \\
                -F {params.features} \\
                -g {params.genome} \\
                {params.fasta_option} \\
                -s \"${{sequencing}}\" \\
                -t {threads} \\
                -B {params.cna_bin_size} \\
                -u {params.left_flank} \\
                -d {params.right_flank} \\
                {params.regions_option} \\
                --noDA \\
                --labelFile \"${{tmp}}/labels.csv\"

            # Extraction only, checked rather than assumed. --noDA is the
            # boundary in cfDNAanalyzer's driver: the guard it opens closes on
            # the driver's last line, so every classifier fit, cross-validation
            # and prediction it can perform sits inside it and writes to one of
            # the two directories below. Their absence is what makes this rule
            # feature extraction and nothing else, so it is asserted here: if a
            # future image ever ships a driver that ignores the flag, this fails
            # the job rather than quietly publishing inferred output.
            for inferred in Feature_Processing_and_Selection Machine_Learning; do
                if [ -e \"${{tmp}}/cda/out/${{inferred}}\" ]; then
                    echo \"Error: cfDNAanalyzer produced ${{inferred}} despite\" \\
                         \"--noDA. This pipeline extracts features only and\" \\
                         \"will not publish inferred output.\" >&2
                    exit 1
                fi
            done

            # Publish only the matrices this run asked for. Each is normally a
            # two-line CSV, a header and this sample's row.
            #
            # A matrix can be absent entirely, though, and absorbing that is
            # what keeps a sample cfDNAanalyzer could not score for one feature
            # from failing the whole job.
            #
            # PFE is the case this was written for, and it is worth being
            # precise about how its CSV disappears, because it is not the
            # all-NA clean-up: that leaves the header behind, exactly as one
            # would hope. It is upstream's PFE standardization, which is three
            # unguarded lines (cfDNAanalyzer.sh, in the `PFE` branch after
            # Data_transformation.py):
            #
            #     Rscript .../standard_PFE.R Features/PFE.csv PFE_standard.csv
            #     rm Features/PFE.csv
            #     mv PFE_standard.csv Features/PFE.csv
            #
            # standard_PFE.R row-scales the matrix, so on the header-only CSV a
            # dropped sample leaves it halts in `colnames<-` ('names' attribute
            # [1] must be the same length as the vector [0]) and writes no
            # PFE_standard.csv. The `rm` runs regardless, the `mv` then has
            # nothing to move, and since the driver does not set -e neither
            # failure stops it. The header-only file is deleted and nothing
            # replaces it.
            #
            # The consequence to keep in mind: that `rm` is unconditional, so
            # *any* failure of standard_PFE.R destroys PFE for the sample, not
            # just this one. The warning below therefore says the matrix is
            # missing rather than asserting why -- quality control is the
            # likely reason, not a proven one.
            #
            # An absent matrix is published as an empty file: the output
            # Snakemake was promised exists, and the gather reads it
            # as a sample with no rows, which is what it is. The sample is then
            # simply missing from the merged matrix, as it would have been had
            # cfDNAanalyzer scored every sample at once.
            #
            # Features/ itself missing is a different thing -- it means the
            # driver never got as far as assembling any matrix -- so that is
            # asserted rather than absorbed, and no amount of quality control
            # can explain it.
            if [ ! -d \"${{tmp}}/cda/out/Features\" ]; then
                echo \"Error: cfDNAanalyzer wrote no Features directory, so it\" \\
                     \"assembled no feature matrices at all. This is a failed\" \\
                     \"run rather than a sample dropped by quality control.\" >&2
                exit 1
            fi
            mkdir -p \"{params.sample_dir}\"
            for matrix in {params.matrices}; do
                if [ -f \"${{tmp}}/cda/out/Features/${{matrix}}.csv\" ]; then
                    cp \"${{tmp}}/cda/out/Features/${{matrix}}.csv\" \\
                        \"{params.sample_dir}/${{matrix}}.csv\"
                else
                    echo \"Warning: cfDNAanalyzer left no ${{matrix}}.csv, so it\" \\
                         \"measured nothing for this sample -- most likely it\" \\
                         \"failed that feature's quality control, though for\" \\
                         \"PFE a standard_PFE.R failure looks the same.\" \\
                         \"Publishing it empty; the sample will be absent from\" \\
                         \"the merged matrix.\" >&2
                    : > \"{params.sample_dir}/${{matrix}}.csv\"
                fi
            done

            # The per-site-list nucleosome profile tables, which only exist
            # when the NP feature was among the ones selected.
            if [ -d \"${{tmp}}/cda/out/Features/NP_site_list\" ]; then
                cp -r \"${{tmp}}/cda/out/Features/NP_site_list\" \\
                    \"{params.sample_dir}/NP_site_list\"
            fi
            """)


    if cda_ordinary_matrix_names:
        rule cda_merge_matrix:
            """
            Merge one ordinary cfDNAanalyzer matrix across all samples. Equal
            schemas are byte-streamed; unequal schemas are aligned by name one
            sample at a time. One matrix per job bounds memory and lets the
            independent gathers run concurrently.
            @Input:
                One per-sample CSV for this matrix
            @Output:
                One project-level feature matrix CSV
            """
            input:
                csvs                    = lambda wildcards: [
                                            join(
                                              cda_sample_dir,
                                              sample,
                                              wildcards.matrix + '.csv'
                                            )
                                            for sample in sample_stems
                                          ],
            output:
                csv                     = join(cda_features_dir, "{matrix}.csv"),
            wildcard_constraints:
                matrix                  = '|'.join(cda_ordinary_matrix_names),
            container:
                config['images']['cfdnaanalyzer']
            resources:
                partition = allocated("partition", "cda_merge_matrix", cluster),
                mem       = allocated("mem",  "cda_merge_matrix", cluster),
                time      = allocated("time", "cda_merge_matrix", cluster),
                gres      = allocated("gres", "cda_merge_matrix", cluster),
            threads:
                int(allocated("threads", "cda_merge_matrix", cluster))
            params:
                rname                   = "cda_merge_matrix",
                tmpdir                  = tmpdir,
                python_script           = join(bin_dir, 'merge_cda_features.py'),
                sources                 = lambda wildcards: ' '.join(
                                            join(
                                              cda_sample_dir,
                                              sample,
                                              wildcards.matrix + '.csv'
                                            )
                                            for sample in sample_stems
                                          ),
                samples                 = ' '.join(sample_stems),
            shell:
                dedent("""
                if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
                tmp=$(mktemp -d -p \"{params.tmpdir}\")
                trap 'rm -rf \"${{tmp}}\"' EXIT

                python {params.python_script} matrix \\
                    --name {wildcards.matrix} \\
                    --sources {params.sources} \\
                    --samples {params.samples} \\
                    --scratch \"${{tmp}}\" \\
                    --output {output.csv}
                """)


    if 'NP' in cda_features:
        rule cda_merge_np_site_lists:
            """Merge the per-site-list nucleosome profile tables."""
            input:
                sites                   = expand(
                                            join(cda_sample_dir, "{sample}", "NP_site_list"),
                                            sample=sample_stems
                                          ),
            output:
                sites                   = directory(cda_np_site_lists),
            container:
                config['images']['cfdnaanalyzer']
            resources:
                partition = allocated("partition", "cda_merge_np_site_lists", cluster),
                mem       = allocated("mem",  "cda_merge_np_site_lists", cluster),
                time      = allocated("time", "cda_merge_np_site_lists", cluster),
                gres      = allocated("gres", "cda_merge_np_site_lists", cluster),
            threads:
                int(allocated("threads", "cda_merge_np_site_lists", cluster))
            params:
                rname                   = "cda_merge_np_site_lists",
                tmpdir                  = tmpdir,
                python_script           = join(bin_dir, 'merge_cda_features.py'),
                sample_dirs             = ' '.join(
                                            join(cda_sample_dir, sample)
                                            for sample in sample_stems
                                          ),
                samples                 = ' '.join(sample_stems),
                output_dir              = cda_features_dir,
            shell:
                dedent("""
                if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
                tmp=$(mktemp -d -p \"{params.tmpdir}\")
                trap 'rm -rf \"${{tmp}}\"' EXIT

                python {params.python_script} site-lists \\
                    --sample-dirs {params.sample_dirs} \\
                    --samples {params.samples} \\
                    --scratch \"${{tmp}}\" \\
                    --output {params.output_dir}
                """)


    if cda_regional_emr in cda_matrix_names:
        rule cda_emr_schema:
            """
            Stream the sorted regional-EMR headers into their union schema.
            Only one feature name per sample is retained in memory.

            cda_dense_merge is the streaming helper, compiled into the
            cfDNAanalyzer image and resolved on PATH; see
            docker/cfDNAanalyzer/cda_dense_merge.cpp for its two modes.
            """
            input:
                csvs                    = expand(
                                            join(
                                              cda_sample_dir,
                                              "{sample}",
                                              cda_regional_emr + '.csv'
                                            ),
                                            sample=sample_stems
                                          ),
            output:
                schema                  = temp(cda_emr_schema),
                header                  = temp(cda_emr_header),
            container:
                config['images']['cfdnaanalyzer']
            resources:
                partition = allocated("partition", "cda_emr_schema", cluster),
                mem       = allocated("mem",  "cda_emr_schema", cluster),
                time      = allocated("time", "cda_emr_schema", cluster),
                gres      = allocated("gres", "cda_emr_schema", cluster),
            threads:
                int(allocated("threads", "cda_emr_schema", cluster))
            params:
                rname                   = "cda_emr_schema",
                tmpdir                  = tmpdir,
                python_script           = join(bin_dir, 'merge_cda_features.py'),
                sources                 = ' '.join(
                                            join(
                                              cda_sample_dir,
                                              sample,
                                              cda_regional_emr + '.csv'
                                            )
                                            for sample in sample_stems
                                          ),
                samples                 = ' '.join(sample_stems),
            shell:
                dedent("""
                if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
                tmp=$(mktemp -d -p \"{params.tmpdir}\")
                trap 'rm -rf \"${{tmp}}\"' EXIT

                python {params.python_script} source-list \\
                    --sources {params.sources} \\
                    --samples {params.samples} \\
                    --output \"${{tmp}}/sources.tsv\"

                if [ -s \"${{tmp}}/sources.tsv\" ]; then
                    cda_dense_merge schema \\
                        --list \"${{tmp}}/sources.tsv\" \\
                        --schema \"${{tmp}}/schema.txt\" \\
                        --header \"${{tmp}}/header.csv\"
                else
                    : > \"${{tmp}}/schema.txt\"
                    printf 'sample,label\\n' > \"${{tmp}}/header.csv\"
                fi

                mkdir -p \"$(dirname {output.schema})\"
                mv \"${{tmp}}/schema.txt\" {output.schema}
                mv \"${{tmp}}/header.csv\" {output.header}
                """)


        rule cda_emr_row:
            """
            Align one regional-EMR sample row to the union schema. Splitting by
            sample makes the 15-million-column transform parallel and bounded.

            The row is written in place rather than staged through scratch: it
            is the only thing this job produces, it is a temp() output only the
            assemble rule consumes, and cda_dense_merge either writes it whole
            or exits non-zero, leaving Snakemake to remove the partial file.
            """
            input:
                csv                     = join(
                                            cda_sample_dir,
                                            "{sample}",
                                            cda_regional_emr + '.csv'
                                          ),
                schema                  = cda_emr_schema,
            output:
                row                     = temp(join(cda_emr_rows_dir, "{sample}.csv")),
            container:
                config['images']['cfdnaanalyzer']
            resources:
                partition = allocated("partition", "cda_emr_row", cluster),
                mem       = allocated("mem",  "cda_emr_row", cluster),
                time      = allocated("time", "cda_emr_row", cluster),
                gres      = allocated("gres", "cda_emr_row", cluster),
            threads:
                int(allocated("threads", "cda_emr_row", cluster))
            params:
                rname                   = "cda_emr_row",
            shell:
                dedent("""
                mkdir -p \"$(dirname {output.row})\"

                if [ ! -s {input.csv} ]; then
                    : > {output.row}
                else
                    cda_dense_merge row \\
                        --source {input.csv} \\
                        --sample {wildcards.sample} \\
                        --schema {input.schema} \\
                        --output {output.row}
                fi
                """)


        rule cda_emr_assemble:
            """Assemble the regional-EMR header and aligned sample rows."""
            input:
                header                  = cda_emr_header,
                rows                    = expand(
                                            join(cda_emr_rows_dir, "{sample}.csv"),
                                            sample=sample_stems
                                          ),
            output:
                csv                     = join(
                                            cda_features_dir,
                                            cda_regional_emr + '.csv'
                                          ),
            container:
                config['images']['cfdnaanalyzer']
            resources:
                partition = allocated("partition", "cda_emr_assemble", cluster),
                mem       = allocated("mem",  "cda_emr_assemble", cluster),
                time      = allocated("time", "cda_emr_assemble", cluster),
                gres      = allocated("gres", "cda_emr_assemble", cluster),
            threads:
                int(allocated("threads", "cda_emr_assemble", cluster))
            params:
                rname                   = "cda_emr_assemble",
                tmpdir                  = tmpdir,
                python_script           = join(bin_dir, 'merge_cda_features.py'),
                rows                    = ' '.join(
                                            join(cda_emr_rows_dir, sample + '.csv')
                                            for sample in sample_stems
                                          ),
            shell:
                dedent("""
                if [ ! -d \"{params.tmpdir}\" ]; then mkdir -p \"{params.tmpdir}\"; fi
                tmp=$(mktemp -d -p \"{params.tmpdir}\")
                trap 'rm -rf \"${{tmp}}\"' EXIT

                python {params.python_script} assemble \\
                    --name {cda_regional_emr} \\
                    --header {input.header} \\
                    --rows {params.rows} \\
                    --scratch \"${{tmp}}\" \\
                    --output {output.csv}
                """)


    # The gather is intentionally split into the independent rules above.
    # The previous single-job implementation loaded all matrices into pandas at
    # once and could not complete the production regional-EMR join.
