# Common helper functions shared across the entire workflow
def provided(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    """
    if not condition:
        # If condition is False,
        # returns an empty list
        # to prevent rule from
        # running
        samplelist = []

    return samplelist


def ignore(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is met (i.e. True) then it will not try to run that rule. This function is the
    inverse to provided().
    """
    if condition:
        # If condition is True,
        # returns an empty list
        # to prevent rule from
        # running
        samplelist = []

    return samplelist


def references(config, reflist):
    """
    Checks if a set of required reference files were provided. Some rules depend
    on a set of required reference files that may only exist for specific reference
    genomes. An example of this would be blasklists arriba. The blacklist are manually
    curated and only exist for a few reference genomes (mm10, hg38, hg19).
    If one of the required reference files does not exist, then it will return
    an empty list.
    """
    _all = True
    for ref in reflist:
        try:
            tmp = config["references"][ref]
        # Check if ref exists in config
        except KeyError:
            _all = False
            break
        # Check if ref is empty key string
        if not tmp:
            _all = False

    return _all


def allocated(resource, rule, lookup, default="__default__"):
    """Pulls resource information for a given rule. If a rule does not have any information
    for a given resource type, then it will pull from the default. Information is pulled from
    definitions in the cluster.json (which is used a job submission). This ensures that any
    resources used at runtime mirror the resources that were allocated.
    :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
    :param rule <str>: rule to lookup its information
    :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
    :param default <str>: default information to use if rule information cannot be found
    :return allocation <str>:
        allocation information for a given resource type for a given rule
    """
    try:
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]

    return allocation


def str_bool(s):
    """Converts a string to boolean. It is dangerous to try to
    typecast a string into a boolean value using the built-in
    `bool()` function. This function avoids any issues that can
    arise when using `bool()`.
    Example:
      boolean('True') returns True
      boolean('False') returns False
      boolean('asdas') raises TypeError
    """
    val = s.lower()
    if val in ["true", "1", "y", "yes"]:
        return True
    elif val in ["false", "0", "n", "no", ""]:
        return False
    else:
        # Provided value could not be
        # type casted into a boolean
        raise TypeError(f"Fatal: cannot type cast {val} into a boolean")


def interval_label(bases):
    """Renders a window width in bases as the short SI-prefixed label used to
    name the generated interval BED. The frontend (--interval) validates and
    normalizes the user's base unit into a plain base count, so this is the one
    place a size becomes a filename and every equivalent spelling of a size
    ('2MB', '2m', '2000kb', '2000000') resolves to the same name.

    The largest prefix that divides the width evenly wins, so a width that is
    not a whole number of the bigger unit falls back to a smaller one rather
    than being rounded, i.e. 1500000 -> '1500kb', not '1mb' or '1.5mb'.
    Example:
      interval_label(1000000) returns '1mb'
      interval_label(5000)    returns '5kb'
      interval_label(500)     returns '500bp'
    :param bases <int>: window width, in bases
    :return <str>: SI-prefixed label, safe to embed in a filename
    """
    for multiplier, suffix in ((1000000000, "gb"), (1000000, "mb"), (1000, "kb")):
        if bases >= multiplier and bases % multiplier == 0:
            return f"{bases // multiplier}{suffix}"

    return f"{bases}bp"


def joint_option(prefix, valueslist):
    """Joins a list while adding a common prefix.
    Example:
      joint_option('-i', [1,2,3])
      '-i 1 -i 2 -i 3'
    """
    s = ""
    for v in valueslist:
        s += f"{prefix} {v} "
    s = s.rstrip()
    return s


# cfDNAanalyzer's feature extractors, mapped to the feature matrices each one
# contributes to a run's Features/ directory. cfDNAanalyzer writes one CSV per
# matrix, every one shaped samples x measurements with a leading sample,label
# column pair, and a feature can contribute more than one (EM produces both a
# motif frequency table and an MDS table, for instance). The mapping mirrors the
# feature_to_keys table in cfDNAanalyzer.sh, which is what actually decides the
# CSVs a run writes, so a change on either side has to be made on both. The
# order is the order its documentation lists the features in: genome-wide,
# region-specific, then TSS-based.
CDA_FEATURE_MATRICES = {
    # Genome-wide features
    "CNA": ("CNA",),
    "EM": ("EM_motifs_frequency", "EM_motifs_mds"),
    "FP": ("FP_fragmentation_profile",),
    # Region-specific features
    "NOF": ("NOF_meanfuziness", "NOF_occupancy"),
    "NP": ("NP_mean_coverage", "NP_central_coverage", "NP_amplitude"),
    "WPS": ("WPS_long", "WPS_short"),
    "OCF": ("OCF",),
    "EMR": (
        "EMR_aggregated_motif_frequency",
        "EMR_aggregated_mds",
        "EMR_region_motif_frequency",
        "EMR_region_mds",
    ),
    "FPR": ("FPR_fragmentation_profile_regions",),
    # Transcription start site features
    "PFE": ("PFE",),
    "TSSC": ("TSSC_average_coverage",),
}

# Features cfDNAanalyzer scores over a BED of regions handed to its -b option.
# It refuses to start at all when any of these is requested without one, so the
# region BED decides whether they can run rather than just what they measure.
CDA_REGION_FEATURES = ("NOF", "NP", "WPS", "OCF", "EMR", "FPR")

# Features that read the reference FastA (cfDNAanalyzer's -f option) to look up
# the sequence at fragment ends or, for NP, to GC-correct coverage.
CDA_FASTA_FEATURES = ("EM", "EMR", "NP")

# The genome builds cfDNAanalyzer accepts for its -g option. Unlike every other
# reference this pipeline hands it, the build is a name rather than a path: its
# CNA, FP, PFE and TSSC extractors read GC/mappability tracks, 100kb bin
# definitions and gene annotations that ship inside its own installation, indexed
# by this name, and its driver refuses to start on any other value. The same
# tuple is CDA_GENOME_BUILDS in src/utils.py, which validates the key on the
# command line.
CDA_GENOME_BUILDS = ("hg19", "hg38")


def cda_genome_build(config, genome_files):
    """The build name cfDNAanalyzer is run with, or None when it cannot be run
    at all. A bundled build is named for the assembly it is, so the selected
    genome is the answer; a custom build declares which of cfDNAanalyzer's two
    supported assemblies its coordinates match with a 'cda_genome' key, since its
    own bundled references are only available for those two.

    Returning None is what drops every cfDNAanalyzer feature from a run against
    a build cfDNAanalyzer has no references for, rather than having each of its
    jobs fail inside the container on a build name it rejects.
    :param config <dict>: the merged pipeline config (config.json)
    :param genome_files <dict>: reference files for the selected genome build
    :return <str or None>: 'hg19', 'hg38', or None if neither applies
    """
    build = genome_files.get("cda_genome") or config["options"]["genome"]

    return build if build in CDA_GENOME_BUILDS else None


def cda_regions(config, genome_files):
    """Resolves the BED whose regions cfDNAanalyzer's region-specific features
    are scored over. --cda-regions names one explicitly; without it the genome
    build's TSS interval BED stands in, since that is the region definition the
    rest of the pipeline already measures against. Returns None when neither is
    available, which is what drops the region-specific features from the run.
    :param config <dict>: the merged pipeline config (config.json)
    :param genome_files <dict>: reference files for the selected genome build
    :return <str or None>: path to the source region BED, if there is one
    """
    return config["options"].get("cda_regions") or genome_files.get("tss_interval")


def cda_analyses(config, genome_files):
    """The cfDNAanalyzer features a run can actually extract, i.e. the features
    --cda-features asked for minus any whose reference files the selected genome
    build does not supply. This mirrors how the rest of the pipeline gates its
    analyses: a feature that has nothing to run against is dropped rather than
    failing the run.

    A build cfDNAanalyzer has no bundled references for drops every feature, not
    just some, since the build name is an argument to its driver rather than
    something any one extractor reads.
    :param config <dict>: the merged pipeline config (config.json)
    :param genome_files <dict>: reference files for the selected genome build
    :return <list>: requested features this run has the references for
    """
    requested = config["options"].get("cda_features") or []
    if not cda_genome_build(config, genome_files):
        return []

    has_regions = bool(cda_regions(config, genome_files))
    has_fasta = bool(genome_files.get("reference_fa"))

    return [
        feature
        for feature in requested
        if (has_regions or feature not in CDA_REGION_FEATURES)
        and (has_fasta or feature not in CDA_FASTA_FEATURES)
    ]


def cda_matrices(features):
    """Lists the feature matrices cfDNAanalyzer writes for a set of features,
    in the catalogue's order so the targets are stable run to run.
    Example:
      cda_matrices(['EM', 'OCF']) returns
      ['EM_motifs_frequency', 'EM_motifs_mds', 'OCF']
    :param features <list>: cfDNAanalyzer feature names
    :return <list>: matrix names, each of which names a <matrix>.csv
    """
    matrices = []
    for feature, feature_matrices in CDA_FEATURE_MATRICES.items():
        if feature in features:
            matrices.extend(feature_matrices)

    return matrices
