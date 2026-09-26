#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
import hashlib
import json
import os
import re
import subprocess
import sys
from shutil import copytree, copyfileobj


class Colors():
    """Class encoding for ANSI escape sequeces for styling terminal text.
    Any string that is formatting with these styles must be terminated with
    the escape sequence, i.e. `Colors.end`.
    """
    # Escape sequence
    end = '\33[0m'
    # Formatting options
    bold   = '\33[1m'
    italic = '\33[3m'
    url    = '\33[4m'
    blink  = '\33[5m'
    higlighted = '\33[7m'
    # Text Colors
    black  = '\33[30m'
    red    = '\33[31m'
    green  = '\33[32m'
    yellow = '\33[33m'
    blue   = '\33[34m'
    pink  = '\33[35m'
    cyan  = '\33[96m'
    white = '\33[37m'
    # Background fill colors
    bg_black  = '\33[40m'
    bg_red    = '\33[41m'
    bg_green  = '\33[42m'
    bg_yellow = '\33[43m'
    bg_blue   = '\33[44m'
    bg_pink  = '\33[45m'
    bg_cyan  = '\33[46m'
    bg_white = '\33[47m'


def cat(files, output_file):
    """Concatenates mutiple files into one file. It operates similar
    to the cat unix command.
    @param files list[<str>]:
        List of files to concatenate together
    @param output_file <str>:
        Name of concatenated output file
    @return output_file <str>
    """
    # shutil.copyfileobj() automatically reads input
    # files chunk by chunk, which is much more memory
    # efficent when working with very large files
    with open(output_file,'wb') as ofh:
        for f in files:
            with open(f,'rb') as ifh:
                copyfileobj(ifh, ofh)
    return output_file


def md5sum(filename, first_block_only = False, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks/chunks defined by the blocksize parameter. This is
    a safer option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param first_block_only <bool>:
        Calculate md5 checksum of the first block/chunk only
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        if first_block_only:
            # Calculate MD5 of first block or chunck of file.
            # This is a useful heuristic for when potentially
            # calculating an MD5 checksum of thousand or
            # millions of file.
            hasher.update(buf)
            return hasher.hexdigest()
        while len(buf) > 0:
            # Calculate MD5 checksum of entire file
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


def permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the user is authorized to access
    a file/directory. Checks for existence, readability, writability and executability via:
    os.F_OK (tests existence), os.R_OK (tests read), os.W_OK (tests write), os.X_OK (tests exec).
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        Returns abs path if it exists and permissions are correct
    """
    if not exists(path):
        parser.error("Path '{0}' does not exists! Failed to provide vaild input.".format(path))
    if not os.access(path, *args, **kwargs):
        parser.error("Path '{0}' exists, but cannot read path due to permissions!".format(path))

    return os.path.abspath(path)


def quality_threshold(parser, value, *args, **kwargs):
    """Checks that a quality-score threshold is a non-negative integer. Used by
    the --mapscore and --baseqscore options, which are both '>=' thresholds; a
    negative threshold would silently keep everything rather than doing what
    the user asked, so it is rejected instead.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param value <str>:
        Value provided on the command line
    @return threshold <int>:
        The threshold as a non-negative integer
    """
    try:
        threshold = int(value)
    except ValueError:
        parser.error(
            "Quality threshold '{0}' is not an integer! Please provide a "
            "non-negative whole number, i.e. a Phred score.".format(value)
        )
    if threshold < 0:
        parser.error(
            "Quality threshold '{0}' is negative! Please provide a "
            "non-negative whole number, where 0 disables the filter.".format(value)
        )

    return threshold


# SI prefixes accepted by interval_size(), mapped to their multiplier in bases.
# Only whole-base multipliers make sense for a genomic coordinate, so the table
# stops at gigabases. A bare number (no prefix) is read as a base count.
SI_BASE_MULTIPLIERS = {
    '': 1,
    'b': 1,
    'bp': 1,
    'k': 1_000,
    'kb': 1_000,
    'm': 1_000_000,
    'mb': 1_000_000,
    'g': 1_000_000_000,
    'gb': 1_000_000_000,
}

# Splits an SI-prefixed base unit into its magnitude and its unit, e.g.
# '2MB' -> ('2', 'MB'), '500' -> ('500', ''). The magnitude is allowed to be
# fractional so '1.5mb' can be given, as long as it lands on a whole base.
SI_BASE_UNIT_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*([a-zA-Z]*)\s*$')


def interval_size(parser, value, *args, **kwargs):
    """Parses an SI-prefixed base unit into a whole number of bases. Used by the
    --interval option, which sets the width of the genomic windows the pipeline
    tiles the reference into. Accepts a bare base count or an SI prefix in either
    case, i.e. 2MB, 2mb, 2m, 2000kb and 2000000 are all two megabases.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param value <str>:
        Value provided on the command line
    @return bases <int>:
        The interval width as a positive whole number of bases
    """
    match = SI_BASE_UNIT_RE.match(str(value))
    if not match:
        parser.error(
            "Interval size '{0}' is not a base unit! Please provide a number "
            "with an optional SI prefix, i.e. 2MB for two megabases, 5kb for "
            "five kilobases, or 500 for five hundred bases.".format(value)
        )

    magnitude, unit = match.group(1), match.group(2).lower()
    if unit not in SI_BASE_MULTIPLIERS:
        parser.error(
            "Interval size '{0}' has an unrecognized unit '{1}'! Supported "
            "units are: {2} (or no unit at all, for a plain base count).".format(
                value, match.group(2),
                ', '.join(u for u in SI_BASE_MULTIPLIERS if u)
            )
        )

    bases = float(magnitude) * SI_BASE_MULTIPLIERS[unit]
    if bases != int(bases):
        # A window has to start and end on a base, so a fractional result is
        # the user asking for something that cannot be tiled, i.e. 1.5bp.
        parser.error(
            "Interval size '{0}' does not resolve to a whole number of bases "
            "({1})! Please provide a size that lands on a base boundary.".format(
                value, bases
            )
        )
    if int(bases) <= 0:
        parser.error(
            "Interval size '{0}' is not positive! Please provide a size "
            "greater than zero.".format(value)
        )

    return int(bases)


# Reference file keys a genome build may define, each mapped to a short note on
# what the pipeline loses without it. Every key is optional apart from the ones
# in GENOME_REQUIRED_KEYS below: the workflow gates each analysis on the files it
# needs being present for the selected build, so a build that omits a key simply
# does not run the analyses that read it. The same keys are read by
# workflow/rules/common.smk, so a key added here has to be read there too or it
# will validate on the command line and then be ignored by the workflow.
GENOME_REFERENCE_KEYS = {
    'chrom_sizes':   'delfi, adjust-wps and cleavage-profile',
    'ref2bit':       'end-motifs, mds, interval-end-motifs and delfi',
    'reference_fa':  'the genomic interval BED, and so frag-length-intervals, '
                     'interval-end-motifs and delfi; FastQ alignment',
    'bwamem2_index': 'FastQ alignment (bwa-mem2)',
    'dict':          'the reference contig filter, which is what checks BAM '
                     'input really was aligned to this build',
    'tss':           'wps and cleavage-profile',
    'tss_interval':  'coverage, adjust-wps and agg-bw',
    'blacklist':     "delfi's blacklist exclusion",
    'gap':           "delfi's structural gap exclusion",
    'cda_genome':    'every cfDNAanalyzer feature',
}

# Reference keys a build cannot leave out. Unlike the rest, the coverage
# analysis is an unconditional target of every run (the merged coverage workbook
# and the MultiQC report both gather over it), so a build with no TSS interval
# BED to quantify over has nothing to fall back to and would fail mid-run.
GENOME_REQUIRED_KEYS = ('tss_interval',)

# Keys whose value is not a path and so is neither resolved nor permission
# checked as one. cda_genome names which of cfDNAanalyzer's two supported builds
# its bundled references should be read from, for a custom build that is
# coordinate compatible with one of them.
GENOME_NON_PATH_KEYS = ('cda_genome',)

# The genome builds cfDNAanalyzer accepts for its -g option. Its driver rejects
# anything else outright, and several of its features index bundled per-build
# reference files by this name, so it is the one part of a custom build that
# cannot simply be pointed at a user file.
CDA_GENOME_BUILDS = ('hg19', 'hg38')

# Names a custom genome build may go by. The name is not cosmetic: it is what
# names the generated interval BED (intervals/{genome}_{size}_intervals.bed) and
# is recorded in config.json, so it has to be safe to embed in a filename.
GENOME_NAME_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._+-]*$')


def is_genome_config(value):
    """Decides whether a --genome value names a custom genome config file rather
    than one of the builds bundled in config/genome.json. A '.json' suffix is the
    whole test: build aliases are bare names ('hg38'), so the two can never be
    confused, and a value that was meant to be a file but is not named like one
    fails as an unrecognized alias with a message that points at this.
    @param value <str>:
        Value provided to --genome
    @return <bool>:
        True when the value should be read as a genome config file
    """
    return str(value).lower().endswith('.json')


def read_genome_config(path):
    """Reads a custom genome config file into the build name and reference file
    set the pipeline injects into config['references'].

    The file holds one genome build, in the shape an entry of
    config/genome.json has, so a build can be developed by copying an entry out
    of that file. Three spellings are accepted, because all three are what
    "one entry" plausibly means:

      1. the reference keys on their own, which is the canonical form:
             {"chrom_sizes": "...", "ref2bit": "...", ...}
      2. the entry with its build name, exactly as it appears in genome.json:
             {"mm10": {"chrom_sizes": "...", ...}}
      3. the whole genome.json shape, for a single build:
             {"references": {"mm10": {"chrom_sizes": "...", ...}}}

    Form 1 takes the build name from an optional "name" key, falling back to the
    config file's own basename; forms 2 and 3 take it from the key naming the
    build. Relative reference paths are resolved against the working directory,
    the same way every other path option of the frontend is, since the workflow
    runs from the output directory and cannot resolve them itself.

    Raises ValueError, rather than reporting the problem itself, so the same
    validation can be surfaced as an argparse error at parse time and as a fatal
    at config-build time.
    @param path <str>:
        Path to the genome config file
    @return (name, references) <tuple[str, dict]>:
        The build name and its reference files, with every path made absolute
    """
    try:
        with open(path) as fh:
            data = json.load(fh)
    except json.JSONDecodeError as e:
        raise ValueError(
            "Genome config '{0}' is not valid JSON: {1}".format(path, e)
        )
    except OSError as e:
        raise ValueError("Genome config '{0}' cannot be read: {1}".format(path, e))

    if not isinstance(data, dict) or not data:
        raise ValueError(
            "Genome config '{0}' is not a non-empty JSON object! Please provide "
            "the reference files of one genome build, in the shape of an entry "
            "of config/genome.json.".format(path)
        )

    # Unwrap forms 3 then 2 down to the reference keys themselves. A build's
    # reference values are all strings and a wrapper's are all objects, so the
    # nesting is what tells the forms apart rather than any declared version.
    if 'references' in data and isinstance(data['references'], dict):
        data = data['references']
    name = None
    if data and all(isinstance(value, dict) for value in data.values()):
        if len(data) != 1:
            raise ValueError(
                "Genome config '{0}' defines {1} genome builds ({2})! Please "
                "provide exactly one; --genome selects a single build.".format(
                    path, len(data), ', '.join(sorted(data))
                )
            )
        name, data = next(iter(data.items()))

    references = dict(data)
    # An explicit name wins over the one the file or its wrapper key implies,
    # so a config can be renamed or shared without changing what its output is
    # labelled with.
    name = references.pop('name', None) or name \
        or os.path.basename(path).rsplit('.', 1)[0]

    if not isinstance(name, str) or not GENOME_NAME_RE.match(name):
        raise ValueError(
            "Genome build name '{0}', from genome config '{1}', is not a valid "
            "name! It names the generated interval BED, so it must start with a "
            "letter or digit and hold only letters, digits, '.', '_', '+' or "
            "'-'. Set a \"name\" key in the config to choose one "
            "explicitly.".format(name, path)
        )

    if not references:
        raise ValueError(
            "Genome config '{0}' defines no reference files! Please provide at "
            "least: {1}.".format(path, ', '.join(GENOME_REQUIRED_KEYS))
        )

    unknown = [key for key in references if key not in GENOME_REFERENCE_KEYS]
    if unknown:
        raise ValueError(
            "Genome config '{0}' has unrecognized key(s) {1}! Supported keys "
            "are: {2}. An unrecognized key is rejected rather than ignored "
            "because a misspelled one would silently disable the analyses that "
            "read it.".format(
                path,
                ', '.join("'{0}'".format(key) for key in sorted(unknown)),
                ', '.join(sorted(GENOME_REFERENCE_KEYS))
            )
        )

    missing = [key for key in GENOME_REQUIRED_KEYS if not references.get(key)]
    if missing:
        raise ValueError(
            "Genome config '{0}' is missing required key(s) {1}! Every other "
            "reference file gates the analyses that read it and may be left "
            "out, but this one is read by an analysis every run performs "
            "({2}).".format(
                path,
                ', '.join("'{0}'".format(key) for key in missing),
                ', '.join(
                    GENOME_REFERENCE_KEYS[key] for key in missing
                )
            )
        )

    for key, value in references.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                "Genome config '{0}' key '{1}' is not a non-empty string! "
                "Please give it a value, or leave the key out entirely to skip "
                "the analyses that read it ({2}).".format(
                    path, key, GENOME_REFERENCE_KEYS[key]
                )
            )

    cda_genome = references.get('cda_genome')
    if cda_genome and cda_genome not in CDA_GENOME_BUILDS:
        raise ValueError(
            "Genome config '{0}' sets cda_genome to '{1}', which cfDNAanalyzer "
            "does not support! It must be one of: {2}. This key says which of "
            "cfDNAanalyzer's own per-build references to read, so it is only "
            "meaningful when your build is coordinate compatible with one of "
            "them; leave it out to skip cfDNAanalyzer.".format(
                path, cda_genome, ', '.join(CDA_GENOME_BUILDS)
            )
        )

    # Absolute paths from here on. resolve_additional_bind_paths() derives the
    # container bind points from these values and assumes they are absolute, so
    # a relative reference path would otherwise be left unbound and unreadable
    # inside the image.
    for key, value in references.items():
        if key in GENOME_NON_PATH_KEYS:
            continue
        references[key] = os.path.abspath(os.path.expanduser(value.strip()))

    unreadable = [
        "{0} ({1})".format(key, references[key])
        for key in references
        if key not in GENOME_NON_PATH_KEYS
        and not os.access(references[key], os.R_OK)
    ]
    if unreadable:
        raise ValueError(
            "Genome config '{0}' points at reference file(s) that do not exist "
            "or cannot be read:\n  {1}".format(path, '\n  '.join(unreadable))
        )

    return name, references


def genome_build(parser, value, repo_path, *args, **kwargs):
    """Checks a --genome value, which is either the alias of a build bundled in
    the pipeline's config/genome.json or the path to a custom genome config file
    (see read_genome_config()). The value is returned as given, i.e. resolution
    into the reference file set happens later, in src.run.setup(), against the
    genome.json that was copied into the output directory.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param value <str>:
        Value provided on the command line
    @param repo_path <str>:
        Path to the pipeline's installation, holding config/genome.json
    @return <str>:
        A bundled build alias, or an absolute path to a genome config file
    """
    if is_genome_config(value):
        path = permissions(parser, value, os.R_OK)
        try:
            read_genome_config(path)
        except ValueError as e:
            parser.error(str(e))
        return path

    bundled = os.path.join(repo_path, 'config', 'genome.json')
    try:
        with open(bundled) as fh:
            builds = sorted(json.load(fh).get('references', {}))
    except (OSError, json.JSONDecodeError):
        # Without the bundled config there is nothing to check the alias
        # against. Let it through rather than failing on the frontend's own
        # installation: src.run.setup() checks it again, against the copy in the
        # output directory, which is the one the workflow actually reads.
        return value

    if value not in builds:
        parser.error(
            "Genome build '{0}' is not one of the bundled builds ({1})! Please "
            "provide one of those, or the path to a genome config JSON file "
            "describing your own reference files.".format(
                value, ', '.join(builds)
            )
        )

    return value


# The feature extractors cfDNAanalyzer exposes through its -F option, grouped
# the way its documentation groups them. Every one of these is a paired-end
# analysis apart from CNA and TSSC, which is not a restriction worth encoding
# here because the pipeline only supports paired-end input to begin with. The
# same names are the keys of CDA_FEATURE_MATRICES in workflow/scripts/common.py,
# which maps each feature to the output matrices it produces, so a feature added
# here has to be added there too or the workflow will not know what it writes.
CDA_FEATURES = (
    # Genome-wide: copy number alteration, end motif, fragmentation profile
    'CNA', 'EM', 'FP',
    # Region-specific: nucleosome occupancy/fuzziness, nucleosome profile,
    # windowed protection score, orientation-aware fragmentation, regional end
    # motif, regional fragmentation profile
    'NOF', 'NP', 'WPS', 'OCF', 'EMR', 'FPR',
    # Transcription start site: promoter fragmentation entropy, TSS coverage
    'PFE', 'TSSC',
)


def cda_feature_list(parser, value, *args, **kwargs):
    """Parses the comma-separated feature list given to --cda-features into the
    cfDNAanalyzer feature names the workflow selects its rules with. Feature
    names are matched case-insensitively, and the two collective names are
    accepted as shorthand: 'all' selects every feature and 'none' selects none,
    which is the default and leaves cfDNAanalyzer out of the run entirely.

    The returned list is ordered and de-duplicated by CDA_FEATURES rather than
    by the order the features were typed in, so that every spelling of the same
    request produces the same pipeline targets.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param value <str>:
        Value provided on the command line, i.e. 'CNA,EM,OCF'
    @return features <list[str]>:
        The selected cfDNAanalyzer features, canonically named and ordered
    """
    requested = [feature.strip().upper() for feature in str(value).split(',')]
    requested = [feature for feature in requested if feature]

    if not requested or requested == ['NONE']:
        return []
    if requested == ['ALL']:
        return list(CDA_FEATURES)

    # 'all'/'none' are collective, so mixing either with a named feature is
    # ambiguous rather than additive, i.e. 'none,CNA' has no useful reading.
    # Checked before the unknown-name check below, which would otherwise claim
    # 'ALL' is an unrecognized feature and bury the actual mistake.
    if 'ALL' in requested or 'NONE' in requested:
        parser.error(
            "cfDNAanalyzer feature list '{0}' mixes 'all' or 'none' with named "
            "features! Please provide either one of those on its own, or a "
            "comma-separated list of feature names.".format(value)
        )

    unknown = [feature for feature in requested if feature not in CDA_FEATURES]
    if unknown:
        parser.error(
            "cfDNAanalyzer feature(s) {0} are not recognized! Please provide a "
            "comma-separated list of any of: {1}; or 'all' for every feature, "
            "or 'none' to skip cfDNAanalyzer.".format(
                ', '.join("'{0}'".format(feature) for feature in unknown),
                ', '.join(CDA_FEATURES)
            )
        )

    return [feature for feature in CDA_FEATURES if feature in requested]


def cda_bin_size(parser, value, *args, **kwargs):
    """Checks that a copy number alteration bin size is one cfDNAanalyzer can
    run. Its CNA extractor reads pre-computed GC and mappability correction
    tracks that ship with the ichorCNA it bundles, and those exist at four bin
    sizes only, so any other size fails inside the container rather than at the
    command line.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param value <str>:
        Value provided on the command line
    @return binsize <int>:
        The bin size, in kilobases
    """
    supported = (10, 50, 500, 1000)
    try:
        binsize = int(value)
    except ValueError:
        parser.error(
            "cfDNAanalyzer CNA bin size '{0}' is not an integer! Please "
            "provide one of: {1} (kilobases).".format(
                value, ', '.join(str(size) for size in supported)
            )
        )
    if binsize not in supported:
        parser.error(
            "cfDNAanalyzer CNA bin size '{0}' is not supported! The bundled "
            "GC and mappability tracks only exist for these bin sizes, in "
            "kilobases: {1}.".format(
                value, ', '.join(str(size) for size in supported)
            )
        )

    return binsize


def standard_input(parser, path, *args, **kwargs):
    """Checks for standard input when provided or permissions using permissions().
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        If path exists and user can read from location
    """
    # Checks for standard input
    if not sys.stdin.isatty():
        # Standard input provided, set path as an
        # empty string to prevent searching of '-'
        path = ''
        return path

    # Checks for positional arguments as paths
    path = permissions(parser, path, *args, **kwargs)

    return path


def tool_version(tool, cmd, strict = False):
    """Gets a tool's version using a known command.
    @param cmd list[<str>]:
        Command to run to get version of software.
    @param tool <str>:
        Name of software to check version.
    @param strict <bool>:
        If True, will exit with a fatal error if the software is not installed.
    @return version <str>:
        Version of the tool, default: ''. If strict is True, will exit with a fatal error.
    """
    # Get version information for a tool
    c = Colors
    version = ''
    try:
        # Merge standard error to standard output
        # some tools print version information to
        # stderr instead of stdout
        version = subprocess.check_output(cmd, stderr=subprocess.STDOUT).strip().decode('utf-8')
    except Exception as e:
        err("\n{0}{1}Warning: could not get version of {2} using: '{3}'{4}".format(
            c.bg_black, c.yellow, tool, ' '.join(cmd), c.end)
        )
        if strict:
            fatal(
                "{0}{1}Error: Please ensure {2} is installed and in $PATH.{3}".format(
                    c.bg_red, c.white, tool, c.end
                )
            )
    return version


def check_snakemake_version():
    """Checks the version of snakemake in the users $PATH.
    The pipeline supports snakemake versions less than 8.0.0. Version 8.0.0
    introduced a set of breaking changes that are not compatible with the
    current pipeline, so this is strictly enforced until we move to profiles.
    Fails fatally if an unsupported version is found.
    @return snakemake_version <str>:
        The full version string of the snakemake found in $PATH
    """
    snakemake_version = tool_version('snakemake', ['snakemake', '--version'], strict=True)

    parsed_version = re.search(
        r'^(?P<prefix>v)?(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)',
        snakemake_version.split()[-1]
    )

    if int(parsed_version.group('major')) >= 8:
        fatal(
            'Error: Snakemake version "{}" is not supported! '
            'Please use a version less than "8.0.0".'.format(snakemake_version)
        )

    return snakemake_version


def exists(testpath):
    """Checks if file exists on the local filesystem.
    @param parser <argparse.ArgumentParser() object>:
        argparse parser object
    @param testpath <str>:
        Name of file/directory to check
    @return does_exist <boolean>:
        True when file/directory exists, False when file/directory does not exist
    """
    does_exist = True
    if not os.path.exists(testpath):
        does_exist = False # File or directory does not exist on the filesystem

    return does_exist


def ln(files, outdir):
    """Creates symlinks for files to an output directory.
    @param files list[<str>]:
        List of filenames
    @param outdir <str>:
        Destination or output directory to create symlinks
    """
    # Create symlinks for each file in the output directory
    for file in files:
        ln = os.path.join(outdir, os.path.basename(file))
        if not exists(ln):
            os.symlink(os.path.abspath(os.path.realpath(file)), ln)


def which(cmd, path=None):
    """Checks if an executable is in $PATH
    @param cmd <str>:
        Name of executable to check
    @param path <list>:
        Optional list of PATHs to check [default: $PATH]
    @return <boolean>:
        True if exe in PATH, False if not in PATH
    """
    if path is None:
        path = os.environ["PATH"].split(os.pathsep)

    for prefix in path:
        filename = os.path.join(prefix, cmd)
        executable = os.access(filename, os.X_OK)
        is_not_directory = os.path.isfile(filename)
        if executable and is_not_directory:
            return True
    return False


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def require(cmds, suggestions, path=None):
    """Enforces an executable is in $PATH
    @param cmds list[<str>]:
        List of executable names to check
    @param suggestions list[<str>]:
        Name of module to suggest loading for a given index
        in param cmd.
    @param path list[<str>]]:
        Optional list of PATHs to check [default: $PATH]
    """
    error = False
    for i in range(len(cmds)):
        available = which(cmds[i])
        if not available:
            c = Colors
            error = True
            err("""\n{0}{1}Fatal: {2} is not in $PATH and is required during runtime!{3}
            └── Possible solution: please 'module load {4}' and run again!""".format(
                c.bg_red, c.white, cmds[i], c.end, suggestions[i])
            )

    if error: fatal()

    return


def safe_copy(source, target, resources = []):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            copytree(os.path.join(source, resource), destination)


def git_commit_hash(repo_path):
    """Gets the git commit hash of the repo.
    @param repo_path <str>:
        Path to git repo
    @return githash <str>:
        Latest git commit hash
    """
    try:
        githash = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT, cwd = repo_path
        ).strip().decode('utf-8')
        # Typecast to fix python3 TypeError (Object of type bytes is not JSON serializable)
        # subprocess.check_output() returns a byte string
        githash = str(githash)
    except Exception as e:
        # Github releases are missing the .git directory,
        # meaning you cannot get a commit hash, set the
        # commit hash to indicate its from a GH release
        githash = 'github_release'
    return githash


def join_jsons(templates):
    """Joins multiple JSON files to into one data structure
    Used to join multiple template JSON files to create a global config dictionary.
    @params templates <list[str]>:
        List of template JSON files to join together
    @return aggregated <dict>:
        Dictionary containing the contents of all the input JSON files
    """
    # Get absolute PATH to templates in git repo
    repo_path = os.path.dirname(os.path.abspath(__file__))
    aggregated = {}

    for file in templates:
        with open(os.path.join(repo_path, file), 'r') as fh:
            aggregated.update(json.load(fh))

    return aggregated


def check_cache(parser, cache, *args, **kwargs):
    """Check if provided SINGULARITY_CACHE is valid. Singularity caches cannot be
    shared across users (and must be owned by the user). Singularity strictly enforces
    0700 user permission on on the cache directory and will return a non-zero exitcode.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param cache <str>:
        Singularity cache directory
    @return cache <str>:
        If singularity cache dir is valid
    """
    c = Colors()
    if not exists(cache):
        # Cache directory does not exist on filesystem
        os.makedirs(cache)
    elif os.path.isfile(cache):
        # Cache directory exists as file, raise error
        parser.error("""\n\t{0}Fatal: Failed to provided a valid singularity cache!{1}
        The provided --singularity-cache already exists on the filesystem as a file.
        Please run {2} again with a different --singularity-cache location.
        """.format(c.red, c.end, sys.argv[0]))
    elif os.path.isdir(cache):
        # Provide cache exists as directory
        # Check that the user owns the child cache directory
        # May revert to os.getuid() if user id is not sufficent
        if exists(os.path.join(cache, 'cache')) and os.stat(os.path.join(cache, 'cache')).st_uid != os.getuid():
                # User does NOT own the cache directory, raise error
                parser.error("""\n\t{0}Fatal: Failed to provided a valid singularity cache!{1}
                The provided --singularity-cache already exists on the filesystem with a different owner.
                Singularity strictly enforces that the cache directory is not shared across users.
                Please run {0} again with a different --singularity-cache location.
                """.format(c.red, c.end, sys.argv[0]))

    return cache


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary
    """
    # Iterate over all values of given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over
            # all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type then
            # yield the value
            yield value


def hashed(l):
    """Returns an MD5 checksum for a list of strings. The list is sorted to
    ensure deterministic results prior to generating the MD5 checksum. This
    function can be used to generate a batch id from a list of input files.
    It is worth noting that path should be removed prior to calculating the
    checksum/hash.
    @Input:
        l list[<str>]: List of strings to hash
    @Output:
        h <str>: MD5 checksum of the sorted list of strings
    Example:
        $ echo -e '1\n2\n3' > tmp
        $ md5sum tmp
        # c0710d6b4f15dfa88f600b0e6b624077  tmp
        hashed([1,2,3])   # returns c0710d6b4f15dfa88f600b0e6b624077
    """
    # Sort list to ensure deterministic results
    l = sorted(l)
    # Convert everything to strings
    l = [str(s) for s in l]
    # Calculate an MD5 checksum of results
    h = hashlib.md5()
    # encode method ensure cross-compatiability
    # across python2 and python3
    h.update("{0}\n".format("\n".join(l)).encode())
    h = h.hexdigest()
    return h


if __name__ == '__main__':
    # Calculate MD5 checksum of entire file
    print('{0}  {1}'.format(md5sum(sys.argv[0]), sys.argv[0]))
    # Calcualte MD5 cehcksum of 512 byte chunck of file,
    # which is similar to following unix command:
    # dd if=utils.py bs=512 count=1 2>/dev/null | md5sum
    print('{0}  {1}'.format(md5sum(sys.argv[0], first_block_only = True, blocksize = 512), sys.argv[0]))
