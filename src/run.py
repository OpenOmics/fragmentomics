#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import os, re, json, sys, subprocess

# Local imports
from utils import (git_commit_hash,
    join_jsons,
    fatal,
    which,
    exists,
    err)

from . import __pipeline__, version as __version__


def init(repo_path, output_path, links=[], required=['workflow', 'resources', 'config']):
    """
    Initialize the output directory. If user provides a output
    directory path that already exists on the filesystem as a file 
    (small chance of happening but possible), a OSError is raised. If the
    output directory PATH already EXISTS, it will not try to create the directory.

    @param repo_path <str>:
        Path to installation source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param links list[<str>]:
        List of files to symlink into output_path
    @param required list[<str>]:
        List of folder to copy over into output_path
    """
    if not exists(output_path):
        # Pipeline output directory
        # does not exist on filesystem
        os.makedirs(output_path)

    elif exists(output_path) and os.path.isfile(output_path):
        # Provided Path for pipeline 
        # output directory exists as file
        raise OSError("""\n\tFatal: Failed to create provided pipeline output directory!
        User provided --output PATH already exists on the filesystem as a file.
        Please run {} again with a different --output PATH.
        """.format(sys.argv[0])
        )

    # Copy over templates are other required resources
    copy_safe(source = repo_path, target = output_path, resources = required)

    # Create renamed symlinks for each rawdata 
    # file provided as input to the pipeline
    inputs = stage_work_files(links, output_path)

    return inputs


def copy_safe(source, target, resources = []):
    """
    Given a list paths it will recursively copy each to the
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
            copytree(os.path.join(source, resource), destination)
    
    return


def stage_work_files(input_data, output_dir):
    """
    TODO

    @param input_data <list[<str>]>:
        List of input files to symlink to target location
    @param target <str>:
        Target path to copy templates and required resources
    @return input_files list[<str>]:
        List of renamed input file
    """
    supported_files = [
        'bam',
        'cram',
        'sam'
    ]

    # validate
    invalid_files = []
    for _file in input_data:
        check = False
        for _ext in supported_files:
            if _file.lower().endswith(f'.{_ext.lower()}'):
                check = True
                break
        if not check:
            invalid_files.append(os.path.basename(_file))

    if invalid_files:
        raise ValueError(
            f"The pipeline {__pipeline__} only accepts {'/'.join(supported_files)} file types!" + \
            f"\t> These files are invalid: {', '.join(invalid_files)}"
        )

    # sym-link
    linked_files = link_files_to_output(input_data, output_dir)
    
    return linked_files


def santitize_fn(unclean_fn):
    """
    TODO
    """
    allowed_chars = re.compile(r"[^a-zA-Z0-9_.-]")
    sanitized = allowed_chars.sub("_", unclean_fn)
    sanitized = sanitized.strip()
    if sanitized.startswith('.'):
        sanitized = '_' + sanitized[1:]
    if sanitized.endswith('.'):
        sanitized = sanitized[:-1] + '_'
    if '.sorted' in sanitized:
        sanitized = sanitized.replace('.sorted', '')
    return sanitized


def link_files_to_output(files, output_dir):
    """
    TODO
    """
    linked_files = []
    for file in files:
        filename = os.path.basename(file)
        renamed = os.path.join(output_dir, santitize_fn(filename))
        if not exists(renamed):
            os.symlink(os.path.abspath(os.path.realpath(file)), renamed)
        linked_files.append(renamed)
    return linked_files


def setup(sub_args, ifiles, repo_path, output_path, config_extra=None):
    """Setup the pipeline for execution and creates config file from templates
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param repo_path <str>:
        Path to installation or source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param config_extra <str>:
        Function for injecting additional processing on to `config`
    @return config <dict>:
         Config dictionary containing metadata to run the pipeline
    """
    # Check for mixed inputs,
    # inputs which are a mixture
    # of FastQ and BAM files 
    mixed_inputs(ifiles)

    # Resolves PATH to reference file 
    # template or a user generated 
    # reference genome built via build 
    # subcommand
    genome_config = os.path.join(output_path,'config','genome.json')
    # if sub_args.genome.endswith('.json'):
        # Provided a custom reference genome generated by build pipline
        # genome_config = os.path.abspath(sub_args.genome)

    required = {
        # Base configuration file
        "base": os.path.join(output_path,'config','config.json'),
        # Template for project-level information
        "project": os.path.join(output_path,'config','containers.json'),
        # Template for genomic reference files
        # User provided argument --genome is used to select the template
        "genome": genome_config,
        # Template for tool information
        "tools": os.path.join(output_path,'config', 'modules.json'),
    }

    # Create the global or master config 
    # file for pipeline, config.json
    config = join_jsons(required.values()) # uses templates in config/*.json 
    config['cluster'] = json.load(open(os.path.join(output_path, 'config', 'cluster.json')))
    config['project'] = {}
    config = add_user_information(config)
    config = add_rawdata_information(sub_args, config, ifiles)

    # Resolves if an image needs to be pulled 
    # from an OCI registry or a local SIF exists
    config = image_cache(sub_args, config, output_path)

    # Add other runtime info for debugging
    config['project']['version'] = __version__
    config['project']['workpath'] = os.path.abspath(sub_args.output)
    git_hash = git_commit_hash(repo_path)
    config['project']['git_commit_hash'] = git_hash   # Add latest git commit hash
    config['project']['pipeline_path'] = repo_path    # Add path to installation

    # Add all cli options for data provenance
    for opt, v in vars(sub_args).items():
        if opt == 'func':
            # Pass over sub command's handler
            continue
        elif not isinstance(v, (list, dict)):
            # CLI value can be converted to a string
            v = str(v)
        config['options'][opt] = v

    if config_extra:
        if not callable(config_extra):
            raise ValueError('`config_extra` needs to be a callable function')
        config = config_extra(config)

    return config


def get_fastq_screen_paths(fastq_screen_confs, match = 'DATABASE', file_index = -1):
    """Parses fastq_screen.conf files to get the paths of each fastq_screen database.
    This path contains bowtie2 indices for reference genome to screen against.
    The paths are added as singularity bind points.
    @param fastq_screen_confs list[<str>]:
        Name of fastq_screen config files to parse
    @param match <string>:
        Keyword to indicate a line match [default: 'DATABASE']
    @param file_index <int>:
        Index of line line containing the fastq_screen database path
    @return list[<str>]:
        Returns a list of fastq_screen database paths
    """
    databases = []
    for file in fastq_screen_confs:
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith(match):
                        db_path = line.strip().split()[file_index]
                        databases.append(db_path)
    return databases


def bind(sub_args, config):
    """Resolves bindpaths for singularity/docker images.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param configfile dict[<any>]:
        Config dictionary generated by setup command.
    @return bindpaths list[<str>]:
        List of singularity/docker bind paths 
    """
    bindpaths = []
    # Bind input file paths, working
    # directory, and other reference
    # genome paths
    rawdata_bind_paths = config['project']['datapath'].split(',')
    working_directory = config['project']['workpath']
    genome_resources = os.path.commonpath(config['references'][config['options']['genome_alias']].values())

    bin_directory = config['binpath']
    bindpaths = [working_directory, bin_directory, genome_resources] + rawdata_bind_paths

    corrected_bind_paths = []
    targets = []
    for _path in bindpaths:
        if _path == os.sep:
            continue
        elif ':' in _path:
            this_target = _path.split(':')[1]
            if this_target in targets:
                raise ValueError('Conflicting singularity binds!')
            targets.append(this_target)
            corrected_bind_paths.append(_path)
            continue
        if os.path.abspath(os.path.realpath(_path)) != _path:
            bind_str = f'{os.path.abspath(os.path.realpath(_path))}:{_path}:rw'
            if _path in targets:
                raise ValueError('Conflicting singularity binds!')
            corrected_bind_paths.append(bind_str)

    return corrected_bind_paths


def mixed_inputs(ifiles):
    """Check if a user has provided a set of input files which contain a
    mixture of FastQ and BAM files. The pipeline does not support processing
    a mix of FastQ and BAM files.
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    """
    bam_files, fq_files = [], []
    fastqs = False
    bams = False
    for file in ifiles:
        if file.endswith('.R1.fastq.gz') or file.endswith('.R2.fastq.gz'):
            fastqs = True 
            fq_files.append(file)
        elif file.endswith('.bam'):
            bams = True
            bam_files.append(file)

    if fastqs and bams:
        # User provided a mix of FastQs and BAMs
        raise TypeError("""\n\tFatal: Detected a mixture of --input data types. 
            A mixture of BAM and FastQ files were provided; however, the pipeline
            does NOT support processing a mixture of input FastQ and BAM files.
            Input FastQ Files:
                {}
            Input BAM Files:
                {}        
            Please do not run the pipeline with a mixture of FastQ and BAM files.
            This feature is currently not supported within '{}', and it is not
            recommended to process samples in this way either. If this is a priority
            for your project, please run the set of FastQ and BAM files separately 
            (in two separate output directories). If you feel like this functionality
            should exist, feel free to open an issue on Github.
            """.format(" ".join(fq_files), " ".join(bam_files), sys.argv[0])
        )

def add_user_information(config):
    """Adds username and user's home directory to config.
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    # Get PATH to user's home directory
    # Method is portable across unix-like 
    # OS and Windows
    home = os.path.expanduser("~")

    # Get username from home directory PATH
    username = os.path.split(home)[-1]

    # Update config with home directory and 
    # username
    config['project']['userhome'] = home
    config['project']['username'] = username

    return config


def add_sample_metadata(input_files, config, group=None):
    """Adds sample metadata such as sample basename, label, and group information.
    If sample sheet is provided, it will default to using information in that file.
    If no sample sheet is provided, it will only add sample basenames and labels.
    @params input_files list[<str>]:
        List containing pipeline input fastq files
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @params group <str>:
        Sample sheet containing basename, group, and label for each sample
    @return config <dict>:
        Updated config with basenames, labels, and groups (if provided)
    """
    config['samples'] = []
    for file in input_files:
        # Split sample name on file extension
        sample = re.split('(?i)\\.(?:sorted\\.)?(?:bam|sam|cram)', os.path.basename(file))[0]
        if sample not in config['samples']:
            config['samples'].append(sample)

    return config


def add_rawdata_information(sub_args, config, ifiles):
    """Adds information about rawdata provided to pipeline.
    Determines whether the dataset is paired-end or single-end and finds the set of all
    rawdata directories (needed for -B option when running singularity). If a user provides
    paired-end data, checks to see if both mates (R1 and R2) are present for each sample.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @params ifiles list[<str>]:
        List containing pipeline input files (renamed symlinks)
    @params config <dict>:
        Config dictionary containing metadata to run pipeline
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    # Finds the set of rawdata directories to bind
    rawdata_paths = get_rawdata_bind_paths(input_files = sub_args.input)
    config['project']['datapath'] = ','.join(rawdata_paths)

    # Add each sample's basename
    config = add_sample_metadata(input_files = ifiles, config = config)

    return config


def image_cache(sub_args, config, repo_path):
    """Adds Docker Image URIs, or SIF paths to config if singularity cache option is provided.
    If singularity cache option is provided and a local SIF does not exist, a warning is
    displayed and the image will be pulled from URI in 'config/containers.json'.
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @params config <file>:
        Docker Image config file
    @param repo_path <str>:
        Path to installation or source code and its templates
    @return config <dict>:
         Updated config dictionary containing user information (username and home directory)
    """
    images = os.path.join(repo_path, 'config','containers.json')

    # Read in config for docker image uris 
    with open(images, 'r') as fh:
        data = json.load(fh)
    # Check if local sif exists 
    for image, uri in data['images'].items():
        if sub_args.sif_cache:
            sif = os.path.join(sub_args.sif_cache, '{}.sif'.format(os.path.basename(uri).replace(':', '_')))
            if not exists(sif):
                # If local sif does not exist on in cache, 
                # print warning and default to pulling from 
                # URI in config/containers.json
                print('Warning: Local image "{}" does not exist in singularity cache'.format(sif), file=sys.stderr)
            else:
                # Change pointer to image from Registry URI 
                # to local SIF
                data['images'][image] = sif

    config.update(data)

    return config


def get_rawdata_bind_paths(input_files):
    """
    Gets rawdata bind paths of user provided fastq files.
    @params input_files list[<str>]:
        List containing user-provided input fastq files
    @return bindpaths <set>:
        Set of rawdata bind paths
    """
    bindpaths = []
    for file in input_files:
        # Get directory of input file
        rawdata_src_path = os.path.dirname(os.path.abspath(file))
        if rawdata_src_path not in bindpaths:
            bindpaths.append(rawdata_src_path)

    return bindpaths


def dryrun(outdir, config='config.json', snakefile=os.path.join('workflow', 'Snakefile')):
    """Dryruns the pipeline to ensure there are no errors prior to runnning.
    @param outdir <str>:
        Pipeline output PATH
    @return dryrun_output <str>:
        Byte string representation of dryrun command
    """
    try:
        # Setting cores to dummy high number so
        # displays the true number of cores a rule
        # will use, it uses the min(--cores CORES, N)
        dryrun_output = subprocess.check_output([
            'snakemake', '-npr',
            '-s', str(snakefile),
            '--use-singularity',
            '--rerun-incomplete',
            '--cores', str(256),
            '--configfile={}'.format(config)
        ], cwd = outdir,
        stderr=subprocess.STDOUT)
    except OSError as e:
        # Catch: OSError: [Errno 2] No such file or directory
        #  Occurs when command returns a non-zero exit-code
        if e.errno == 2 and not which('snakemake'):
            # Failure caused because snakemake is NOT in $PATH
            err('\n\x1b[6;37;41mError: Are snakemake AND singularity in your $PATH?\x1b[0m')
            fatal('\x1b[6;37;41mPlease check before proceeding again!\x1b[0m')
        else:
            # Failure caused by unknown cause, raise error
            raise e
    except subprocess.CalledProcessError as e:
        print(e, e.output.decode("utf-8"))
        raise(e)

    return dryrun_output


def runner(mode, outdir, alt_cache, logger, additional_bind_paths = None, 
    threads=2,  jobname='pl:master', submission_script='run.sh',
    tmp_dir = '/lscratch/$SLURM_JOBID/'):
    """Runs the pipeline via selected executor: local or slurm.
    If 'local' is selected, the pipeline is executed locally on a compute node/instance.
    If 'slurm' is selected, jobs will be submited to the cluster using SLURM job scheduler.
    Support for additional job schedulers (i.e. PBS, SGE, LSF) may be added in the future.
    @param outdir <str>:
        Pipeline output PATH
    @param mode <str>:
        Execution method or mode:
            local runs serially a compute instance without submitting to the cluster.
            slurm will submit jobs to the cluster using the SLURM job scheduler.
    @param additional_bind_paths <str>:
        Additional paths to bind to container filesystem (i.e. input file paths)
    @param alt_cache <str>:
        Alternative singularity cache location
    @param logger <file-handle>:
        An open file handle for writing
    @param threads <str>:
        Number of threads to use for local execution method
    @param masterjob <str>:
        Name of the master job
    @return masterjob <subprocess.Popen() object>:
    """
    # Add additional singularity bind PATHs
    # to mount the local filesystem to the 
    # containers filesystem, NOTE: these 
    # PATHs must be an absolute PATHs
    outdir = os.path.abspath(outdir)
    # Add any default PATHs to bind to 
    # the container's filesystem, like 
    # tmp directories, /lscratch
    addpaths = []
    temp = os.path.dirname(tmp_dir.rstrip('/'))
    if temp == os.sep:
        temp = tmp_dir.rstrip('/')
    # if outdir not in additional_bind_paths.split(','):
    #     addpaths.append(outdir)
    if temp not in additional_bind_paths.split(','):
        addpaths.append(temp)
    bindpaths = ','.join(addpaths)

    # Set ENV variable 'SINGULARITY_CACHEDIR' 
    # to output directory
    my_env = {}; my_env.update(os.environ)
    cache = os.path.join(outdir, ".singularity")
    my_env['SINGULARITY_CACHEDIR'] = cache
    if alt_cache:
        # Override the pipeline's default 
        # cache location
        my_env['SINGULARITY_CACHEDIR'] = alt_cache
        cache = alt_cache

    if additional_bind_paths:
        # Add Bind PATHs for outdir and tmp dir
        if bindpaths:
            bindpaths = ",{}".format(bindpaths)
        bindpaths = "{}{}".format(additional_bind_paths,bindpaths)

    if not exists(os.path.join(outdir, 'logfiles')):
        # Create directory for logfiles
        os.makedirs(os.path.join(outdir, 'logfiles'))

    # Create .singularity directory for 
    # installations of snakemake without
    # setuid which creates a sandbox in
    # the SINGULARITY_CACHEDIR
    if not exists(cache):
        # Create directory for sandbox 
        # and image layers
        os.makedirs(cache)

    # Run on compute node or instance
    # without submitting jobs to a scheduler
    if mode == 'local':
        # Run pipeline's main process
        # Look into later: it maybe worth 
        # replacing Popen subprocess with a direct
        # snakemake API call: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
        masterjob = subprocess.Popen([
                'snakemake', '-pr', '--rerun-incomplete',
                '--use-singularity',
                '--singularity-args', "'-B {}'".format(bindpaths),
                '--cores', str(threads),
                '--configfile=config.json'
            ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)

    # Submitting jobs to cluster via SLURM's job scheduler
    elif mode == 'slurm':
        # Run pipeline's main process
        # Look into later: it maybe worth 
        # replacing Popen subprocess with a direct
        # snakemake API call: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
        # CLUSTER_OPTS="'sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} \
        #   -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURMDIR/slurm-%j_{params.rname}.out \
        #   -o $SLURMDIR/slurm-%j_{params.rname}.out'"
        # sbatch --parsable -J "$2" --gres=lscratch:500 --time=5-00:00:00 --mail-type=BEGIN,END,FAIL \
        #   --cpus-per-task=32 --mem=96g --output "$3"/logfiles/snakemake.log --error "$3"/logfiles/snakemake.log \
        # snakemake --latency-wait 120 -s "$3"/workflow/Snakefile -d "$3" \
        #   --use-singularity --singularity-args "'-B $4'" --configfile="$3"/config.json \
        #   --printshellcmds --cluster-config "$3"/resources/cluster.json \
        #   --cluster "${CLUSTER_OPTS}" --keep-going --restart-times 3 -j 500 \
        #   --rerun-incomplete --stats "$3"/logfiles/runtime_statistics.json \
        #   --keep-remote --local-cores 30 2>&1 | tee -a "$3"/logfiles/master.log
        masterjob = subprocess.Popen([
                str(submission_script), mode,
                '-j', jobname, '-b', str(bindpaths),
                '-o', str(outdir), '-c', str(cache),
                '-t', "'{}'".format(tmp_dir)
            ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)

    return masterjob
