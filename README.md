<div align="center">
   
  <h1>fragmentomics 🔬</h1>
  
  **_A pipeline for analyzing cell-free DNA profiles_**

  [![tests](https://github.com/OpenOmics/fragmentomics/workflows/tests/badge.svg)](https://github.com/OpenOmics/fragmentomics/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/fragmentomics/workflows/docs/badge.svg)](https://github.com/OpenOmics/fragmentomics/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/fragmentomics?color=brightgreen)](https://github.com/OpenOmics/fragmentomics/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/fragmentomics)](https://github.com/OpenOmics/fragmentomics/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, fragmentomics. Its long-term goals: to provide a standardized, reproducible, and scalable framework for cfDNA fragmentomics analysis, enabling comprehensive characterization of fragmentation patterns and supporting biomarker discovery across research and clinical applications.!
  </i>
</div>

## Overview

Welcome to fragmentomics! Before getting started, we highly recommend reading through [fragmentomics's documentation](https://openomics.github.io/fragmentomics/).

The **`./fragmentomics`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>fragmentomics <b>run</b></code>](https://openomics.github.io/fragmentomics/usage/run/): Run the fragmentomics pipeline with your input files.
 * [<code>fragmentomics <b>unlock</b></code>](https://openomics.github.io/fragmentomics/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>fragmentomics <b>install</b></code>](https://openomics.github.io/fragmentomics/usage/install/): Download reference files locally.
 * [<code>fragmentomics <b>cache</b></code>](https://openomics.github.io/fragmentomics/usage/cache/): Cache software containers locally.

**fragmentomics** is a comprehensive workflow for analyzing fragmentation patterns in cell-free DNA sequencing data. It runs [FinaleToolkit<sup>1</sup>](https://github.com/epifluidlab/FinaleToolkit), a python-based fragmentomics package, to process sorted and indexed BAM files generated from whole-genome sequencing, whole-genome bisulfite sequencing, or ChIP-seq experiments. The pipeline characterizes multiple cfDNA fragmentation features, including fragment-length distributions, genomic coverage, fragment-end motifs, motif diversity scores, window protection scores, GC-corrected DELFI short-to-long fragment ratios, and cleavage profiles. Together, these measurements provide a detailed view of cfDNA fragmentation patterns and generate standardized outputs that can support downstream biomarker discovery, disease classification, and other research applications. It relies on technologies like [Singularity<sup>2</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>3</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

As input, it accepts a set of paired-end Illumina FastQ or BAM files (FastQ and BAM cannot be mixed in a single run) and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users. Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/fragmentomics/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/fragmentomics/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/fragmentomics/issues).

## Dependencies

**Requires:** `singularity>=3.5`  `snakemake<=7.32.4`

At the current moment, the pipeline uses a mixture of enviroment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation

Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/fragmentomics.git
# Change your working directory
cd fragmentomics/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./fragmentomics -h
```

## Contribute

This site is a living document, created for and by members like you. fragmentomics is maintained by the members of [OpenOmics](https://openomics.github.io) and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/fragmentomics).


## Cite

If you use this software, please cite it as below:  

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
Citation coming soon!
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Citation coming soon!
```

</details>

## References

<sup>**1.** James Wenhan Li, Ravi Bandaru, Kundan Baliga, Yaping Liu. FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit. Bioinformatics Advances, 2025, vbaf236.</sup>  
<sup>**2.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**3.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  