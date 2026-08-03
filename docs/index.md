<div align="center">

  <h1 style="font-size: 250%">fragmentomics 🔬</h1>

  <b><i>A pipeline for analyzing cell-free DNA profiles</i></b><br> 
  <a href="https://github.com/OpenOmics/fragmentomics/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/fragmentomics/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/fragmentomics/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/fragmentomics/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/fragmentomics/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/fragmentomics?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/fragmentomics/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/fragmentomics">
  </a>

  <p>
    This is the home of the pipeline, fragmentomics. Its long-term goals: to provide a standardized, reproducible, and scalable framework for cfDNA fragmentomics analysis, enabling comprehensive characterization of fragmentation patterns and supporting biomarker discovery across research and clinical applications.!
  </p>

</div>  


## Overview

Welcome to fragmentomics's documentation! This guide is the main source of documentation for users that are getting started with the [long pipeline name](https://github.com/OpenOmics/fragmentomics/). 

The **`./fragmentomics`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">fragmentomics <b>run</b></code>](usage/run.md)   
    Run the fragmentomics pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">fragmentomics <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">fragmentomics <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">fragmentomics <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**fragmentomics** is a comprehensive workflow for analyzing fragmentation patterns in cell-free DNA sequencing data. It runs [FinaleToolkit<sup>1</sup>](https://github.com/epifluidlab/FinaleToolkit), a python-based fragmentomics package, to process sorted and indexed BAM files generated from whole-genome sequencing, whole-genome bisulfite sequencing, or ChIP-seq experiments. The pipeline characterizes multiple cfDNA fragmentation features, including fragment-length distributions, genomic coverage, fragment-end motifs, motif diversity scores, window protection scores, GC-corrected DELFI short-to-long fragment ratios, and cleavage profiles. Together, these measurements provide a detailed view of cfDNA fragmentation patterns and generate standardized outputs that can support downstream biomarker discovery, disease classification, and other research applications.. It relies on technologies like [Singularity<sup>2</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>3</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

As input, it accepts a set of paired-end Illumina FastQ or BAM files (FastQ and BAM cannot be mixed in a single run) and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users. Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/fragmentomics/issues).

## Contribute 

This site is a living document, created for and by members like you. fragmentomics is maintained by the members of [OpenOmics](https://openomics.github.io/) and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/fragmentomics).

## Citation 

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    Citation coming soon!
    ```

=== "APA"

    ```
    Citation coming soon!
    ```

## References

<sup>**1.** James Wenhan Li, Ravi Bandaru, Kundan Baliga, Yaping Liu. FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit. Bioinformatics Advances, 2025, vbaf236.</sup>  
<sup>**2.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**3.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
