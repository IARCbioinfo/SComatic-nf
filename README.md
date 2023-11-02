# SComatic-nf
## Empty template for nextflow pipelines (short description)

[![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/template-nf)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1404)
[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

![Workflow representation](template-nf.png)

## Description
Detects somatic single-nucleotide mutations in high-throughput single-cell genomics and transcriptomics data sets, such as single-cell RNA-seq and single-cell ATAC-seq.
SComatic runs sequentially in 4 steps:

Step 1: Splitting alignment file(bam) in cell type specific bams using precomputed cell type annotations. <br>
Step 2: Collecting base count information at each position of individual cell type.  <br>
Step 3: Merging base count matrices of all cell.  <br>
Step 4: Detection of somatic mutations. Consists of 2 steps:  <br>
	&emsp;Step4.1. Applies a set of hard filters and Beta binomial tests to discount sites affected by recurrent technical artefacts as somatic mutations.  <br>
	&emsp;Step4.2. Additional filters based on external datasets (RNA editing and Panel of Normals), and flags clustered mutations. High quality mutations are marked with the label PASS in the FILTER column of the output file.


## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. [SComatic](https://github.com/cortes-ciriano-lab/SComatic)
   
You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input
  | Type      | Description     |
  |-----------|---------------|
  | --bam    | BAM file to be analysed (*bai must be available in the same folder). |
  | --meta    | Metadata file mapping cell barcodes to cell type. |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --scomat_path    |            /Users/lipika/SComatic | Scomatic installation folder path |
| --ref    |            ref.fa | genome reference files (with index) |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --param3   |            xx | ...... |
| --param4    |            xx | ...... |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| --flag2    |      .... |


## Usage
  ```
  ...
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |


## Detailed description (optional section)
...

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/template-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | contrib1*    |            xx | Developer to contact for support (link to specific gitter chatroom) |
  | contrib2    |            xx | Developer |
  | contrib3    |            xx | Tester |

## References (optional)

## FAQ (optional)
