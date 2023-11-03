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

	&nbsp; Step4.1. Applies a set of hard filters and Beta binomial tests to discount sites affected by recurrent technical artefacts as somatic mutations.  <br>
 
	&nbsp; Step4.2. Additional filters based on external datasets (RNA editing and Panel of Normals), and flags clustered mutations. High quality mutations are marked with the label PASS in the FILTER column of the output file.


## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. [SComatic](https://github.com/cortes-ciriano-lab/SComatic)
   
You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input
  | Type      | Description     |
  |-----------|---------------|
  | --bam_folder    | Folder containing BAM files (*bai must be available in the same folder). |
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
| --cpu   |            2 | Number of CPUs |
| --mem    |            20 | memory |
| --output_folder    |            SComatic-nf-results | Output folder |
| --nTrim    |            5 | Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read |
| --maxNM   |            5 | Maximum number of mismatches permitted to consider reads for analysis |
| --maxNH    |            1 | Maximum number of alignment hits permitted to consider reads for analysis |
| --chrom    |            all | Chromosome to be analysed |
| --minbq    |            30 | Minimum base quality permited for the base counts |
| --nprocs    |            1 | Number of processes |
| --pon    |            30 | Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |



## Usage
 To use SComatic on your bamFile.bam, having metadata.tsv file mapping cell barcodes to cell type and reference genome ref.fa used in alignment, use this command
  ```
  nextflow run iarcbioinfo/SComatic-nf --bam bamFile.bam --meta metadata.tsv --ref ref.fa --scomat_path path/to/Scomatfolder
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | SplitBamCellTypes/sample.*.bam    | Folder containing cell-type-specific BAM files (step1 output)  |
  | Step2_BaseCellCounts/sample.*.tsv    | Folder containing base count information for each cell type and for every position in the genome (step2 output) |
  | Step3_BaseCellCountsMerged/sample.BaseCellCounts.AllCellTypes.tsv    | Folder containing merged base count file of all cell types. (step3 output)   |
  | Step4_VariantCalling/sample.calling.step*.tsv  (*=1,2)  | Folder containing two files files (1*.tsv: SNV called after applying filters for removing technical artefacts, 2*.tsv: Further filtered for RNA editing and PoN). (step4 output)   |
  
  


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
