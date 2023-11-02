#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (C) 2023 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "------------------------------------------------------------"
log.info "  Nextflow Pipeline: Calling somatic variants with SComatic "
log.info "------------------------------------------------------------"
log.info "Copyright (C) 2023 IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run SComatic-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--scomat_path   <DIR>   SComatic install folder"
//    log.info "--input_folder2   <DIR>   Folder containing other VCF files"
    log.info ""
    log.info "Optional arguments:"
    log.info "--output_folder      <DIR>   Output folder"
//    log.info "--vcfSuffix1_snvs   <STRING>  Suffix (without the extension) of files in"
//    log.info "                              input_folder1 containing SNVs"
//    log.info "--vcfSuffix1_indels <STRING>  Suffix (without the extension) of files in"
//    log.info "                              input_folder1 containing indels and MNVs"
//    log.info "--vcfSuffix2_snvs   <STRING>  Suffix (without the extension) of files in"
//    log.info "                              input_folder2 containing SNVs"
//    log.info "--vcfSuffix2_indels <STRING>  Suffix (without the extension) of files in"
//    log.info "                              input_folder2 containing indels and MNVs"
//    log.info "--ext               <STRING>  Extension of variant calling files"

// Cluster options
params.mem = 2
params.cpu = 2
exit 0
} else {
    /* Software information */
    log.info "help:               ${params.help}"
    log.info "scomat_path:    ${params.scomat_path}"
//    log.info "input_folder2:   ${params.input_folder2}"
    log.info "output_folder:       ${params.output_folder}"
}


// Define the pipeline parameters

// results directory- step 0 // also SComatic path
params.output_folder="SComatic"

//common parameters
params.sample="Example"

//step1 paramters
//params.scomat_path="SComatic"
params.bam="Example.scrnaseq.bam"
params.bai="Example.scrnaseq.bam.bai"
params.meta="Example.cell_barcode_annotations.tsv"
params.nTrim=5
params.maxNM=5
params.maxNH=1

//step2 parameters
//params.output_step_dir2="SComatic/results/Step2_BaseCellCounts"
params.ref="chr10.fa"
params.chrom='all'
params.minbq=30
params.nprocs=1

//step4 parameters
//params.output_step_dir3="SComatic/results/Step3_BaseCellCountsMerged"
//params.output_step_dir4="SComatic/results/Step4_VariantCalling"
params.editing="SComatic/RNAediting/AllEditingSites.hg38.txt"
params.pon="SComatic/PoNs/PoN.scRNAseq.hg38.tsv"


process step1 {
	maxForks 1
	
    input:
    path scomat_folder
    tuple val(ID), path(files)
    path meta
	val nTrim
	val maxNM
	val maxNH
	//path output_step1

	output:
	path("Step1_BamCellTypes")

	publishDir "${params.output_folder}", mode: 'copy'

	script:
	"""
	mkdir -p Step1_BamCellTypes
	python ${scomat_folder}/scripts/SplitBam/SplitBamCellTypes.py --bam ${files[0]} \
		--meta ${meta} \
		--n_trim ${nTrim} \
        	--max_nM ${maxNM} \
        	--max_NH ${maxNH} \
        	--outdir Step1_BamCellTypes
    	"""
}


process step2 {
	maxForks 1
	
    input:
    //	path output_dir2
	path y
	path output_step1
	path ref
	val chrom
	val minbq
	val nprocs
	//val ready

	output: 
    path("Step2_BaseCellCounts")

	publishDir "${params.output_folder}", mode: 'copy'

    	script:
    	"""
	echo 'Hello!'
	for bam in \$(ls -d $output_step1/*bam);do
		echo 'World!'
  		echo \$bam

 		cell_type=\$(basename \$bam | awk -F '.' '{print \$(NF-1)}')
        	echo \$cell_type

  		# Temp folder
  		temp=\$(echo Step2_BaseCellCounts/temp_\${cell_type})
		echo \$temp
  		mkdir -p \$temp

        	# Command line to submit to cluster
        	python ${y}/scripts/BaseCellCounter/BaseCellCounter.py --bam \$bam \
            	--ref $ref \
            	--chrom $chrom \
            	--out_folder Step2_BaseCellCounts \
            	--min_bq $minbq\
            	--tmp_dir \${temp} \
            	--nprocs $nprocs
		
		rm -rf \${temp}
        	
    	done
    	"""
}

process step3 {
	maxForks 1
	
	input:
	path y
	path output_step2
	val sample
	//val ready

	output: 
    path("Step3_BaseCellCountsMerged")

	publishDir "${params.output_folder}", mode: 'copy'

	script:
	"""
	mkdir -p Step3_BaseCellCountsMerged

	python ${y}/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder $output_step2 \
	--outfile Step3_BaseCellCountsMerged/${sample}.BaseCellCounts.AllCellTypes.tsv
	"""

}

process step4 {
	
	input:
	path y
	path output_step3
	val sample
	//path output_dir4
	path ref
	path editing
	path pon
	//val ready


	output: 
    path("Step4_VariantCalling")

	publishDir "${params.output_folder}", mode: 'copy'

	script:
	"""
	mkdir -p Step4_VariantCalling

	python ${y}/scripts/BaseCellCalling/BaseCellCalling.step1.py \
	--infile ${output_step3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
	--outfile Step4_VariantCalling/${sample} \
        --ref $ref

	python ${y}/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile Step4_VariantCalling/${sample}.calling.step1.tsv \
          --outfile Step4_VariantCalling/${sample} \
          --editing $editing \
          --pon $pon
	"""

}


workflow {

	//outputDirectory=Channel.fromPath(params.res_dir)
    	//create_result_directory(outputDirectory)

    	//scomat_path=Channel.fromPath(params.res_dir)
	
	//step1_dir=Channel.fromPath(params.output_step_dir1)
	bam=Channel.fromPath(params.bam)
	bai=Channel.fromPath(params.bai)
	bamFile= Channel.fromFilePairs(params.bam_folder+"*.{bam,bam.bai}")
					.view()
	//tuple(
        //file("/data/mesomics/work/mesomics2/lipikal/nextflow-practice/SComatic/example_data/Example.scrnaseq.bam"),
        //file("/data/mesomics/work/mesomics2/lipikal/nextflow-practice/SComatic/example_data/Example.scrnaseq.bam.bai"))	
    	
	metaFile=Channel.fromPath(params.meta)
    	nTrim=Channel.of(params.nTrim)
    	maxNM=Channel.of(params.maxNM)
    	maxNH=Channel.of(params.maxNH)

    step1(params.scomat_path, bamFile, metaFile, nTrim, maxNM, maxNH)
	
	//step2_dir = Channel.fromPath(params.output_step_dir2)
	Ref = Channel.fromPath(params.ref)
	minbq = Channel.of(params.minbq)
	nprocs= Channel.of(params.nprocs)
	chrom= Channel.of(params.chrom)
	
	step2(params.scomat_path, step1.out, Ref, chrom, minbq, nprocs)

	sample = Channel.of(params.sample)
	
	step3(params.scomat_path, step2.out, sample)

	//step3_dir=Channel.fromPath(params.output_step_dir3)
	//step4_dir=Channel.fromPath(params.output_step_dir4)
	edit=Channel.fromPath(params.editing)
	PON=Channel.fromPath(params.pon)
	
	step4(params.scomat_path, step3.out, sample, Ref, edit, PON)
	
}
