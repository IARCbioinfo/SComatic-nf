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
    log.info "nextflow iarcbioinfo/SComatic-nf --bam_folder bamFile.bam --meta metadata.tsv --ref ref.fa --scomat_path path/to/Scomatfolder --annovar_path path/to/annovar --hdb path/to/humandb --hg_build GRCh38"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--scomat_path   <DIR>   SComatic install folder"
    log.info "--bam_folder    <DIR>   Folder containing BAM files"
    log.info "--meta    <FILE>   Metadata file mapping cell barcodes to cell type."
    log.info "--ref    <FILE>   genome reference files (with index)"
    log.info "--annovar_path   <DIR>   path to annovar"
    log.info "--hdb    <DIR>   path to human database for annotation (refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome- required)"
    log.info "--hg_build   <STRING>   genome build"

    log.info ""
    log.info "Optional arguments:"
    log.info "--output_folder      <DIR>   Output folder"
    log.info "--cpu      <VAL>   Number of CPUs"
    log.info "--mem      <VAL>   Memory allocated"
    log.info "--nTrim      <VAL>   Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read"
    log.info "--maxNM      <VAL>   Maximum number of mismatches permitted to consider reads for analysis"
    log.info "--maxNH     <VAL>   Maximum number of alignment hits permitted to consider reads for analysis"
    log.info "--chrom      <STRING/VAL>   Chromosome to be analysed"
    log.info "--minbq     <VAL>   Minimum base quality permited for the base counts"
    log.info "--nprocs      <VAL>   Number of processes"
    log.info "--pon      <VAL>   Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts"
    log.info "--min_signatures      <VAL>   Minimum number of Mutational signatures"
    log.info "--max_signatures      <VAL>   Maximum number of Mutational signatures"

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
params.min_signatures=1
params.max_signatures=10

//step1 paramters
params.scomat_path="SComatic"
params.bam_folder="Example.scrnaseq.bam"
//params.bai="Example.scrnaseq.bam.bai"
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

//annovar parameters


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
	for bam in \$(ls -d $output_step1/*bam);do

 		cell_type=\$(basename \$bam | awk -F '.' '{print \$(NF-1)}')

  		# Temp folder
  		temp=\$(echo Step2_BaseCellCounts/temp_\${cell_type})
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

process annovar {
	maxForks 1
	
    input:
    path output_step4
    path annovar_path
    path hdb
    val sample
    

	output:
	path("Step5_annovar")

	publishDir "${params.output_folder}", mode: 'copy'

	script:
	"""
	mkdir -p Step5_annovar

	grep -v '#' ${output_step4}/${sample}.calling.step2.tsv | tr '\t' '-' | awk -F'-' -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$5,\$0}' > Step5_annovar/${sample}.variants.avinput
        perl ${annovar_path}/table_annovar.pl Step5_annovar/${sample}.variants.avinput \
                ${hdb} -buildver hg38 \
                -out Step5_annovar/${sample}.variants.annovar \
                -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
                -nastring . -csvout -polish --otherinfo
    	"""
}

process preprocessing {
	maxForks 1
	
	input:
	path output_step4
	val sample

	output:
	path("sigprofiler-input")

	publishDir "${params.output_folder}", mode: 'copy'
	
	script:
	"""
	mkdir -p sigprofiler-input
	
	
	tail -n +29 ${output_step4}/${sample}.calling.step2.tsv | awk -F'\t' 'BEGIN {OFS="\t"} {print \$1, \$2, \$3, \$4, \$5, \$7, \$24}' > tmp.tsv
	grep '^chr' tmp.tsv > tmp1.tsv

	head -n 1 tmp.tsv > column-names.tsv
	cat column-names.tsv tmp1.tsv >> tmp4.tsv

	
	sed 's/^chr//' tmp4.tsv > tmp5.tsv
  	awk 'NR > 1 {{gsub(/chr/, "", \$1); split(\$5, altArr, ","); split(\$6, cellTypesArr, ","); split(\$7, cellTypeFilterArr, ","); for (i in altArr) { print \$1, \$2, \$3, \$4, altArr[i], cellTypesArr[i], cellTypeFilterArr[i] } }}' tmp5.tsv > tmp6.tsv
	
	
	head -n 1 tmp.tsv > column-names.tsv
	cat column-names.tsv tmp6.tsv >> tmp7.tsv
	awk '{if (\$7=="PASS") print \$0}' tmp7.tsv > tmp8.tsv	
	
	
	head -n 1 tmp.tsv > column-names.tsv
	\\cat column-names.tsv tmp8.tsv > tmp9.tsv	
	awk '{{\$3=""; sub("\t\t","\t")} print \$0}' tmp8.tsv > tmp10.tsv
	\\cat tmp10.tsv | sed "1 s/Cell_type_Filter/FILTER/" >  tmp11.tsv	
	awk '{print \$1, \$2, \$3, \$4, \$6, \$5}' tmp10.tsv >  tmp12.tsv
	unique_cell_types=\$(awk '{print \$6}' tmp12.tsv | sort -u) 
	
	
	for cell_type in \$unique_cell_types; do

		awk -v cell_type=\$cell_type '\$6==cell_type' tmp12.tsv > tmp13.tsv
        	awk 'BEGIN{ OFS="	"} {\$NF=""; print \$1, \$2, \$3, \$4, \$5 }' tmp13.tsv > tmp15.tsv
		awk 'BEGIN{ OFS="	"}{print \$1, \$2, \$5, \$3, \$4 }' tmp15.tsv > tmp16.tsv

		echo "#CHROM POS FILTER REF ALT" > header.txt
		cat header.txt tmp16.tsv > tmp17.tsv
		# Add header row
   		# Write the data frame to a file with the cell type name
   		file_name="${sample}_\$(echo "\$cell_type" | sed 's/ /_/g').vcf"
   		cat  tmp17.tsv > sigprofiler-input/\$file_name
	done
	"""
}

process Step5_sigprofiler{
	maxForks 1
	input:
	val hg_build
	path output_preprocessing
	val min_signatures
	val max_signatures

	output:
	path("sigprofiler-results")

        publishDir "${params.output_folder}", mode: 'copy'
		
	script:
	"""
	#!/home/lipikal/miniconda3/envs/sigpro/bin/python
	
	from SigProfilerExtractor import sigpro as sig
	def main_function():
    		sig.sigProfilerExtractor("vcf", "sigprofiler-results", "$output_preprocessing",reference_genome="$hg_build", minimum_signatures=$min_signatures, maximum_signatures=$max_signatures)
		
	if __name__=="__main__":
		main_function()	
	"""
}

workflow {

	//outputDirectory=Channel.fromPath(params.res_dir)
    	//create_result_directory(outputDirectory)

    	//scomat_path=Channel.fromPath(params.res_dir)
	
	//step1_dir=Channel.fromPath(params.output_step_dir1)
	//bam=Channel.fromPath(params.bam)
	//bai=Channel.fromPath(params.bai)
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

	annovar(step4.out, params.annovar_path, params.hdb, sample)
	
	preprocessing(step4.out, sample)

	Step5_sigprofiler(params.hg_build, preprocessing.out, params.min_signatures, params.max_signatures)
	
	
}
