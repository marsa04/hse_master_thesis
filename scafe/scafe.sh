# run scafe container with mounted local folders
docker run -it -v /home:/home marsa

# look on cr-outs 
ls -lh ../home/Shared/m_sabirov/cr-outs/

# get full path to folder
realpath ../home/Shared/m_sabirov/cr-outs/
path=/home/Shared/m_sabirov/cr-outs

# create list with sample names
samples=(K2_1 K2_2 K2_3 K2_4 K2_5 K2_6 Rat1 Rat2 Rat3 Rat4 Rat5 Rat8)
echo ${samples[@]}

# check paths to bams and barcodes
for sample in ${samples[@]}; do 
	echo $path/${sample}_cellout/outs/possorted_genome_bam.bam 
	echo $path/${sample}_cellout/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; done

# run scafe.workflow.sc.solo 
for sample in ${samples[@]}; do 
	scafe.workflow.sc.solo \
	--max_thread=10 \
	--run_bam_path=$path/${sample}_cellout/outs/possorted_genome_bam.bam \
	--run_cellbarcode_path=$path/${sample}_cellout/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	--genome=mRatBN7.2 \
	--run_tag=$sample \
	--run_outDir=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.sc.solo.outs/ ; done




######################### AGREGATION ###################################

## create lib_list_path.txt 
samples=(K2_1 K2_2 K2_3 K2_4 K2_5 K2_6 Rat1 Rat2 Rat3 Rat4 Rat5 Rat8)
path=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.sc.solo.outs/bam_to_ctss

for sample in ${samples[@]}; do 
	echo $sample
	echo $path/$sample/bed/${sample}.collapse.ctss.bed.gz
	echo $path/$sample/bed/${sample}.unencoded_G.collapse.ctss.bed.gz ; done

# run scafe.workflow.cm.aggregate

scafe.workflow.cm.aggregate \
--max_thread=10 \
--lib_list_path=/home/Shared/m_sabirov/scafe.outs/lib_list_path.txt \
--genome=mRatBN7.2 \
--run_tag=agrsv.tamed \
--run_outDir=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.cm.agregate.outs/

# The output of scafe.workflow.cm.aggregate is library agnostic, it only output the CRE definition 
# and you have to run scafe.tool.sc.count using the aggregated CREs to generate count matrix per library.



################# scafe.tool.sc.solo ################################
# re run tool count with aggreagte data https://github.com/chung-lab/SCAFE/issues/18

samples=(K2_1 K2_2 K2_3 K2_4 K2_5 K2_6 Rat1 Rat2 Rat3 Rat4 Rat5 Rat8)
path=/home/Shared/m_sabirov/cr-outs

for sample in ${samples[@]}; do 
	scafe.tool.sc.count \
	--countRegion_bed_path=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.cm.agregate.outs/annotate/agrsv.tamed/bed/agrsv.tamed.CRE.coord.bed.gz \
	--cellBarcode_list_path=$path/${sample}_cellout/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
	--ctss_bed_path=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.sc.solo.outs/bam_to_ctss/$sample/bed/${sample}.CB.ctss.bed.gz \
	--genome=mRatBN7.2 \
	--ctss_scope_bed_path=/home/Shared/m_sabirov/scafe.outs/scafe.workflow.cm.agregate.outs/filter/agrsv.tamed/bed/agrsv.tamed.tssCluster.default.filtered.bed.gz \
	--outputPrefix=$sample \
	--outDir=/home/Shared/m_sabirov/scafe.outs/scafe.tool.sc.solo.aggregate.outs ; done







