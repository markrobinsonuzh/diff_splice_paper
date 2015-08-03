#!/bin/bash

## Define paths to software and reference files
BASEDIR=/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte
SRATOOLKIT=$BASEDIR/software/sratoolkit.2.5.0-1-ubuntu64
REFERENCEDIR=$BASEDIR/hsapiens/reference_files
RSEM=$BASEDIR/software/rsem-1.2.21
NULLSIMULATION_NODE=$BASEDIR/hsapiens/no_diffexpression/null_simulation
NONNULLSIMULATION_NODE=$BASEDIR/hsapiens/no_diffexpression/non_null_simulation
NULLSIMULATION_DE=$BASEDIR/hsapiens/with_diffexpression/null_simulation
NONNULLSIMULATION_DE=$BASEDIR/hsapiens/with_diffexpression/non_null_simulation
ASTALAVISTA=$BASEDIR/software/astalavista-3.2/bin
FIGDIR=$BASEDIR/hsapiens/figures
RCODEGEN=$BASEDIR/software/Rcode
ROUT=$BASEDIR/hsapiens/Rout
DEXSEQ=/home/Shared/Rlib/release-3.1-lib/DEXSeq/python_scripts
TOPHAT=$BASEDIR/software/tophat-2.0.14.Linux_x86_64
CUFFLINKS=$BASEDIR/software/cufflinks-2.2.1.Linux_x86_64
MISOENV=$BASEDIR/software/miso/misoenv2
MISOINDEXDIR=$REFERENCEDIR/miso_index/
MISOINDEXDIR_MISSING20=$REFERENCEDIR/INCOMPLETE_MISSING20/miso_index/
KALLISTO=$BASEDIR/software/kallisto_linux-v0.42.1
bedGraphToBigWig=$BASEDIR/software/bedGraphToBigWig
rMATS=$BASEDIR/software/rMATS.3.0.9

## ------------------------- INPUT PREPARATION ----------------------------- ##

## The basis for the simulation is generated from a sample downloaded from the SRA. 
## This sample was sequenced with an Illumina HiSeq, with a paired-end protocol,
## and with a read length of 101 bp.
## Extract the fastq files from the sra archive
$SRATOOLKIT/fastq-dump --split-files \
-Q 33 -O $REFERENCEDIR $REFERENCEDIR/SRR493366.sra

## Extract only lines corresponding to primary assembly chromosomes from gtf file
grep -v "^H" $REFERENCEDIR/Homo_sapiens.GRCh37.71.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf

## Extract only lines corresponding to protein coding genes from gtf file
grep "protein_coding" $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf

## Prepare the reference files for RSEM
$RSEM/rsem-prepare-reference \
--gtf $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf --bowtie2 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71

## Estimate the model file with RSEM
$RSEM/rsem-calculate-expression \
--paired-end --bowtie2 --seed 123 \
$REFERENCEDIR/SRR493366_1.fastq $REFERENCEDIR/SRR493366_2.fastq \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
$REFERENCEDIR/rsem_model/SRR493366

## Plot some characteristics of the estimated model
$RSEM/rsem-plot-model $REFERENCEDIR/rsem_model/SRR493366 $FIGDIR/RSEM_model.pdf

## Modify the quality score distribution in the RSEM model file so that the probability 
## of quality score 2 is 0. Otherwise, the quality scores of the simulated data may be very low.
## Also make sure that the transition probabilities into this state are 0.
## ->> $REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ.model

## Estimate and plot the isoform percentage distributions from the RSEM results
R CMD BATCH --no-restore --no-save "--args referencefile='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' outdir='$FIGDIR'" $RCODEGEN/isopct_distribution.R $ROUT/isopct_distribution_human.Rout

## Run ASTALAVISTA to classify splicing events
$ASTALAVISTA/astalavista -t asta \
-i $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf -e [ASE,ASI,DSP,VST]
gunzip $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf.gz

## Duplicate the protein_coding gtf file into one with the same name as the genome
## fasta file, for index building (not really necessary)
scp $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf

## Build TopHat index
bowtie2-build -f $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly
ln -s $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa

## Build TopHat transcriptome index
tophat -G $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf \
--transcriptome-index=$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly

## Build kallisto index from the fasta file created in the TopHat transcriptome index
$KALLISTO/kallisto index \
--index=$REFERENCEDIR/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa

## Create conversion table from transcript index to transcript name (to interpret kallisto output)
grep "^>" $REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa | sed -e 's/>//' | cut -d" " -f1,2 > $REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt

## Prepare gtf file for cuffdiff
$CUFFLINKS/cuffcompare -s $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
-CG -r $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf \
-o $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_cuffdiff \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf

## Prepare gtf file for cuffdiff, where Ensembl gene IDs are used instead of symbols
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf' output_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl.gtf'" $RCODEGEN/generate_gtf_cuffdiff.R $ROUT/generate_gtf_cuffdiff_human.Rout
$CUFFLINKS/cuffcompare -s $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
-CG -r $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl.gtf \
-o $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_cuffdiff \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl.gtf

## Prepare gtf file for cuffdiff, where all exons are renamed to CDS, and all CDS removed,
## and Ensembl gene IDs are used instead of symbols
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf' output_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_exontocds.gtf'" $RCODEGEN/generate_gtf_exontocds.R $ROUT/generate_gtf_exontocds_human.Rout
$CUFFLINKS/cuffcompare -s $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
-CG -r $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_exontocds.gtf \
-o $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_exontocds_cuffdiff \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_exontocds.gtf

## Prepare flattened annotations (for DEXSeq)
python $DEXSEQ/dexseq_prepare_annotation.py \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff

python $DEXSEQ/dexseq_prepare_annotation.py --aggregate='no' \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq)
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign.gtf' ignore_strand=TRUE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf.Rout
sed -i -e 's/[*]/./' $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq) - strand specific overlaps
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.strandspec.gtf' ignore_strand=FALSE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf_strandspec.Rout
sed -i -e 's/[*]/./' $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.strandspec.gff

## Fix original gtf file to count on (real) exon level with featureCounts
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene.gtf'" $RCODEGEN/generate_renamed_gtf_for_featurecounts.R $ROUT/generate_renamed_gtf_for_featurecounts.Rout

## Prepare gtf file with 'chr' as chromosome prefix (for Casper)
sed -e 's/^/chr/' $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr.gtf

## -------------------------- DATA SIMULATION ------------------------------ ##

## Estimate the simulation parameters for the individual samples
R CMD BATCH --no-restore --no-save "--args path_to_generate_rsem_files='$RCODEGEN/generate_rsem_files_function.R' seed=123 isoform_results_file='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' nbr_per_group=3 meandisp.file='$REFERENCEDIR/Pickrell.Cheung.Mu.Phi.Estimates.rds' outdirbase='$BASEDIR/hsapiens/no_diffexpression' librarysize=40000000 keepchr=NULL nbr_diff_spliced=1000 nbr_diff_expr=0 fold_changes=NULL" $RCODEGEN/generate_rsem_files_human_run.R $ROUT/generate_rsem_files_human_run_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_generate_rsem_files='$RCODEGEN/generate_rsem_files_function.R' seed=123 isoform_results_file='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' nbr_per_group=3 meandisp.file='$REFERENCEDIR/Pickrell.Cheung.Mu.Phi.Estimates.rds' outdirbase='$BASEDIR/hsapiens/with_diffexpression' librarysize=40000000 keepchr=NULL nbr_diff_spliced=1000 nbr_diff_expr=1000 fold_changes='expon'" $RCODEGEN/generate_rsem_files_human_run.R $ROUT/generate_rsem_files_human_run_de.Rout

R CMD BATCH --no-restore --no-save "--args path_to_sim_details='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' output_pdf='$FIGDIR/simulation_details_nodiffexp_nonnull.pdf'" $RCODEGEN/plot_simulation_details_2.R $ROUT/plot_simulation_details_2_node_nonnull.Rout
R CMD BATCH --no-restore --no-save "--args path_to_sim_details='$NONNULLSIMULATION_DE/3_truth/simulation_details.txt' output_pdf='$FIGDIR/simulation_details_withdiffexp_nonnull.pdf'" $RCODEGEN/plot_simulation_details.R $ROUT/plot_simulation_details_withde_nonnull.Rout

## Generate truth files
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NULLSIMULATION_DE/3_truth/simulation_details.txt' out.file='$NULLSIMULATION_DE/3_truth/truth_human_null.txt' astalavista.file=NULL gtf.file=NULL flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_null_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NULLSIMULATION_NODE/3_truth/simulation_details.txt' out.file='$NULLSIMULATION_NODE/3_truth/truth_human_null.txt' astalavista.file=NULL gtf.file=NULL flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_null_node.Rout

R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NONNULLSIMULATION_DE/3_truth/simulation_details.txt' out.file='$NONNULLSIMULATION_DE/3_truth/truth_human_non_null.txt' astalavista.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' out.file='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' astalavista.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_nonnull_node.Rout

## Simulate reads for 6 samples (40 million pairs/sample),
## for non-null situation, with and without differential expression
for n in 1 2 3 4 5 6
do
	$RSEM/rsem-simulate-reads \
	$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
	$REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ.model \
	$NONNULLSIMULATION_NODE/1_reads/rsem_files/sample${n}.txt \
	0.05 40000000 $NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n} \
	--seed 123

	$RSEM/rsem-simulate-reads \
	$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
	$REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ.model \
	$NONNULLSIMULATION_DE/1_reads/rsem_files/sample${n}.txt \
	0.05 40000000 $NONNULLSIMULATION_DE/1_reads/reads/sample${n}/sample_${n} \
	--seed 123
done	

## -------------------------- READ ALIGNMENT ------------------------------- ##

## Align with TopHat
for n in 1 2 3 4 5 6
do
	$TOPHAT/tophat2 -p 6 --no-coverage-search \
	--transcriptome-index=$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	-o $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n} \
	$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_2.fq
	
	samtools index $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam
done

for n in 1 2 3 4 5 6
do 
	$TOPHAT/tophat2 -p 6 --no-coverage-search \
	--transcriptome-index=$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	-o $NONNULLSIMULATION_DE/1_reads/tophat/sample${n} \
	$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	$NONNULLSIMULATION_DE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_DE/1_reads/reads/sample${n}/sample_${n}_2.fq
	
	samtools index $NONNULLSIMULATION_DE/1_reads/tophat/sample${n}/accepted_hits.bam
done

## Subsample the aligned read files to approximately 60% of the original number, sort and index
for n in 1 2 3 4 5 6
do
	samtools view -b \
	-o $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits.bam \
	-s 0.6 $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam
	
	samtools sort -T $NONNULLSIMULATION_NODE/1_reads/tophat/tmp_${n} \
	-o $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits_sorted.bam \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits.bam
	
	samtools index $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits_sorted.bam
done

## Get chromosome lengths from SAM header
samtools view -H $NONNULLSIMULATION_NODE/1_reads/tophat/sample1/accepted_hits.bam | \
grep "@SQ" | cut -f2,3 | sed -e 's/SN://' | sed -e 's/LN://' > \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_chromosome_lengths.txt

## Convert bam files to bigWig
for n in 1 2 3 4 5 6
do
	bedtools genomecov -split \
	-ibam $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
	-bg > $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bedGraph
	
	bedGraphToBigWig $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bedGraph \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_chromosome_lengths.txt \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bw
	
	rm -rf $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bedGraph
done

## ------------------------------- DEXSEQ ---------------------------------- ##
## Count reads for DEXSeq
for n in 1 2 3 4 5 6
do 
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff \
	$NONNULLSIMULATION_DE/1_reads/tophat/sample${n}/accepted_hits.bam \
	$NONNULLSIMULATION_DE/2_counts/dexseq/dexseq${n}.txt
	
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
	$NONNULLSIMULATION_NODE/2_counts/dexseq/dexseq${n}.txt
done

## Run DEXSeq/HTSeq analysis
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/dexseq' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff' out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_htseq' method_name='dexseq_htseq'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq' method_name='dexseq_htseq'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node.Rout

## ------------------------------ cuffdiff --------------------------------- ##
$CUFFLINKS/cuffdiff -o $NONNULLSIMULATION_NODE/2_counts/cuffdiff -p 6 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_cuffdiff.combined.gtf \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample1/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample2/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample3/accepted_hits.bam \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample4/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample5/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample6/accepted_hits.bam

R CMD BATCH --no-restore --no-save "--args path_to_cuffdiff_result='$NONNULLSIMULATION_NODE/2_counts/cuffdiff/cds.diff' path_to_gtf_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf' output_file='$NONNULLSIMULATION_NODE/4_results/cuffdiff.txt' method_name='cuffdiff'" $RCODEGEN/clean_cuffdiff.R $ROUT/clean_cuffdiff_human.Rout

## ---------------------- cuffdiff with Ensembl IDs ------------------------ ##
## TODO: Run (done, just check)
$CUFFLINKS/cuffdiff -o $NONNULLSIMULATION_NODE/2_counts/cuffdiff_ensembl -p 6 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_cuffdiff.combined.gtf \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample1/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample2/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample3/accepted_hits.bam \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample4/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample5/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample6/accepted_hits.bam

R CMD BATCH --no-restore --no-save "--args path_to_cuffdiff_result='$NONNULLSIMULATION_NODE/2_counts/cuffdiff_ensembl/cds.diff' output_file='$NONNULLSIMULATION_NODE/4_results/cuffdiff_ensembl.txt' method_name='cuffdiff_ensembl'" $RCODEGEN/clean_cuffdiff_ensembl.R $ROUT/clean_cuffdiff_ensembl_human.Rout

## ---------------------- cuffdiff with exons as cds ----------------------- ##
## TODO: Rerun (done, just check)
$CUFFLINKS/cuffdiff -o $NONNULLSIMULATION_NODE/2_counts/cuffdiff_exontocds -p 6 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly_ensembl_exontocds_cuffdiff.combined.gtf \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample1/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample2/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample3/accepted_hits.bam \
$NONNULLSIMULATION_NODE/1_reads/tophat/sample4/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample5/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample6/accepted_hits.bam

R CMD BATCH --no-restore --no-save "--args path_to_cuffdiff_result='$NONNULLSIMULATION_NODE/2_counts/cuffdiff_exontocds/cds.diff' output_file='$NONNULLSIMULATION_NODE/4_results/cuffdiff_exontocds.txt' method_name='cuffdiff_exontocds'" $RCODEGEN/clean_cuffdiff_ensembl.R $ROUT/clean_cuffdiff_ensembl_human_exontocds.Rout

## ------------------------------- rMATS ----------------------------------- ##
## Must use samtools 0.1.19
export PATH=$BASEDIR/software/samtools-0.1.19:$PATH
python $rMATS/RNASeq-MATS.py \
-b1 $NONNULLSIMULATION_NODE/1_reads/tophat/sample1/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample2/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample3/accepted_hits.bam \
-b2 $NONNULLSIMULATION_NODE/1_reads/tophat/sample4/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample5/accepted_hits.bam,$NONNULLSIMULATION_NODE/1_reads/tophat/sample6/accepted_hits.bam \
-gtf $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf \
-o $NONNULLSIMULATION_NODE/2_counts/rMATS -t paired -len 101 \
-r1 175,175,175 -r2 175,175,175 -sd1 56,56,56 -sd2 56,56,56 -c 0 -analysis U
export PATH=$(echo $PATH | sed -e 's;:\?/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/software/samtools-0.1.19:;;' -e 's;$/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/software/samtools-0.1.19:\?;;')

R CMD BATCH --no-restore --no-save "--args mats_dir='$NONNULLSIMULATION_NODE/2_counts/rMATS/MATS_output' output_file='$NONNULLSIMULATION_NODE/4_results/rMATS.txt' method_name='rMATS'" $RCODEGEN/rmats_mergeres_run.R $ROUT/rmats_mergeres_run_human.Rout

## -------------- DEXSEQ with non-merged flattened file -------------------- ##
## Count reads for DEXSeq
for n in 1 2 3 4 5 6
do 
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff \
	$NONNULLSIMULATION_DE/1_reads/tophat/sample${n}/accepted_hits.bam \
	$NONNULLSIMULATION_DE/2_counts/dexseq_nomerge/dexseq${n}.txt
	
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
	$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge/dexseq${n}.txt
	
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits.bam \
	$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge_sub0.6/dexseq${n}.txt
done

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_htseq_nomerge' method_name='dexseq_htseq_nomerge'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_de_nomerge.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq_nomerge' method_name='dexseq_htseq_nomerge'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_nomerge.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq_nomerge_sub0.6' method_name='dexseq_htseq_nomerge_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_nomerge_sub0.6.Rout
Rold CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_BM_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_BM_htseq_nomerge' method_name='dexseq_BM_htseq_nomerge'" $RCODEGEN/dexseq_BM_run.R $ROUT/dexseq_BM_run_human_nonnull_node_nomerge.Rout
Rold CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_1.8_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_1.8_htseq_nomerge' method_name='dexseq_1.8_htseq_nomerge'" $RCODEGEN/dexseq_1.8_run.R $ROUT/dexseq_1.8_run_human_nonnull_node_nomerge.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq_nomerge.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_htseq_nomerge/nonnnull_node_dexseq_htseq_nomerge' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='HTSeq/DEXSeq'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_htseq_nomerge_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_BM_htseq_nomerge.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_bm_htseq_nomerge/nonnnull_node_dexseq_bm_htseq_nomerge' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='HTSeq/DEXSeq_BM'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_bm_results_htseq_nomerge_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_1.8_htseq_nomerge.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_1.8_htseq_nomerge/nonnnull_node_dexseq_1.8_htseq_nomerge' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='HTSeq/DEXSeq_1.8'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_1.8_results_htseq_nomerge_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_htseq_nomerge' method_name='voom_diffsplice_htseq_nomerge'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_de_nomerge.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_htseq_nomerge' method_name='voom_diffsplice_htseq_nomerge'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_node_nomerge.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_htseq_nomerge_sub0.6' method_name='voom_diffsplice_htseq_nomerge_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_node_nomerge_sub0.6.Rout

## Summarise DEXSeq exon p-values in other ways
R CMD BATCH --no-restore --no-save "--args path_to_flattened_gff='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' path_to_rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq_nomerge.Rdata' out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_htseq_nomerge_diffclass' method_name='dexseq_htseq_nomerge_diffclass'" $RCODEGEN/summarise_pvalues_dexseq.R $ROUT/summarise_pvalues_dexseq_nonnull_node.Rout

## -------------------------------- MISO ----------------------------------- ##

perl $BASEDIR/software/gtf_to_gff.pl \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf > \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gff3

## Activate the local environment for MISO
source $MISOENV/bin/activate

index_gff --index $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gff3 $MISOINDEXDIR

for n in 1 2 3 4 5 6
do
	miso --run $MISOINDEXDIR $NONNULLSIMULATION_DE/1_reads/tophat/sample${n}/accepted_hits.bam \
	-p 6 --output-dir $NONNULLSIMULATION_DE/2_counts/miso/output/sample${n} \
	--read-len 101 --paired-end 175 56
	
	miso --run $MISOINDEXDIR $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
	-p 6 --output-dir $NONNULLSIMULATION_NODE/2_counts/miso/output/sample${n} \
	--read-len 101 --paired-end 175 56

	miso --run $MISOINDEXDIR $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_sub0.6/accepted_hits_sorted.bam \
	-p 6 --output-dir $NONNULLSIMULATION_NODE/2_counts/miso_sub0.6/output/sample${n} \
	--read-len 101 --paired-end 175 56
done

## Deactivate the local environment
deactivate

## Summarize miso output
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_DE/2_counts/miso/output/' counts_output_dir='$NONNULLSIMULATION_DE/2_counts/miso'" $RCODEGEN/miso_run.R $ROUT/miso_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso/output/' counts_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso'" $RCODEGEN/miso_run.R $ROUT/miso_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso_sub0.6/output/' counts_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso_sub0.6'" $RCODEGEN/miso_run.R $ROUT/miso_run_human_nonnull_node_sub0.6.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/miso' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_miso' method_name='dexseq_miso'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_miso' method_name='dexseq_miso'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_miso_sub0.6' method_name='dexseq_miso_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_node_sub0.6.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_miso.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_miso/nonnnull_node_dexseq_miso' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='MISO'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_miso_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/miso' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_miso' method_name='voom_diffsplice_miso'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_miso_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_miso' method_name='voom_diffsplice_miso'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_miso_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_miso_sub0.6' method_name='voom_diffsplice_miso_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_miso_node_sub0.6.Rout

## Subset MISO counts to only those assignable to at least one isoform and apply DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_input_files='$NONNULLSIMULATION_NODE/2_counts/miso' path_to_output_files='$NONNULLSIMULATION_NODE/2_counts/miso_assignable'" $RCODEGEN/subset_miso_counts_assignable.R $ROUT/subset_miso_counts_assignable_human.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso_assignable' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_miso_assignable' method_name='dexseq_miso_assignable'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_assignable_node.Rout

## Quantify transcripts instead of transcript combinations and apply DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso/output/' counts_output_dir='$NONNULLSIMULATION_NODE/2_counts/miso_transcripts'" $RCODEGEN/miso_transcripts_run.R $ROUT/miso_transcripts_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/miso_transcripts' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_miso_transcripts' method_name='dexseq_miso_transcripts'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_transcripts_node.Rout

## ---------------------- TopHat junction counts --------------------------- ##
R CMD BATCH --no-restore --no-save "--args gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' tophat.output.folder='$NONNULLSIMULATION_DE/1_reads/tophat' pattern='sample[0-9]$' output.dir='$NONNULLSIMULATION_DE/2_counts/tophat_junc'" $RCODEGEN/tophatjunction_run.R $ROUT/tophatjunction_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' tophat.output.folder='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]$' output.dir='$NONNULLSIMULATION_NODE/2_counts/tophat_junc'" $RCODEGEN/tophatjunction_run.R $ROUT/tophatjunction_run_human_nonnull_node.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/tophat_junc' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_tophat_junc' method_name='dexseq_tophat_junc'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_tophat_junc_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/tophat_junc' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_tophat_junc' method_name='dexseq_tophat_junc'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_tophat_junc_node.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_tophat_junc.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_tophat_junc/nonnnull_node_dexseq_tophat_junc' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='TopHat_junctions'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_tophat_junc_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/tophat_junc' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_tophat_junc' method_name='voom_diffsplice_tophat_junc'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_tophat_junc_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/tophat_junc' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_tophat_junc' method_name='voom_diffsplice_tophat_junc'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_tophat_junc_node.Rout

## ------------------------------- casper ---------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_casper_fcn='$RCODEGEN/casper_function.R' path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr.gtf' path_to_tophat='$NONNULLSIMULATION_DE/1_reads/tophat' pattern='sample[0-9]$' genome='hg19' read_length=101 output.dir='$NONNULLSIMULATION_DE/2_counts/casper'" $RCODEGEN/casper_run.R $ROUT/casper_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_casper_fcn='$RCODEGEN/casper_function.R' path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr.gtf' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]$' genome='hg19' read_length=101 output.dir='$NONNULLSIMULATION_NODE/2_counts/casper'" $RCODEGEN/casper_run.R $ROUT/casper_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_casper_fcn='$RCODEGEN/casper_function.R' path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr.gtf' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_sub0.6$' genome='hg19' read_length=101 output.dir='$NONNULLSIMULATION_NODE/2_counts/casper_sub0.6'" $RCODEGEN/casper_run.R $ROUT/casper_run_human_nonnull_node_sub0.6.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/casper' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_casper' method_name='dexseq_casper'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_casper_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/casper' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_casper' method_name='dexseq_casper'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_casper_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/casper_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_casper_sub0.6' method_name='dexseq_casper_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_casper_node_sub0.6.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_casper.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_casper/nonnnull_node_dexseq_casper' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='casper'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_casper_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/casper' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_casper' method_name='voom_diffsplice_casper'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_casper_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/casper' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_casper' method_name='voom_diffsplice_casper'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_casper_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/casper_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_casper_sub0.6' method_name='voom_diffsplice_casper_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_casper_node_sub0.6.Rout

## --------------------------- SplicingGraph ------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_splicinggraph_fcn='$RCODEGEN/splicinggraph_function.R' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output.dir='$NONNULLSIMULATION_DE/2_counts/splicinggraph' path_to_tophat='$NONNULLSIMULATION_DE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/splicinggraph_run.R $ROUT/splicinggraph_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_splicinggraph_fcn='$RCODEGEN/splicinggraph_function.R' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/splicinggraph' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/splicinggraph_run.R $ROUT/splicinggraph_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_splicinggraph_fcn='$RCODEGEN/splicinggraph_function.R' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/splicinggraph_sub0.6' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_sub0.6$'" $RCODEGEN/splicinggraph_run.R $ROUT/splicinggraph_run_human_nonnull_node_sub0.6.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/splicinggraph' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_splicinggraph' method_name='dexseq_splicinggraph'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_splicinggraph_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/splicinggraph' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_splicinggraph' method_name='dexseq_splicinggraph'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_splicinggraph_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/splicinggraph_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_splicinggraph_sub0.6' method_name='dexseq_splicinggraph_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_splicinggraph_node_sub0.6.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_splicinggraph.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_splicinggraph/nonnnull_node_dexseq_splicinggraph' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='SplicingGraph'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_splicinggraph_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/splicinggraph' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_splicinggraph' method_name='voom_diffsplice_splicinggraph'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_splicinggraph_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/splicinggraph' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_splicinggraph' method_name='voom_diffsplice_splicinggraph'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_splicinggraph_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/splicinggraph_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_splicinggraph_sub0.6' method_name='voom_diffsplice_splicinggraph_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_splicinggraph_node_sub0.6.Rout

## --------------------------- featureCounts ------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign.gtf' output.dir='$NONNULLSIMULATION_DE/2_counts/featurecounts_flat' path_to_tophat='$NONNULLSIMULATION_DE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat_sub0.6' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_sub0.6$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_run_human_nonnull_node_sub0.6.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/featurecounts_flat' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_featurecounts_flat' method_name='dexseq_featurecounts_flat'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_flat_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_featurecounts_flat' method_name='dexseq_featurecounts_flat'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_flat_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_featurecounts_flat_sub0.6' method_name='dexseq_featurecounts_flat_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_flat_node_sub0.6.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_featurecounts_flat.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_featurecounts_flat/nonnnull_node_dexseq_featurecounts_flat' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='featureCounts'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_featurecounts_flat_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/featurecounts_flat' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_featurecounts_flat' method_name='voom_diffsplice_featurecounts_flat'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_flat_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_featurecounts_flat' method_name='voom_diffsplice_featurecounts_flat'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_flat_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_flat_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_featurecounts_flat_sub0.6' method_name='voom_diffsplice_featurecounts_flat_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_flat_node_sub0.6.Rout

## ------------------- featureCounts, not flattened ------------------------ ##
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene.gtf' output.dir='$NONNULLSIMULATION_DE/2_counts/featurecounts_noflat' path_to_tophat='$NONNULLSIMULATION_DE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_noflat_run_human_nonnull_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_noflat_run_human_nonnull_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat_sub0.6' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_sub0.6$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_noflat_run_human_nonnull_node_sub0.6.Rout

R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/featurecounts_noflat' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_featurecounts_noflat' method_name='dexseq_featurecounts_noflat'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_noflat_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_featurecounts_noflat' method_name='dexseq_featurecounts_noflat'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_noflat_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat_sub0.6' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_featurecounts_noflat_sub0.6' method_name='dexseq_featurecounts_noflat_sub0.6'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_noflat_node_sub0.6.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/featurecounts_noflat' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_featurecounts_noflat' method_name='voom_diffsplice_featurecounts_noflat'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_noflat_de.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_featurecounts_noflat' method_name='voom_diffsplice_featurecounts_noflat'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_noflat_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/featurecounts_noflat_sub0.6' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_featurecounts_noflat_sub0.6' method_name='voom_diffsplice_featurecounts_noflat_sub0.6'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_featurecounts_noflat_node_sub0.6.Rout

## ------------------------------ KALLISTO --------------------------------- ##

for n in 1 2 3 4 5 6
do
	$KALLISTO/kallisto quant \
	--index=$REFERENCEDIR/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	--output-dir=$NONNULLSIMULATION_NODE/2_counts/kallisto/sample${n} \
	--seed=123 --plaintext \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_2.fq
	
	$KALLISTO/kallisto quant \
	--index=$REFERENCEDIR/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	--output-dir=$NONNULLSIMULATION_DE/2_counts/kallisto/sample${n} \
	--seed=123 --plaintext \
	$NONNULLSIMULATION_DE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_DE/1_reads/reads/sample${n}/sample_${n}_2.fq
done

## Summarise kallisto output
R CMD BATCH --no-restore --no-save "--args path_to_kallisto_output='$NONNULLSIMULATION_NODE/2_counts/kallisto' path_to_conversion_table='$REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt' path_to_isoform_results='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' pattern='sample[0-9]$' output_dir='$NONNULLSIMULATION_NODE/2_counts/kallisto'" $RCODEGEN/kallisto_summarise_run.R $ROUT/kallisto_summarise_run_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_kallisto_output='$NONNULLSIMULATION_DE/2_counts/kallisto' path_to_conversion_table='$REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt' path_to_isoform_results='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' pattern='sample[0-9]$' output_dir='$NONNULLSIMULATION_DE/2_counts/kallisto'" $RCODEGEN/kallisto_summarise_run.R $ROUT/kallisto_summarise_run_de.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/kallisto' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/dexseq_kallisto' method_name='dexseq_kallisto'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_kallisto_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/kallisto' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_DE/4_results/dexseq_kallisto' method_name='dexseq_kallisto'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_kallisto_de.Rout

R CMD BATCH --no-restore --no-save "--args path_to_Rdata_object='$NONNULLSIMULATION_NODE/4_results/dexseq_kallisto.Rdata' output_basename='$FIGDIR/nonnull_node_dexseq_kallisto/nonnnull_node_dexseq_kallisto' path_to_truth_table='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null.txt' method_name='kallisto'" $RCODEGEN/plot_dexseq_results.R $ROUT/plot_dexseq_results_kallisto_nonnull_node.Rout

## Run voom/diffSplice
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/kallisto' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_NODE/4_results/voom_diffsplice_kallisto' method_name='voom_diffsplice_kallisto'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_kallisto_node.Rout
R CMD BATCH --no-restore --no-save "--args path_to_voom_diffsplice_fcn='$RCODEGEN/voom_diffsplice_function.R' path_to_count_files='$NONNULLSIMULATION_DE/2_counts/kallisto' conditions=rep(c('c1','c2'),each=3) out_basename='$NONNULLSIMULATION_DE/4_results/voom_diffsplice_kallisto' method_name='voom_diffsplice_kallisto'" $RCODEGEN/voom_diffsplice_run.R $ROUT/voom_diffsplice_run_human_nonnull_kallisto_de.Rout

## ---------------------- OUTPUT CHARACTERIZATION -------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_count_tables='$NONNULLSIMULATION_NODE/2_counts' count_methods=c('casper','dexseq','dexseq_nomerge','featurecounts_flat','featurecounts_noflat','miso_assignable','splicinggraph','tophat_junc','kallisto') output_file='$FIGDIR/counting_characterization_nonnull_node_human.pdf'" $RCODEGEN/characterize_counting_methods.R $ROUT/characterize_counting_methods.Rout

## ========================================================================= ##
##                         INCOMPLETE ANNOTATION                             ##
## ========================================================================= ##

## ------------------------- INPUT PREPARATION ----------------------------- ##
## Determine which transcripts to remove and generate a new gtf file. 
## Remove 20% of the differentially spliced ones,
## and 20% of the non-differentially spliced ones.
R CMD BATCH --no-restore --no-save "--args path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' path_to_truth='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' remove_fraction=0.2 output_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' seed=123 output_missingtx='$REFERENCEDIR/INCOMPLETE_MISSING20/excluded_transcripts.txt'" $RCODEGEN/generate_incomplete_gtf.R $ROUT/generate_incomplete_gtf_human.Rout

## Prepare annotations and indices as above
## Duplicate the protein_coding gtf file into one with the same name as the genome
## fasta file, for index building
scp $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf

## Build TopHat transcriptome index
tophat -G $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf \
--transcriptome-index=$REFERENCEDIR/INCOMPLETE_MISSING20/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly

## Build kallisto index from the fasta file created in the TopHat transcriptome index
$KALLISTO/kallisto index \
--index=$REFERENCEDIR/INCOMPLETE_MISSING20/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/INCOMPLETE_MISSING20/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa

## Create conversion table from transcript index to transcript name (to interpret kallisto output)
grep "^>" $REFERENCEDIR/INCOMPLETE_MISSING20/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa | sed -e 's/>//' | cut -d" " -f1,2 > $REFERENCEDIR/INCOMPLETE_MISSING20/KallistoIndex/TranscriptID_conversion.txt

## Prepare flattened annotations (for DEXSeq)
python $DEXSEQ/dexseq_prepare_annotation.py \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened_missing20.gff

python $DEXSEQ/dexseq_prepare_annotation.py --aggregate='no' \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_missing20.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq)
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' output_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign_missing20.gtf' ignore_strand=TRUE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf_missing20.Rout
sed -i -e 's/[*]/./' $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign_missing20.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq) - strand specific overlap
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' output_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.strandspec_missing20.gtf' ignore_strand=FALSE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf_strandspec_missing20.Rout
sed -i -e 's/[*]/./' $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.strandspec_missing20.gff

## Fix original gtf file to count on (real) exon level with featureCounts
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' output_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene_missing20.gtf'" $RCODEGEN/generate_renamed_gtf_for_featurecounts.R $ROUT/generate_renamed_gtf_for_featurecounts_missing20.Rout

## Prepare gtf file with 'chr' as chromosome prefix (for Casper)
sed -e 's/^/chr/' $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf \
> $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr_missing20.gtf
sed -i -e 's/chr##/##/' $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr_missing20.gtf

## Prepare gff3 file for MISO
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' output_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20_miso.gtf'" $RCODEGEN/fix_gtf_for_miso.R $ROUT/fix_gtf_for_miso_human.Rout
perl $BASEDIR/software/gtf_to_gff.pl \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20_miso.gtf > \
$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gff3

## Generate new truth file
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' out.file='$NONNULLSIMULATION_NODE/3_truth/truth_human_non_null_missing20.txt' astalavista.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file='$REFERENCEDIR/INCOMPLETE_MISSING20/excluded_transcripts.txt'" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_nonnull_node_missing20.Rout

## -------------------------- READ ALIGNMENT ------------------------------- ##
## Align with TopHat
for n in 1 2 3 4 5 6
do
	$TOPHAT/tophat2 -p 12 --no-coverage-search \
	--transcriptome-index=$REFERENCEDIR/INCOMPLETE_MISSING20/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	-o $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_incompleteannot_missing20 \
	$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_2.fq
	
	samtools index $NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_incompleteannot_missing20/accepted_hits.bam
done

## ---------------- DEXSEQ with default flattened file --------------------- ##
## Count reads for DEXSeq
for n in 1 2 3 4 5 6
do 
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened_missing20.gff \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_incompleteannot_missing20/accepted_hits.bam \
	$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/dexseq_missing20/dexseq${n}.txt
done

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/dexseq_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened_missing20.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_htseq_missing20' method_name='dexseq_htseq_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_missing20.Rout

## -------------- DEXSEQ with non-merged flattened file -------------------- ##
## Count reads for DEXSeq
for n in 1 2 3 4 5 6
do 
	python $DEXSEQ/dexseq_count.py \
	-p yes -r pos -s no -f bam \
	$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_missing20.gff \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_incompleteannot_missing20/accepted_hits.bam \
	$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/dexseq_nomerge_missing20/dexseq${n}.txt
done

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/dexseq_nomerge_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_missing20.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_htseq_nomerge_missing20' method_name='dexseq_htseq_nomerge_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_nomerge_missing20.Rout

## ---------------------- TopHat junction counts --------------------------- ##
R CMD BATCH --no-restore --no-save "--args gtf.file='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' tophat.output.folder='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_incompleteannot_missing20$' output.dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/tophat_junc_missing20'" $RCODEGEN/tophatjunction_run.R $ROUT/tophatjunction_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/tophat_junc_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_tophat_junc_missing20' method_name='dexseq_tophat_junc_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_tophat_junc_node_missing20.Rout

## ------------------------------ KALLISTO --------------------------------- ##
for n in 1 2 3 4 5 6
do
	$KALLISTO/kallisto quant \
	--index=$REFERENCEDIR/INCOMPLETE_MISSING20/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
	--output-dir=$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/kallisto_missing20/sample${n} \
	--seed=123 --plaintext \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_1.fq \
	$NONNULLSIMULATION_NODE/1_reads/reads/sample${n}/sample_${n}_2.fq
done

## Summarise kallisto output
R CMD BATCH --no-restore --no-save "--args path_to_kallisto_output='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/kallisto_missing20' path_to_conversion_table='$REFERENCEDIR/INCOMPLETE_MISSING20/KallistoIndex/TranscriptID_conversion.txt' path_to_isoform_results='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' pattern='sample[0-9]$' output_dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/kallisto_missing20'" $RCODEGEN/kallisto_summarise_run.R $ROUT/kallisto_summarise_run_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/kallisto_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_kallisto_missing20' method_name='dexseq_kallisto_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_kallisto_node_missing20.Rout

## --------------------------- SplicingGraph ------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_splicinggraph_fcn='$RCODEGEN/splicinggraph_function.R' gtf.file='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/splicinggraph_missing20' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_incompleteannot_missing20$'" $RCODEGEN/splicinggraph_run.R $ROUT/splicinggraph_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/splicinggraph_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_splicinggraph_missing20' method_name='dexseq_splicinggraph_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_splicinggraph_node_missing20.Rout

## --------------------------- featureCounts ------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flatman.ign_missing20.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/featurecounts_flat_missing20' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_incompleteannot_missing20$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/featurecounts_flat_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_featurecounts_flat_missing20' method_name='dexseq_featurecounts_flat_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_flat_node_missing20.Rout

## ------------------- featureCounts, not flattened ------------------------ ##
R CMD BATCH --no-restore --no-save "--args path_to_featurecounts_fcn='$RCODEGEN/featureCounts_function.R' gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.exongene_missing20.gtf' output.dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/featurecounts_noflat_missing20' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_incompleteannot_missing20$'" $RCODEGEN/featureCounts_run.R $ROUT/featureCounts_noflat_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/featurecounts_noflat_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_featurecounts_noflat_missing20' method_name='dexseq_featurecounts_noflat_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_featurecounts_noflat_node_missing20.Rout

## ------------------------------- casper ---------------------------------- ##
R CMD BATCH --no-restore --no-save "--args path_to_casper_fcn='$RCODEGEN/casper_function.R' path_to_gtf='$REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_chr_missing20.gtf' path_to_tophat='$NONNULLSIMULATION_NODE/1_reads/tophat' pattern='sample[0-9]_incompleteannot_missing20$' genome='hg19' read_length=101 output.dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/casper_missing20'" $RCODEGEN/casper_run.R $ROUT/casper_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/casper_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_casper_missing20' method_name='dexseq_casper_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_casper_node_missing20.Rout

## -------------------------------- MISO ----------------------------------- ##
## Activate the local environment for MISO
source $MISOENV/bin/activate

index_gff --index $REFERENCEDIR/INCOMPLETE_MISSING20/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_missing20.gff3 $MISOINDEXDIR_MISSING20

for n in 1 2 3 4 5 6
do
	miso --run $MISOINDEXDIR_MISSING20 \
	$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}_incompleteannot_missing20/accepted_hits.bam \
	-p 16 --output-dir $NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20/output/sample${n} \
	--read-len 101 --paired-end 175 56
done

## Deactivate the local environment
deactivate

## Summarize miso output
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20/output/' counts_output_dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20'" $RCODEGEN/miso_run.R $ROUT/miso_run_human_nonnull_node_missing20.Rout

## Run DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_miso_missing20' method_name='dexseq_miso_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_node_missing20.Rout

## Subset MISO counts to only those assignable to at least one isoform and apply DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_input_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20' path_to_output_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_assignable_missing20'" $RCODEGEN/subset_miso_counts_assignable.R $ROUT/subset_miso_counts_assignable_human_missing20.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_assignable_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_miso_assignable_missing20' method_name='dexseq_miso_assignable_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_assignable_node_missing20.Rout

## Quantify transcripts instead of transcript combinations and apply DEXSeq
R CMD BATCH --no-restore --no-save "--args path_to_miso_fcn='$RCODEGEN/miso_function.R' miso_output_dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_missing20/output/' counts_output_dir='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_transcripts_missing20'" $RCODEGEN/miso_transcripts_run.R $ROUT/miso_transcripts_run_human_nonnull_node_missing20.Rout
R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_MISSING20/miso_transcripts_missing20' conditions=rep(c('c1','c2'),each=3) flattened_file=NULL out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_MISSING20/dexseq_miso_transcripts_missing20' method_name='dexseq_miso_transcripts_missing20'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_miso_transcripts_node_missing20.Rout

## ------------ Compare complete and incomplete annotations  --------------- ##
## Count the number of reads and bins for each method and complete/incomplete annotations
R CMD BATCH --no-restore --no-save "--args path_to_counts='$NONNULLSIMULATION_NODE/2_counts' output_txt_file='$NONNULLSIMULATION_NODE/2_counts/complete_incomplete_counts.txt'" $RCODEGEN/compare_complete_incomplete.R $ROUT/compare_complete_incomplete_human.Rout

## ========================================================================= ##
##                          ADDITIONAL ANALYSES                              ##
## ========================================================================= ##

## ------------ Compare AS event types between different subsets ----------- ##
R CMD BATCH --no-restore --no-save "--args path_to_astalavista='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' path_to_simulation_design='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' output_pdf='$FIGDIR/compare_as_event_types_nonnull_node.pdf'" $RCODEGEN/compare_as_type_events.R $ROUT/compare_as_type_events_nonnull_node_human.Rout

## ----- Filter out all isoforms with < x% abundance before DEXSeq --------- ##
## Generate new simulation_details file containing also the IsoPct2
R CMD BATCH --no-restore --no-save "--args input_file='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' output_file='$NONNULLSIMULATION_NODE/3_truth/simulation_details_isopct2.txt'" $RCODEGEN/add_isopct2.R $ROUT/add_isopct2.Rout

for v in 5 10 15 25
do
	R CMD BATCH --no-restore --no-save "--args path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' path_to_truth='$NONNULLSIMULATION_NODE/3_truth/simulation_details_isopct2.txt' min_abundance_perc=${v} output_gtf='$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_atleast${v}.gtf' output_missingtx='$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/excluded_transcripts.txt'" $RCODEGEN/generate_incomplete_gtf_perc.R $ROUT/generate_incomplete_gtf_perc_human_${v}.Rout
	python $DEXSEQ/dexseq_prepare_annotation.py --aggregate='no' \
	$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_atleast${v}.gtf \
	$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast${v}.gff
	for n in 1 2 3 4 5 6
	do 
		python $DEXSEQ/dexseq_count.py \
		-p yes -r pos -s no -f bam \
		$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast${v}.gff \
		$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
		$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_ATLEAST${v}/dexseq_nomerge_atleast${v}/dexseq${n}.txt
	done
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_ATLEAST${v}/dexseq_nomerge_atleast${v}' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/INCOMPLETE_ATLEAST${v}/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast${v}.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_ATLEAST${v}/dexseq_htseq_nomerge_atleast${v}' method_name='dexseq_htseq_nomerge_atleast${v}'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_nomerge_atleast${v}.Rout
done

## -- Filter out all isoforms with estimated abundance < x% before DEXSeq -- ##
## Use kallisto for the estimation, and exclude all transcripts with estimated relative abundance below x% in all samples

for v in 5 10 15 25
do
	## Get estimated isoform percentages (TPM) from kallisto results and subset gtf file
	R CMD BATCH --no-restore --no-save "--args path_to_kallisto_output='$NONNULLSIMULATION_NODE/2_counts/kallisto' path_to_conversion_table='$REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt' path_to_isoform_results='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' pattern='sample[0-9]$' path_to_gtf='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' min_abundance_perc=${v} output_gtf='$REFERENCEDIR/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_kallistoest_atleast${v}.gtf' output_missingtx='$REFERENCEDIR/INCOMPLETE_KALLISTOEST/excluded_transcripts_atleast${v}.txt'" $RCODEGEN/kallisto_filter_gtf.R $ROUT/kallisto_filter_gtf_human_${v}.Rout
	python $DEXSEQ/dexseq_prepare_annotation.py --aggregate='no' \
	$REFERENCEDIR/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_kallistoest_atleast${v}.gtf \
	$REFERENCEDIR/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_kallistoest_atleast${v}.gff
	for n in 1 2 3 4 5 6
	do 
		python $DEXSEQ/dexseq_count.py \
		-p yes -r pos -s no -f bam \
		$REFERENCEDIR/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_kallistoest_atleast${v}.gff \
		$NONNULLSIMULATION_NODE/1_reads/tophat/sample${n}/accepted_hits.bam \
		$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast${v}/dexseq${n}.txt
	done
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast${v}' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_kallistoest_atleast${v}.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast${v}' method_name='dexseq_htseq_nomerge_kallistoest_atleast${v}'" $RCODEGEN/dexseq_run.R $ROUT/dexseq_run_human_nonnull_node_nomerge_kallistoest_atleast${v}.Rout
done

## -------------- DEXSeq-noaggreg with counting bin filtering -------------- ##
## Filter by total (normalized) bin count across samples
for v in 39.10969 74.73187 108.93658 174.93327
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_count_${v}' method_name='dexseq_htseq_nomerge_binfilt_count_${v}' filter_bin_count=${v} filter_bin_variance=NULL filter_perc=NULL filter_perc_perkb=NULL filter_perc_pertx=NULL filter_perc_perkb_pertx=NULL" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_count_${v}.Rout
done

## Filter by variance across samples
for v in 0.1500974 0.1773393 0.2069248 0.2545928
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_variance_${v}' method_name='dexseq_htseq_nomerge_binfilt_variance_${v}' filter_bin_count=NULL filter_bin_variance=${v} filter_perc=NULL filter_perc_perkb=NULL filter_perc_pertx=NULL filter_perc_perkb_pertx=NULL" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_variance_${v}.Rout
done

## Filter by relative expression
for v in 0.01506810 0.02045075 0.02355890 0.02811174
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_perc_${v}' method_name='dexseq_htseq_nomerge_binfilt_perc_${v}' filter_bin_count=NULL filter_bin_variance=NULL filter_perc=${v} filter_perc_perkb=NULL filter_perc_pertx=NULL filter_perc_perkb_pertx=NULL" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_perc_${v}.Rout
done

## Filter by relative expression after accounting for bin size
for v in 0.009028732 0.012398071 0.014812590 0.018823718
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_${v}' method_name='dexseq_htseq_nomerge_binfilt_perc_perkb_${v}' filter_bin_count=NULL filter_bin_variance=NULL filter_perc=NULL filter_perc_perkb=${v} filter_perc_pertx=NULL filter_perc_perkb_pertx=NULL" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_perc_perkb_${v}.Rout
done

## Filter by relative expression after accounting for number of covering isoforms
for v in 0.01463351 0.01828132 0.02063856 0.02430609
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_${v}' method_name='dexseq_htseq_nomerge_binfilt_perc_pertx_${v}' filter_bin_count=NULL filter_bin_variance=NULL filter_perc=NULL filter_perc_perkb=NULL filter_perc_pertx=${v} filter_perc_perkb_pertx=NULL" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_perc_pertx_${v}.Rout
done

## Filter by relative expression after accounting for number of covering isoforms and length
for v in 0.01032703 0.01357155 0.01594104 0.01988352
do
	R CMD BATCH --no-restore --no-save "--args path_to_dexseq_fcn='$RCODEGEN/dexseq_function_filter.R' path_to_count_files='$NONNULLSIMULATION_NODE/2_counts/dexseq_nomerge' conditions=rep(c('c1','c2'),each=3) flattened_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' out_basename='$NONNULLSIMULATION_NODE/4_results/FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_${v}' method_name='dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_${v}' filter_bin_count=NULL filter_bin_variance=NULL filter_perc=NULL filter_perc_perkb=NULL filter_perc_pertx=NULL filter_perc_perkb_pertx=${v}" $RCODEGEN/dexseq_run_filter.R $ROUT/dexseq_run_human_nonnull_node_nomerge_binfilt_perc_perkb_pertx_${v}.Rout
done

## ------- Plot coverage of new and eliminated FPs after filtering --------- ##
R CMD BATCH --no-restore --no-save $RCODEGEN/compare_FPs_orig_filt0.05.R $ROUT/compare_FPs_orig_filt0.05.Rout
R CMD BATCH --no-restore --no-save "--args path_to_tophat_output='$NONNULLSIMULATION_NODE/1_reads/tophat' path_to_gtf_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' path_to_second_gtf_file='$REFERENCEDIR/INCOMPLETE_ATLEAST5/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast5.gff' path_to_display_file='$FIGDIR/coverage_plots/FP_filt_0.05_not_orig.txt' highlight_type='exon' conditions=rep(c('c1','c2'),each=3) pattern='sample[0-9]$' output_filename='$FIGDIR/coverage_plots/FP_filt_0.05_not_orig.pdf' plot_height=7 bounding_box=FALSE show_name=FALSE ext_before=0.01 ext_after=0.01" $RCODEGEN/plot_coverage.R $ROUT/plot_coverage_human_new_fps_bins.Rout
R CMD BATCH --no-restore --no-save "--args path_to_tophat_output='$NONNULLSIMULATION_NODE/1_reads/tophat' path_to_gtf_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' path_to_second_gtf_file='$REFERENCEDIR/INCOMPLETE_ATLEAST5/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast5.gff' path_to_display_file='$FIGDIR/coverage_plots/FP_orig_not_filt_0.05.txt' highlight_type='exon' conditions=rep(c('c1','c2'),each=3) pattern='sample[0-9]$' output_filename='$FIGDIR/coverage_plots/FP_orig_not_filt_0.05.pdf' plot_height=7 bounding_box=FALSE show_name=FALSE ext_before=0.01 ext_after=0.01" $RCODEGEN/plot_coverage.R $ROUT/plot_coverage_human_elim_fps_bins.Rout
R CMD BATCH --no-restore --no-save "--args path_to_tophat_output='$NONNULLSIMULATION_NODE/1_reads/tophat' path_to_gtf_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' path_to_second_gtf_file='$REFERENCEDIR/INCOMPLETE_ATLEAST5/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast5.gff' path_to_display_file='$FIGDIR/coverage_plots/FP_filt_0.05_not_orig.txt' highlight_type='exon' conditions=rep(c('c1','c2'),each=3) pattern='sample[0-9]$' output_filename='$FIGDIR/coverage_plots/FP_filt_0.05_not_orig_withbb.pdf' plot_height=7 bounding_box=TRUE show_name=FALSE ext_before=0.01 ext_after=0.01" $RCODEGEN/plot_coverage.R $ROUT/plot_coverage_human_new_fps_bins_withbb.Rout
R CMD BATCH --no-restore --no-save "--args path_to_tophat_output='$NONNULLSIMULATION_NODE/1_reads/tophat' path_to_gtf_file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' path_to_second_gtf_file='$REFERENCEDIR/INCOMPLETE_ATLEAST5/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge_atleast5.gff' path_to_display_file='$FIGDIR/coverage_plots/FP_orig_not_filt_0.05.txt' highlight_type='exon' conditions=rep(c('c1','c2'),each=3) pattern='sample[0-9]$' output_filename='$FIGDIR/coverage_plots/FP_orig_not_filt_0.05_withbb.pdf' plot_height=7 bounding_box=TRUE show_name=FALSE ext_before=0.01 ext_after=0.01" $RCODEGEN/plot_coverage.R $ROUT/plot_coverage_human_elim_fps_bins_withbb.Rout









