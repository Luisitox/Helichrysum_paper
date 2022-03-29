### Genome assembly ###
### Hifiasm command
hifiasm -o heli_hifiasm_HiC -t 40 --h1 $HiC_R1_FASTQ --h2 $HiC_R2_FASTQ $CCS1_FASTQ $CCS2_FASTQ $CCS3_FASTQ
### SALSA command
run_pipeline.py -a $CONTIGS_FASTA -l $CONTIGS_FAI -b $ALIGNMENT_BED -e $RESTRICTION_ENZYMES -o $OUT_FOLDER -m yes
### EDTA command
EDTA.pl --genome $GENOME_FASTA --cds $CDS_FASTA --sensitive 1 --ann 1 -t 30

### Processing of True-Seq Illumina data ###
### Rcorrector to correct short reads
perl run_rcorrector.pl -maxcorK 3 -t 6 -1 $R1_FASTQ -2 $R2_FASTQ
python FilterUncorrectabledPEfastq.py -1 $R1_FASTQ_COR -2 $R2_FASTQ_COR
### Trimgalore trimming command
trim_galore --paired --retain_unpaired --phred33 --output_dir . --length 36 -q 5 --stringency 1 -e 0.1 -j 6 $R1_FASTQ_UNFIXRM $R2_FASTQ_UNFIXRM
### Bowtie2 mapping for rRNA removal
bowtie2 --very-sensitive-local --phred33 -x $rRNA_SILVA -1 $R1_FASTQ_TRIMMED -2 $R2_FASTQ_TRIMMED --threads 12 --met-file metrics.txt --al-conc-gz blacklist_paired_aligned.fq.gz --un-conc-gz rRNA_filtered_paired.fq.gz --al-gz blacklist_unpaired_aligned.fq.gz --un-gz rRNA_filtered_unpaired.fq.gz
### STAR mapping of short reads to Helicrysum genome
STAR --genomeDir $GENOME_DIR --readFilesIn $R1_FASTQ_FILTERED $R2_FASTQ_FILTERED --readFilesCommand zcat --runThreadN 20 --twopassMode Basic --outSAMtype BAM SortedByCoordinate 
### Minimap mapping of long reads to Helicrysum genome
minimap2 -t 30 -ax splice:hq -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $ISOSEQ_FASTA > hq_isoforms.fasta.sam 2> hq_isoforms.fasta.sam.log
### Genome guided Trinity
Trinity --genome_guided_bam $ILLUMINA_BAM --long_reads_bam $ISOSEQ_BAM --genome_guided_max_intron 10000 --output GGtrinity_illumina_isoseq --max_memory 200G --CPU 30
### Genome independent Trinity
Trinity --max_memory 64G --CPU 30 --seqType fq --left $R1_FASTQ_FILTERED --right $R2_FASTQ_FILTERED 

### Iso-seq processing ###
isoseq3 refine --require-polya combined_demux.consensusreadset.xml $INPUT_FOLDER/primers.fasta flnc.bam
isoseq3 cluster flnc.bam polished.bam --verbose --use-qvs

### de novo gene models with braker ###
braker.pl --species=helicrysum --genome=$GENOME --bam=$RNASEQ_BAM,$ISOSEQ_BAM --gff3 --softmasking --addUTR=on --cores 30 

### transcript based gene models with PASA ###
cat $TRINITY_FASTA $TRINITYGG_FASTA $ISOSEQ_FASTA > transcripts.fa
seqclean transcripts.fa
perl accession_extractor.pl < $TRINITY_FASTA > tdn.accs
perl Launch_PASA_pipeline.pl --TRANSDECODER -c $ALIGN_CONFIG -C -R -g $GENOME_FASTA -t transcripts.fa.clean -T -u transcripts.fa --ALIGNERS minimap2 --ALT_SPLICE --TDN tdn.accs --CPU 20 --transcribed_is_aligned_orient
### comprehensive transcriptome
build_comprehensive_transcriptome.dbi -c $ALIGN_CONFIG -t transcripts.fa.clean --min_per_ID 95 --min_per_aligned 30

### merge braker and PASA gene models with evidence modeler ###
partition_EVM_inputs.pl --genome $GENOME_FASTA --gene_predictions $BRAKER_GFF --transcript_alignments $PASA_GFF --segmentSize 1000000 --overlapSize 200000 --partition_listing partitions_list.out
write_EVM_commands.pl --genome $GENOME_FASTA --weights $WEIGHTS_FILE --gene_predictions $BRAKER_GFF --transcript_alignments $PASA_GFF --output_file_name evm.out --partitions partitions_list.out > commands.list
parallel -j 30 < commands.list >& evm.log
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome $GENOME_FASTA
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM_PASA.gff3

### final round of PASA###
Load_Current_Gene_Annotations.dbi -c $ALIGN_CONFIG -g $GENOME_FASTA -P $EVMODELER_GFF
perl Launch_PASA_pipeline.pl -c $ANNOT_CONFIG -A -g $GENOME_FASTA -t transcripts.fa.clean --CPU 10

### UMI aware 3' Tran-seq data analysis ###
### extract UMIs from R2 files
umi_tools extract --bc-pattern=NNNN --stdin=$R2_FASTQ--read2-in=$R1_FASTQ --stdout=$SAMPLE_UMI_FASTQ --read2-stdout
### trimming in two steps
trim_galore -j 8 $SAMPLE_UMI_FASTQ
trim_galore -j 8 --polyA $SAMPLE_UMI_TRIMMED_FASTQ
### mapping 
STAR --genomeDir $GENOME_DIR --readFilesIn $SAMPLE_UMI_TRIMMED2_FASTQ --readFilesCommand zcat --runThreadN 10 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $SAMPLE_UMI
### deduplication
umi_tools dedup --stdin=$SAMPLE_UMI_BAM --log=$SAMPLE_UMI_LOG > $SAMPLE_UMI_DEDUPED
