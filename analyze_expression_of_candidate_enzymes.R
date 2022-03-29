library(stringr)
library(CEMiTool)
library(DESeq2)
library(edgeR)
library(gdata)
library(ggplot2)
library(reshape2)
library(seqinr)
library(rhmmer)

### 100 aa minimum
blastp_OACs <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/blastp_OACs_heli.txt")
blastp_PTs <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/blastp_PTs_heli.txt")
blastp_TKSs <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/blastp_TKSs_heli.txt")
blastp_AAEs <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/blastp_AAEs_heli.txt")
blastp_BBEs <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/blastp_BBEs_heli.txt")

### change the column names
colnames(blastp_OACs) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blastp_PTs) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blastp_TKSs) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blastp_AAEs) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blastp_BBEs) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

### create a common data frame adding the name of the dataframe
blastp_all <- combine(blastp_OACs, blastp_PTs, blastp_TKSs, blastp_AAEs, blastp_BBEs)

blastp_all$gene_id <- apply(str_split_fixed(blastp_all$sseqid, "\\.", 4)[, 1:2], 1, paste, collapse=".")

gene_ids_all_candidates <- unique(blastp_all[, c("gene_id", "source")])

### read the annotation table
annotation_table <- fread(file = "/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation/10_integration/Humb_PASA_v3_parsed.transcripts.fa_annotation_report.xls", na.strings = ".")

### keep only entries with identified orfs
annotation_table <- as.data.frame(annotation_table[!is.na(annotation_table$prot_id), ])

### remove unused columns
annotation_table <- annotation_table[, colSums(is.na(annotation_table)) < nrow(annotation_table)]

### rename the commented #gene_id
colnames(annotation_table)[1] <- "gene_id"

### read the pfam data
pfam_data <- read_domtblout("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation/3_homology/3_hmmer_Pfam/3_hmmer_Pfam_ALL.txt")
pfam_data <- pfam_data[, c("query_name", "domain_accession", "domain_name")]
pfam_data <- aggregate(list(pfam_data$domain_accession, pfam_data$domain_name), by = list(pfam_data$query_name), toString)
colnames(pfam_data) <- c("query_name", "domain_accession", "domain_name")
annotation_table <- merge(annotation_table, pfam_data, by.x = "prot_id", by.y = "query_name", all.x = T)

### read the structural info
gff_heli <- as.data.frame(readGFF("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/mapping_to_assembly/softmasked/trinity_transcripts/PASA5/Humb_PASA_v3_parsed.gff3"))
gff_heli_genes <- gff_heli[gff_heli$type == "gene", ]
gff_heli_genes <- gff_heli_genes[, c("seqid", "ID", "start", "end", "strand")]
colnames(gff_heli_genes) <- c("scaffold", "gene_id", "start", "end", "strand")

genome_coordinates <- data.frame(gene_id = gff_heli_genes$gene_id,
                               genome_coordinates = as.character(paste0(gff_heli_genes$scaffold, "-", gff_heli_genes$start, ":", gff_heli_genes$end, "(", gff_heli_genes$strand, ")")))

annotation_table <- merge(annotation_table, genome_coordinates, by = "gene_id", all.x = T)

### add the protein sequence
proteins_fasta <- read.fasta(file = "/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation/Humb_PASA_v3_parsed.transcripts.fa.transdecoder.pep", seqtype = "AA", as.string = T)
proteins_fasta_df <- as.data.frame(unlist(proteins_fasta))
proteins_fasta_df$ID <- rownames(proteins_fasta_df)
colnames(proteins_fasta_df)[1] <- "protein_sequence"

annotation_table <- merge(annotation_table, proteins_fasta_df, by.x = "prot_id", by.y = "ID", all.x = T)
annotation_table$protein_sequence <- as.character(annotation_table$protein_sequence)

### add the cds sequence
cds_fasta <- read.fasta(file = "/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation/Humb_PASA_v3_parsed.transcripts.fa.transdecoder.cds", seqtype = "DNA", as.string = T)
cds_fasta_df <- as.data.frame(unlist(cds_fasta))
cds_fasta_df$ID <- rownames(cds_fasta_df)
colnames(cds_fasta_df)[1] <- "cds_sequence"

annotation_table <- merge(annotation_table, cds_fasta_df, by.x = "prot_id", by.y = "ID", all.x = T)
annotation_table$cds_sequence <- as.character(annotation_table$cds_sequence)

### read go data to do the enrichment analysis for the coexpression
go_data <- annotation_table[, c("gene_id", "gene_ontology_BLASTP")]
go_data <- go_data[!duplicated(go_data$`gene_id`), ]
go_data_parsed <- as.data.frame(str_split_fixed(go_data$gene_ontology_BLASTP, pattern = "`", n = Inf))
go_data_parsed$gene <- go_data$`gene_id`
go_data_parsed_long <- melt(go_data_parsed, value.name = "term", id.vars = "gene")[, c("term", "gene")]

### read the trueseq counts table
counts_trueseq <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/mapping_to_assembly/softmasked/short_reads/counts_GFF_pasa_FINAL.txt", header = T)

### remove redundant part of the colnames
colnames(counts_trueseq) <- gsub("_STARAligned.sortedByCoord.out.bam", "", colnames(counts_trueseq), fixed = T)
colnames(counts_trueseq)[7:ncol(counts_trueseq)] <- c("trueseq_youngleaf", "trueseq_oldleaf", "trueseq_stem", "trueseq_root", "trueseq_leafWOtrichomes")

### read the transeq counts table
counts_transeq <- read.table("/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/mapping_to_assembly/softmasked/transeq/counts_deduped_GFF_pasa_FINAL.txt", header = T)

### remove redundant part of the colnames
colnames(counts_transeq) <- gsub("_UMI2Aligned.deduped.bam", "", colnames(counts_transeq), fixed = T)
colnames(counts_transeq) <- gsub("pctgsalsa_softmasked_", "", colnames(counts_transeq), fixed = T)
colnames(counts_transeq)[7:ncol(counts_transeq)] <- paste0("transeq_", colnames(counts_transeq)[7:ncol(counts_transeq)])

counts_all <- cbind(counts_trueseq, counts_transeq[, 7:ncol(counts_transeq)])

counts_cpm_all <- round(cpm(counts_all[, 7:ncol(counts_all)], normalized.lib.sizes=TRUE, log=FALSE), 2)

colnames(counts_cpm_all) <- paste0(colnames(counts_cpm_all), "_cpm")

rownames(counts_cpm_all) <- counts_all$Geneid

annotation_table_with_counts <- merge(annotation_table, counts_all, by.x = "gene_id", by.y = "Geneid", all.x = T)
annotation_table_with_counts <- merge(annotation_table_with_counts, counts_cpm_all, by.x = "gene_id", by.y = 0, all.x = T)


### MAKE THE BLAST HITS MORE INFORMATIVE
### create columns to keep only gene ids
annotation_table_with_counts$cannabis_orth <- sapply(strsplit(annotation_table_with_counts$GCF_blastp_BLASTP, "^", fixed = T), `[`, 1)
annotation_table_with_counts$arabidopsis_orth <- sapply(strsplit(annotation_table_with_counts$arabidopsis_blastp_BLASTP, "^", fixed = T), `[`, 1)
annotation_table_with_counts$sunflower_orth <- sapply(strsplit(annotation_table_with_counts$sunflower_blastp_BLASTP, "^", fixed = T), `[`, 1)
annotation_table_with_counts$tomato_orth <- sapply(strsplit(annotation_table_with_counts$tomato_blastp_BLASTP, "^", fixed = T), `[`, 1)
annotation_table_with_counts$rice_orth <- sapply(strsplit(annotation_table_with_counts$rice_blastp_BLASTP, "^", fixed = T), `[`, 1)

### read and tidy fasta annotation
cannabis <- as.data.frame(str_split_fixed(readLines(file("/home/labs/aharoni/luisdh/uniprot/GCF_900626175.1_cs10_protein.annot2", encoding = "utf-8")), "\t", 2))
cannabis$V1 <- gsub(">", "", cannabis$V1)
colnames(cannabis) <- c("V1", "cannabis_annot")

arabidopsis <- as.data.frame(str_split_fixed(readLines(file("/home/labs/aharoni/luisdh/uniprot/arabidopsis_UP000006548_3702.annot2", encoding = "utf-8")), "\t", 2))
arabidopsis$V1 <- gsub(">", "", arabidopsis$V1)
colnames(arabidopsis) <- c("V1", "arabidopsis_annot")

sunflower <- as.data.frame(str_split_fixed(readLines(file("/home/labs/aharoni/luisdh/uniprot/sunflower_UP000215914_4232.annot2", encoding = "utf-8")), "\t", 2))
sunflower$V1 <- gsub(">", "", sunflower$V1)
colnames(sunflower) <- c("V1", "sunflower_annot")

tomato <- as.data.frame(str_split_fixed(readLines(file("/home/labs/aharoni/luisdh/uniprot/tomato_UP000004994_4081.annot2", encoding = "utf-8")), "\t", 2))
tomato$V1 <- gsub(">", "", tomato$V1)
colnames(tomato) <- c("V1", "tomato_annot")

rice <- as.data.frame(str_split_fixed(readLines(file("/home/labs/aharoni/luisdh/uniprot/rice_UP000059680_39947.annot2", encoding = "utf-8")), "\t", 2))
rice$V1 <- gsub(">", "", rice$V1)
colnames(rice) <- c("V1", "rice_annot")

annotation_table_with_counts <- merge(annotation_table_with_counts, cannabis, by.x = "cannabis_orth", by.y = "V1", all.x = T)
annotation_table_with_counts <- merge(annotation_table_with_counts, arabidopsis, by.x = "arabidopsis_orth", by.y = "V1", all.x = T)
annotation_table_with_counts <- merge(annotation_table_with_counts, sunflower, by.x = "sunflower_orth", by.y = "V1", all.x = T)
annotation_table_with_counts <- merge(annotation_table_with_counts, tomato, by.x = "tomato_orth", by.y = "V1", all.x = T)
annotation_table_with_counts <- merge(annotation_table_with_counts, rice, by.x = "rice_orth", by.y = "V1", all.x = T)

annotation_table_with_counts <- annotation_table_with_counts[order(annotation_table_with_counts$gene_id), ]

all_analyses_modules <- list()

give_me_coexpression <- function(counts, ### THESE COUNTS MUST INCLUDE THE 7 COLUMNS OF THE FEATURECOUNTS
                                 pval = 0.1, 
                                 corr = "pearson", 
                                 dis_threshold = 0.6, 
                                 annotation, 
                                 out_dir_base, 
                                 number_of_analysis){
        out_dir <- paste0(out_dir_base, corr, "_", pval)
        
        cts <- counts[, 7:ncol(counts)]
        rownames(cts) <- counts$Geneid
        
        ### remove low count genes
        low_count_mask <- rowSums(cts) < ncol(cts)
        cts <- cts[!low_count_mask, ]
        
        
        ### create a description of the experiment
        coldata <- data.frame(row.names = colnames(cts), 
                              treatment = str_split_fixed(
                                      str_split_fixed(colnames(cts), "_", 2)[, 2], 
                                      "\\.", 2)[, 1])
        
        ### create dds object
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = ~ treatment)
        
        dds <- DESeq(dds)
        
        ###### choose a normalization 
        vsd <- varianceStabilizingTransformation(dds)
        norm_counts <- as.data.frame(assay(vsd))
        
        ### add sample annotation for cemitools
        sample_annotation <- data.frame(SampleName = rownames(coldata), Class = coldata$treatment)
        
        ### run cemitools
        cem <- cemitool(norm_counts, annot = sample_annotation, min_ngen = 10, force_beta = T, filter_pval = pval, cor_method = corr, diss_thresh = dis_threshold)
        
        diagnostic_report(cem, directory = out_dir, force = T)
        
        modules_genes <- module_genes(cem)
        
        modules_genes_annotated <- merge(modules_genes, annotation_table, by.x = "genes", by.y = "gene_id", all.x = T)
        modules_genes_annotated <- modules_genes_annotated[!duplicated(modules_genes_annotated$genes), ]
        modules_genes_annotated_candidates <- merge(modules_genes_annotated, gene_ids_all_candidates, by.x = "genes", by.y = "gene_id")
        
        # modules_genes_annotated_pfam <- modules_genes_annotated[grepl(paste(pfam_domains,collapse="|"), modules_genes_annotated$domain_name), ]
        
        ### run the enrichment analyses
        cem <- mod_ora(cem, go_data_parsed_long)
        cem <- mod_gsea(cem)
        
        # plot results
        cem <- plot_profile(cem)
        cem <- plot_ora(cem)
        cem <- plot_gsea(cem)
        
        ### get hubs
        hubs <- get_hubs(cem,200)
        hubs_df <- data.frame(module = str_split_fixed(names(unlist(hubs)), "\\.", 2)[, 1], 
                              gene_id = str_split_fixed(names(unlist(hubs)), "\\.", 2)[, 2])
        
        hubs_df <- merge(modules_genes_annotated, hubs_df, by.x = "genes", by.y = 0)
        
        ### generate files
        generate_report(cem, directory = out_dir, force = T)
        write_files(cem, directory = out_dir, force = T)
        save_plots(cem, "all", directory = out_dir, force = T)
        
        all_analyses_modules[[number_of_analysis]] <<- modules_genes_annotated_candidates
        
        norm_counts <- cbind(rownames(norm_counts), norm_counts)
        
        write.xlsx(list("counts" = counts,
                        "counts_transeq_normalized" = norm_counts,
                        "modules" = modules_genes_annotated,
                        "modules_candidates" = modules_genes_annotated_candidates, 
                        "hub_genes" = hubs_df),
                   file = paste0(out_dir, "/modules_annotated.xlsx"), overwrite = T)
        
}

give_me_coexpression(counts = counts_transeq, 
                     pval = 0.1, 
                     corr = "pearson", 
                     dis_threshold = 0.6, 
                     annotation = annotation_table, 
                     number_of_analysis = 1,
                     out_dir_base = "/home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/coexpression/transeq_")

