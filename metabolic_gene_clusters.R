library(ggplot2)
library(gggenes)
library(rtracklayer)
library(Cairo)
library(openxlsx)

gff_heli <- as.data.frame(readGFF("~/heli/pacbio_genome/hifiasm-HiC/mapping_to_assembly/softmasked/trinity_transcripts/PASA5/Humb_PASA_v3_parsed.gff3"))

gff_heli_genes <- gff_heli[gff_heli$type == "gene", ]

gff_heli_genes <- gff_heli_genes[, c("seqid", "ID", "start", "end", "strand")]

colnames(gff_heli_genes) <- c("scaffold", "gene_ID", "start", "end", "strand")

gff_heli_genes$strand <- gsub("+", 1, gff_heli_genes$strand, fixed = T)
gff_heli_genes$strand <- gsub("-", -1, gff_heli_genes$strand, fixed = T)
gff_heli_genes$strand <- as.numeric(gff_heli_genes$strand)

annotation_table_with_counts <- read.xlsx("~/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/Humb_PASA_v3_parsed_functional_annotation_with_counts.xlsx", sheet = 1)

### OLSs
chromstart <- 133927125
chromend <- 134561248

gff_heli_genes_subsetted <- gff_heli_genes[which(gff_heli_genes$scaffold == "scaffold_1" & gff_heli_genes[,"start"] >= chromstart & gff_heli_genes[,"end"] <= chromend),]

genomic_length <- max(gff_heli_genes_subsetted$end) - min(gff_heli_genes_subsetted$start) 

genomic_intervals_low <- c(min(gff_heli_genes_subsetted$start) - 1000,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*1,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*2,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*3,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*4,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*5)

genomic_intervals_high <- genomic_intervals_low + genomic_length/6
genomic_intervals_high[6] <- chromend+1000

gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[1]] <- "scaffold_1_bin1"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[2] & gff_heli_genes_subsetted$start > genomic_intervals_high[1]] <- "scaffold_1_bin2"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[3] & gff_heli_genes_subsetted$start > genomic_intervals_high[2]] <- "scaffold_1_bin3"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[4] & gff_heli_genes_subsetted$start > genomic_intervals_high[3]] <- "scaffold_1_bin4"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[5] & gff_heli_genes_subsetted$start > genomic_intervals_high[4]] <- "scaffold_1_bin5"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[6] & gff_heli_genes_subsetted$start > genomic_intervals_high[5]] <- "scaffold_1_bin6"

gff_heli_genes_subsetted <- merge(gff_heli_genes_subsetted, annotation_table_with_counts[, c('gene_id', "sprot_Top_BLASTP_hit")], by.x = "gene_ID", by.y = "gene_id", all.x = T)

gff_heli_genes_subsetted$class[grepl("CHSY", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <- "Type III PKS"
gff_heli_genes_subsetted$class[grepl("AP1", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <- "AP-1 complex subunit mu-2"
gff_heli_genes_subsetted$class[grepl("ransposon", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <- "Transposon-related protein"
gff_heli_genes_subsetted$class[is.na(gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <- "Protein of unknown function"
gff_heli_genes_subsetted$class[grepl("MYO9", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <- "Myosin 9"
gff_heli_genes_subsetted$class[grepl("SPSA4", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <-  "Sucrose phosphate synthase"
gff_heli_genes_subsetted$class[grepl("RNHX1", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <-  "Ribonuclase H"
gff_heli_genes_subsetted$class[grepl("Y1765", gff_heli_genes_subsetted$sprot_Top_BLASTP_hit)] <-  "LRR receptor"

dummies <- data.frame(scaffold_2 = unique(gff_heli_genes_subsetted$scaffold_2),
                      start = genomic_intervals_low,
                      end = genomic_intervals_high,
                      gene_ID = "any",
                      scaffold = unique(gff_heli_genes_subsetted$scaffold),
                      strand = -1,
                      class = "Myosin 9")

gff_heli_genes_subsetted$class <- factor(gff_heli_genes_subsetted$class, levels = c("AP-1 complex subunit mu-2", "LRR receptor", "Myosin 9", "Protein of unknown function", "Ribonuclase H", "Sucrose phosphate synthase", "Transposon-related protein", "Type III PKS"))

Cairo(file="~/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/clusters/chalcone_synthase_cluster_scaffold1_v3.png",
      type="png",
      units="cm", 
      width=30, 
      height=15, 
      pointsize=12, 
      dpi=300)

ggplot(gff_heli_genes_subsetted, aes(xmin = start, xmax = end, y = scaffold, fill = class, 
                          forward = strand, label = gene_ID)) +
        geom_gene_arrow() +
        # geom_gene_label() +
        geom_blank(data = dummies) +
        ylab("") +
        facet_wrap(~ scaffold_2, ncol = 1, scales = "free") +
        # scale_fill_brewer(palette = "Set1") +
        scale_fill_manual(values = c("#fbf8cc", "#fde4cf", "#ffcfd2", "#f1c0e8", "#cfbaf0", "#a3c4f3", "#90dbf4", "red")) +
        theme_genes()

dev.off()


### PTs scaffold_5:57604810-58602643(-)
### set the coordinates you have to work on
chromstart <- 57604810
chromend <- 58602643
scaffold_id <- "scaffold_5"

### keep only genes from the gff that are within that range
gff_heli_genes_subsetted <- gff_heli_genes[which(gff_heli_genes$scaffold == scaffold_id & gff_heli_genes[,"start"] >= chromstart & gff_heli_genes[,"end"] <= chromend),]

### calculate the length of the genomic interval to plot
genomic_length <- max(gff_heli_genes_subsetted$end) - min(gff_heli_genes_subsetted$start) 

### divide the interval in subintervals to plot one after the other
genomic_intervals_low <- c(min(gff_heli_genes_subsetted$start) - 1000,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*1,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*2,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*3,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*4,
                           min(gff_heli_genes_subsetted$start) + (genomic_length/6)*5)

genomic_intervals_high <- genomic_intervals_low + genomic_length/6
genomic_intervals_high[6] <- chromend+1000

### assign each gene to one interval
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[1]] <- "scaffold_1_bin1"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[2] & gff_heli_genes_subsetted$start > genomic_intervals_high[1]] <- "scaffold_1_bin2"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[3] & gff_heli_genes_subsetted$start > genomic_intervals_high[2]] <- "scaffold_1_bin3"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[4] & gff_heli_genes_subsetted$start > genomic_intervals_high[3]] <- "scaffold_1_bin4"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[5] & gff_heli_genes_subsetted$start > genomic_intervals_high[4]] <- "scaffold_1_bin5"
gff_heli_genes_subsetted$scaffold_2[gff_heli_genes_subsetted$end <= genomic_intervals_high[6] & gff_heli_genes_subsetted$start > genomic_intervals_high[5]] <- "scaffold_1_bin6"

gff_heli_genes_subsetted <- merge(gff_heli_genes_subsetted, annotation_table_with_counts[, c('gene_id', "sunflower_annot")], by.x = "gene_ID", by.y = "gene_id", all.x = T)

gff_heli_genes_subsetted$class[grepl("naringenin", gff_heli_genes_subsetted$sunflower_annot)] <- "UbiA prenyltransferase"
gff_heli_genes_subsetted$class[grepl("Putative transposase", gff_heli_genes_subsetted$sunflower_annot)] <- "Transposon-related protein"
gff_heli_genes_subsetted$class[grepl("Retrotransposon gag", gff_heli_genes_subsetted$sunflower_annot)] <- "Transposon-related protein"
gff_heli_genes_subsetted$class[is.na(gff_heli_genes_subsetted$sunflower_annot)] <- "Protein of unknown function"
gff_heli_genes_subsetted$class[grepl("Reverse transcriptase Ty1/copia-type domain", gff_heli_genes_subsetted$sunflower_annot)] <- "Transposon-related protein"
gff_heli_genes_subsetted$class[grepl("Putative reverse transcriptase", gff_heli_genes_subsetted$sunflower_annot)] <-  "Transposon-related protein"
gff_heli_genes_subsetted$class[grepl("Longin-like domain protein", gff_heli_genes_subsetted$sunflower_annot)] <-  "Longin-domain containing protein"
gff_heli_genes_subsetted$class[grepl("Uncharacterized", gff_heli_genes_subsetted$sunflower_annot)] <-  "Protein of unknown function"
gff_heli_genes_subsetted$class[grepl("Magnesium transporter", gff_heli_genes_subsetted$sunflower_annot)] <-  "Transposon-related protein"
gff_heli_genes_subsetted$class[grepl("Integrase", gff_heli_genes_subsetted$sunflower_annot)] <-  "Transposon-related protein"
gff_heli_genes_subsetted$class[grepl("GRF-type domain-containing", gff_heli_genes_subsetted$sunflower_annot)] <-  "GRF-type protein"

dummies <- data.frame(scaffold_2 = unique(gff_heli_genes_subsetted$scaffold_2),
                      start = genomic_intervals_low,
                      end = genomic_intervals_high,
                      gene_ID = "any",
                      scaffold = unique(gff_heli_genes_subsetted$scaffold),
                      strand = -1,
                      class = "Myosin 9")

Cairo(file="~/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/cannabinoid_genes/clusters/prenyltransferases_cluster_scaffold5.png",
      type="png",
      units="cm", 
      width=30, 
      height=15, 
      pointsize=12, 
      dpi=300)

ggplot(gff_heli_genes_subsetted, aes(xmin = start, xmax = end, y = scaffold, fill = class, 
                                     forward = strand, label = gene_ID)) +
        geom_gene_arrow() +
        # geom_gene_label() +
        geom_blank(data = dummies) +
        ylab("") +
        facet_wrap(~ scaffold_2, ncol = 1, scales = "free") +
        # scale_fill_brewer(palette = "Paired") +
        scale_fill_manual(values = c("#fbf8cc", "#fde4cf", "#ffcfd2", "#f1c0e8", "#cfbaf0", "blue")) +
        theme_genes()

dev.off()
