library(rtracklayer)
library(HelloRanges)
library(circlize) 
library(ComplexHeatmap)
library(stringr)


### read the selected contigs/scaffolds (small ones removed)
contigs <- read.table("~/heli/data_circos/hifiasm-HiC/selected_contigs.fai")

### read the density of the different features to be plotted 
genes <- read.table("~/heli/data_circos/hifiasm-HiC/genes100K.bedg", col.names = c("chr", "start", "end", "value1"))
TEs <-read.table("~/heli/data_circos/hifiasm-HiC/TEs100K.bedg", col.names = c("chr", "start", "end", "value1"))
True_seq_cov <- read.table("~/heli/data_circos/hifiasm-HiC/TrueSeqCov.bedg", col.names = c("chr", "start", "end", "value1"))
Tran_seq_cov <- read.table("~/heli/data_circos/hifiasm-HiC/TranSeqCov.bedg", col.names = c("chr", "start", "end", "value1"))

### read the position of relevant genes
relevant_genes <- read.table("~/heli/data_circos/hifiasm-HiC/relevant_genes.txt", header = T)

gff_heli <- as.data.frame(readGFF("~/heli/pacbio_genome/hifiasm-HiC/mapping_to_assembly/softmasked/trinity_transcripts/PASA5/Humb_PASA_v3_parsed.gff3"))
gff_heli_genes <- gff_heli[gff_heli$type == "gene", ]
gff_heli_genes <- gff_heli_genes[, c("seqid", "ID", "start", "end", "strand")]
colnames(gff_heli_genes) <- c("scaffold", "gene_ID", "start", "end", "strand")
gff_heli_genes$strand <- gsub("+", 1, gff_heli_genes$strand, fixed = T)
gff_heli_genes$strand <- gsub("-", -1, gff_heli_genes$strand, fixed = T)
gff_heli_genes$strand <- as.numeric(gff_heli_genes$strand)

pos_relevant_genes <- gff_heli_genes[gff_heli_genes$gene_ID %in% relevant_genes$Gene, ]

pos_relevant_genes <- merge(pos_relevant_genes, relevant_genes, by.x = "gene_ID", by.y = "Gene", all.x = T)

pos_relevant_genes <- data.frame("chr" = pos_relevant_genes$scaffold,
                                 "start" = pos_relevant_genes$start,
                                 "end" = pos_relevant_genes$start,
                                 "value1" = 1,
                                 "gene_name" = pos_relevant_genes$Name,
                                 "gene_type" = pos_relevant_genes$Type)
                                 

chromosome_lenghts <- matrix(c(rep(0, nrow(contigs)), contigs$V2), ncol=2)
chromosome_lenghts <- chromosome_lenghts[order(chromosome_lenghts[, 2], decreasing = T), ]
chromosome_names <- paste0("sc", seq(1:nrow(chromosome_lenghts)))

png(filename="~/heli/data_circos/hifiasm-HiC/circlize_v2.png", 
    type="cairo",
    units="cm", 
    width=30, 
    height=30, 
    pointsize=12, 
    res=300)

color1 <- "#d7191c"
color2 <- "#fdae61"
color3 <- "#74add1"
color4 <- "#313695"

circos.par(start.degree = 90)

circos.initialize(factors = unique(genes$chr),
                  xlim = chromosome_lenghts)

### this layer creates the chromosomes
circos.track(ylim=c(0,1), 
             panel.fun=function(x,y) {
                        chr=CELL_META$sector.index
                        xlim=CELL_META$xlim
                        ylim=CELL_META$ylim
                        circos.text(mean(xlim), 
                                    mean(ylim), 
                                    gsub(".*scaffold_", "", CELL_META$sector.index), 
                                    cex=0.8, 
                                    col="grey40",
                                    # facing = "bending.outside",
                                    # facing = "bending", 
                                    facing="bending.inside",
                                    niceFacing=TRUE
                                    )}, 
             bg.col="grey90", 
             bg.border=F, 
             track.height=0.06)

### this layer adds the scale and breakes
brk <- c(0:40)*10^7
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        circos.axis(h="top", 
                    major.at=brk, 
                    labels=round(brk/10^6,1), 
                    labels.cex=0.4,
                    col="grey40",
                    labels.col="grey40",
                    lwd=0.7,
                    labels.facing="clockwise",
                    minor.ticks = 2, major.tick.length = 0.5)
},bg.border=F)


### gene density
cov <- genes
cov$start <- cov$start+1 ### because makewindows overlap in 1
max(cov$value1) ### check the maximum value to make a scale

circos.genomicTrack(data=cov, panel.fun=function(region, value, ...) {
        circos.genomicLines(region, value, type="l", col= color1, border = color1, lwd=0.5, area = T)
        circos.segments(x0=0, x1 = max(contigs$V2), y0=40, y1=40, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1= max(contigs$V2), y0=20, y1=20, lwd=0.6, lty="11", col="grey90")
        # circos.segments(x0=0, x1=nrow(contigs), y0=500, y1=500, lwd=0.6, lty="11", col="grey90")
}, track.height=0.13, bg.border=F, bg.col = "#f7f7f7")

circos.yaxis(at=c(20, 40), labels.cex=0.25, lwd=0, tick.length=0, labels.col="grey40", col="#FFFFFF")


### TEs density
cov <- TEs
cov$start <- cov$start+1 ### because makewindows overlap in 1
max(cov$value1) ### check the maximum value to make a scale

circos.genomicTrack(data=cov, panel.fun=function(region, value, ...) {
        circos.genomicLines(region, value, type="l", col= color2, border = color2, lwd=0.5, area = T)
        # circos.segments(x0=0, x1 = max(contigs$V2), y0=40, y1=40, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=300, y1=300, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=500, y1=500, lwd=0.6, lty="11", col="grey90")
}, track.height=0.13, bg.border=F, bg.col = "#f7f7f7")

circos.yaxis(at=c(0, 300, 500), labels.cex=0.25, lwd=0, tick.length=0, labels.col="grey40", col="#FFFFFF")

# TranSeqCoverage
cov <- Tran_seq_cov
cov$start <- cov$start+1 ### because makewindows overlap in 1
max(cov$value1) ### check the maximum value to make a scale
min(cov$value1)

circos.genomicTrack(data=cov, panel.fun=function(region, value, ...) {
        circos.genomicLines(region, value, type="l", col= color3, border = color3, lwd=0.5, area = T)
        # circos.segments(x0=0, x1=max(contigs$V2), y0=1000, y1=1000, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=5000, y1=5000, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=10000, y1=10000, lwd=0.6, lty="11", col="grey90")
}, ylim=range(cov$value1), track.height=0.1, bg.border=F, bg.col = "#f7f7f7")

# y axis
circos.yaxis(at=c(5000, 10000), labels.cex=0.25, lwd=0, tick.length=0, labels.col="grey40", col="#FFFFFF")

# TrueSeqCoverage
cov <- True_seq_cov
cov$start <- cov$start+1 ### because makewindows overlap in 1
max(cov$value1) ### check the maximum value to make a scale
min(cov$value1)

circos.genomicTrack(data=cov, panel.fun=function(region, value, ...) {
        circos.genomicLines(region, value, type="l", col= color4, border = color4, lwd=0.5, area = T)
        # circos.segments(x0=0, x1=max(contigs$V2), y0=5000, y1=5000, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=10000, y1=10000, lwd=0.6, lty="11", col="grey90")
        circos.segments(x0=0, x1=max(contigs$V2), y0=20000, y1=20000, lwd=0.6, lty="11", col="grey90")
}, ylim=range(cov$value1), track.height=0.1, bg.border=F, bg.col = "#f7f7f7")

# axis
circos.yaxis(at=c(10000, 20000), labels.cex=0.25, lwd=0, tick.length=0, labels.col="grey40", col="#71cd48")

### selected genes
# gene labels
relevant_genes_plot <- pos_relevant_genes[pos_relevant_genes$gene_type %in% c("UGT", "AAT"), ]
# relevant_genes_plot <- pos_relevant_genes

circos.genomicLabels(relevant_genes_plot, 
                     labels.column=5, 
                     cex=0.8, 
                     col="grey40", 
                     line_lwd=0.5, 
                     line_col="grey80", 
                     side="inside", 
                     connection_height=0.04, 
                     labels_height=0.04)

draw.sector(52, 53, rou1 = 1, rou2 = 0.4, clock.wise = F, col = add_transparency("yellow", 0.4), border = NA)
draw.sector(168, 169, rou1 = 1, rou2 = 0.4, clock.wise = F, col = add_transparency("yellow", 0.4), border = NA)


circos.clear()

dev.off()
