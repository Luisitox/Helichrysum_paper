# Helichrysum paper
## Commands used to analyze the data of Helichrysum genome

1) **general_commands.sh** is a shell file including all the commands run to do the assembly and annotate the genome. It is a guideline of the parameters used at each stage and the full process is described in the Methods section of the paper.

2) **annotation_pipeline.sh** is a shell file with the annotation pipeline used to annotate the proteins taking as input a transcripts file. It was written to work in an LSF environment with certain modules and conda environments.

3) **analyze_expression_of_candidate_enzymes.R** is an R script that reads the RNAseq counts, the annotation files obtained in **1.** and **2.**, and blastp searches with query candidate enzymes and finaly creates and run a function for counts normalization and coexpression analysis using DEseq2 and CEMItool respectively. It outputs, besides the typical outputs of CEMItool, an excel file filtering the genes according to the candidate enzymes.

4) **circlize_plot.R** is an R script that reads the density of genes, the density of TEs, the coverage of the RNAseq, and the genomic position of candidate genes and make a circos plot with circlize package.

5) **metabolic_gene_cluster.R** is an R script that takes the structural and functional annotation and makes plots of the genes present in the selected genomic coordinates.
