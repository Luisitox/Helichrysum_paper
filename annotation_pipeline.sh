### THIS SCRIPT TAKES A FASTA FILE OF TRANSCRIPTS, RUN TRANSDECODER, SPLIT THE TRANSCRIPTS AND PROTEIN FILES AND ANNOTATE WITH THE FASTA BLAST DATABASES SUPLIED IN THE FILE "species_for_blast.txt". IT ALSO RUNS THMM AND SIGNALP.

### ACTIVATE THE ENVIRONMENT WITH ALL THE NEEDED PACKAGES
source ~/anaconda3/etc/profile.d/conda.sh;
conda activate trinity-env;

### this module has to be activated independently
module load signalp/4.1
module load rnammer/1.2
module load tmhmm/2.0c
module load ncbi-blast+
module load hmmer


### PLEASE CHANGE THE FOLLOWING FILES ACCORDINGLY
### FOLDER OF THE ANNOTATION:
DIR="~/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation"

### TRANSCRIPT FILE
FASTA_FILE="~/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/Humb_PASA_v3_parsed.transcripts.fa"

grep '>' $FASTA_FILE | cut -d'>' -f2 | cut -d' ' -f1 > col2
grep '>' $FASTA_FILE | cut -d'>' -f2 | cut -d' ' -f2 | sed 's/gene=//' > col1
paste col1 col2 > /home/labs/aharoni/luisdh/heli/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/gene_trans_map
rm col*

GENE_TRANS_FILE="~/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/gene_trans_map"

### SPECIES TO USE FOR THE ANNOTATION
SPECIES_FOR_BLAST="~/pacbio_genome/hifiasm-HiC/protein_annotation/PASA5/annotation/species_for_blast.txt" ### this file has to be a plain list of the full path to the fasta files for annotation. Blastdb has to be present in the same folder.

NAME=$(basename $FASTA_FILE)

### run trinity stats
TrinityStats.pl $FASTA_FILE > trinityStats.log

### run transdecoder
TransDecoder.LongOrfs -t $FASTA_FILE

### split the transcripts to multiple fasta files
mkdir 1_splitted_transcript_fasta
cd 1_splitted_transcript_fasta

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("splitted_transcript_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $FASTA_FILE

cd ..

### split the proteins fasta file
mkdir 2_splitted_protein_fasta
cd 2_splitted_protein_fasta

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("splitted_longest_orfs_%d.pep",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ../${NAME}.transdecoder_dir/longest_orfs.pep

cd ..

mkdir 3_homology
cd 3_homology

# 1_blastx_uniprot

mkdir 1_blastx_uniprot
cd 1_blastx_uniprot

mkdir logs

DB="/home/labs/aharoni/luisdh/uniprot/uniprot_sprot.pep"
FASTA_DIR="${DIR}/1_splitted_transcript_fasta"

for file in $FASTA_DIR/*.fa
do
        bn=$(basename $file | cut -d. -f1)
        echo $bn
        cmd="'blastx -query $file -db $DB  -max_target_seqs 1 -outfmt 6 -evalue 1e-20 -num_threads 10 -out blastx_uniprot-sprot_${bn}.txt'"
        cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=2000]" -J hom_blastx_${bn} -o logs/blastx_uniprot_${bn}_%J.out -e logs/blastx_uniprot_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

cd ..

# 2_blastp_uniprot

mkdir 2_blastp_uniprot
cd 2_blastp_uniprot

mkdir logs

DB="/home/labs/aharoni/luisdh/uniprot/uniprot_sprot.pep"
FASTA_DIR="${DIR}/2_splitted_protein_fasta"

for file in $FASTA_DIR/*.pep
do
        bn=$(basename $file | cut -d. -f1)
        #echo $bn
        cmd="blastp -query $file -db $DB -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 10 -out blastp_uniprot-sprot_outfmt6_${bn}.txt"
        cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=2000]" -J hom_blastp_${bn} -o logs/blastp_${bn}_%J.out -e logs/blastp_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

cd .. 

# 3_hmmer

mkdir 3_hmmer_Pfam

cd 3_hmmer_Pfam

mkdir logs

DB="/home/labs/aharoni/luisdh/uniprot/Pfam-A.hmm"
FASTA_DIR="${DIR}/2_splitted_protein_fasta"

for file in $FASTA_DIR/*.pep
do
        bn=$(basename $file | cut -d. -f1)
        #echo $bn
        cmd="hmmscan --cpu 10 --domtblout pfam_domtblout_${bn}.txt $DB $file"
        cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=1000]" -J hom_hmmer_${bn} -o logs/hmmer_${bn}_%J.out -e logs/hmmer_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

cd ..


### wait until all bsub commands finish
echo "wait until all bsub commands finish"
bwait -w "ended(hom_*)"  

### generate an ordered files of the partial blasts
ls ${DIR}/3_homology/1_blastx_uniprot/*txt | sort -V > ${DIR}/3_homology/1_blastx_uniprot/ordered_files.list
ls ${DIR}/3_homology/2_blastp_uniprot/*txt | sort -V > ${DIR}/3_homology/2_blastp_uniprot/ordered_files.list
ls ${DIR}/3_homology/3_hmmer_Pfam/*txt | sort -V > ${DIR}/3_homology/3_hmmer_Pfam/ordered_files.list

### concatenate the files according to the order
cat $(cat ${DIR}/3_homology/1_blastx_uniprot/ordered_files.list) > ${DIR}/3_homology/1_blastx_uniprot/1_blastx_uniprot_ALL.txt
cat $(cat ${DIR}/3_homology/2_blastp_uniprot/ordered_files.list) > ${DIR}/3_homology/2_blastp_uniprot/2_blastp_uniprot_ALL.txt
cat $(cat ${DIR}/3_homology/3_hmmer_Pfam/ordered_files.list) > ${DIR}/3_homology/3_hmmer_Pfam/3_hmmer_Pfam_ALL.txt


BLASTX_FILE="${DIR}/3_homology/1_blastx_uniprot/1_blastx_uniprot_ALL.txt"
BLASTP_FILE="${DIR}/3_homology/2_blastp_uniprot/2_blastp_uniprot_ALL.txt"
PFAM_FILE="${DIR}/3_homology/3_hmmer_Pfam/3_hmmer_Pfam_ALL.txt"

cd ..

cmd="TransDecoder.Predict -t $FASTA_FILE --retain_pfam_hits $PFAM_FILE --retain_blastp_hits $BLASTP_FILE"

echo $cmd
eval $cmd

FASTA_FILE_TRANSDECODER=$(pwd)/$(basename $FASTA_FILE).transdecoder.pep

### split the proteins fasta file


### run transdecoder
mkdir 4_splitted_protein_fasta_transdecoder
cd 4_splitted_protein_fasta_transdecoder

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("splitted_transdecoder_predict_%d.pep",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $FASTA_FILE_TRANSDECODER

cd ..

FASTA_DIR="${DIR}/4_splitted_protein_fasta_transdecoder"

mkdir 5_blastp_plants
cd 5_blastp_plants

### continue with the annotation of more plant species
### reads a text file with a plain list to the fasta files that have already a blastp database in the same directory
while read line;
do
        name=$(basename $line | cut -d_ -f1)
        mkdir ${name}_blastp
        cd ${name}_blastp
        DB=$line
        for file in $FASTA_DIR/*.pep
        do
                bn=$(basename $file | cut -d. -f1)
                #eme}cho $bn
                cmd="blastp -query $file -db $DB -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 10 -out blastp_${name}_outfmt6_${bn}.txt"
                cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=2000]" -J blastp_${name}_${bn} -o logs/blastp_${name}_${bn}_%J.out -e logs/blastp_${name}_${bn}_%J.err"
                cmd_eval="${cmd_2} ${cmd}"
                echo $cmd_eval
                eval $cmd_eval
                #sleep 1
        done
        cd ..
done < $SPECIES_FOR_BLAST

echo "wait until blastp finishes"
bwait -w "ended(blastp*)"


for dir in */;
do
        bn=$(basename $dir)
        cd $dir
        ### generate an ordered files of the partial blasts
        ls *txt | sort -V > ordered_files.list

        ### concatenate the files according to the order
        cat $(cat ordered_files.list) > ${bn}_ALL.txt
        echo "${bn}_ALL.txt written"
        cd ..
done

cd ..

# signalP

mkdir 6_signalP

cd 6_signalP

mkdir logs

mkdir temp

for file in $FASTA_DIR/*.pep
do
        bn=$(basename $file | cut -d. -f1)
        #echo $bn
        cmd="signalp -f short -T temp -n signalp_${bn}.out $file"
        cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=2000]" -J signalp_${bn} -o logs/signalp_${bn}_%J.out -e logs/signalp_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

cd ..


# thmm

mkdir 7_thmm

cd 7_thmm

mkdir logs

for file in $FASTA_DIR/*.pep
do
        bn=$(basename $file | cut -d. -f1)
        #echo $bn
        cmd="'tmhmm --short < $file  > tmhmm_${bn}.out'"
        cmd_2="bsub -q "new-short" -n 10 -R "rusage[mem=2000]" -J thmmP_${bn} -o logs/thmm_${bn}_%J.out -e logs/thmm_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

cd ..

mkdir 8_targetP

cd 8_targetP

mkdir logs

for file in $FASTA_DIR/*.pep
do
        bn=$(basename $file | cut -d. -f1)
        #echo $bn
        cmd="/home/labs/aharoni/luisdh/bin/targetp-2.0/bin/targetp -gff3 -org pl -prefix targetP_${bn} -fasta $file"
        cmd_2="bsub -q "new-short" -n 2 -R "rusage[mem=2000]" -J targetP_${bn} -o logs/targetP_${bn}_%J.out -e logs/targetP_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done


### wait until all bsub commands finish
echo "wait until signalP, thmmP and targetP finish"
bwait -w "ended(*P_*)"  

cd ..

for dir in {6..7}*/; 
do
        bn=$(basename $dir)
        cd $dir
        ### generate an ordered files of the partial blasts
        ls *out | sort -V > ordered_files.list
        ### concatenate the files according to the order
        cat $(cat ordered_files.list) > ${bn}_ALL.out
        echo "${bn}_ALL.txt written"
        cd ..
done

cd 8_targetP

ls *gff3 | sort -V > ordered_files.list
cat $(cat ordered_files.list) > 8_targetP_ALL.gff3

cd ..

conda deactivate;
conda activate busco-env;

mkdir 9_busco
cd 9_busco
busco -r -m transcriptome -c 20 -i $FASTA_FILE -o ${NAME}_busco -l embryophyta_odb10

cd ..

conda deactivate;
conda activate trinotate;

mkdir 10_integration
cd 10_integration

gene_trans_map=$GENE_TRANS_FILE
transcript_fasta=$FASTA_FILE
transdecoder_pep=$FASTA_FILE_TRANSDECODER
blastx_uniprot="${DIR}/3_homology/1_blastx_uniprot/1_blastx_uniprot_ALL.txt"
blastp_uniprot="${DIR}/3_homology/2_blastp_uniprot/2_blastp_uniprot_ALL.txt"
hmmer_PFAM="${DIR}/3_homology/3_hmmer_Pfam/3_hmmer_Pfam_ALL.txt"
signalp="${DIR}/6_signalP/6_signalP_ALL.out"
thmm="${DIR}/7_thmm/7_thmm_ALL.out"
targetP="${DIR}/8_targetP/8_targetP_ALL.gff3"

Build_Trinotate_Boilerplate_SQLite_db.pl $NAME

Trinotate ${NAME}.sqlite init --gene_trans_map $gene_trans_map --transcript_fasta $transcript_fasta --transdecoder_pep $transdecoder_pep
Trinotate ${NAME}.sqlite LOAD_swissprot_blastx $blastx_uniprot
Trinotate ${NAME}.sqlite LOAD_swissprot_blastp $blastp_uniprot
Trinotate ${NAME}.sqlite LOAD_pfam $hmmer_PFAM
Trinotate ${NAME}.sqlite LOAD_tmhmm $thmm
Trinotate ${NAME}.sqlite LOAD_signalp $signalp


for dir in ${DIR}/5_blastp_plants/*/;
do
        bn=$(basename $dir)
        blastp_file=${dir}/${bn}_ALL.txt
        Trinotate ${NAME}.sqlite LOAD_custom_blast --outfmt6 $blastp_file --prog blastp --dbtype $bn
done


Trinotate ${NAME}.sqlite report > ${NAME}_annotation_report.xls
