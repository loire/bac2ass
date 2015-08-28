#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q short.q
#$ -V

#set -xv

# Commande : ./progname R1.fastq R.fastq Vector_and_contaminent.fa qual

# R1 and R2 are both raw fastq file (can be compressed with gunzip)
# Vector_and_contaminent.fa is a fastafile with probable source of contamination
# qual is a real number corresponding to quality threshold used by seqtk to remove ambiguous bases in fastq file


# Parsing reads name
arrIN=(${1//_/ })
pbase=${arrIN[0]}
base=$(basename $pbase)



#
## Create temporary directory for storing intermediate files
mkdir -p tmp
#
#
## First, clean with seqtk for low quality bases and for illumina adapters:
#
/home/bin/Seqtk/20150323/x64/trimadap $1 > ./tmp/${base}_trim_R1.fastq
/home/bin/Seqtk/20150323/x64/trimadap $2 > ./tmp/${base}_trim_R2.fastq
/home/bin/Seqtk/20150323/x64/seqtk trimfq -q $4 ./tmp/${base}_trim_R1.fastq > ./tmp/${base}_q$4_R1.fastq
/home/bin/Seqtk/20150323/x64/seqtk trimfq -q $4 ./tmp/${base}_trim_R2.fastq > ./tmp/${base}_q$4_R2.fastq

#echo "Using flash to merge overlapping reads"
#
## First, use flash to merge overlapping reads:
#
#flash -M 200 -t 8 -o  ./tmp/${base}_q$4 ./tmp/${base}_q$4_R1.fastq ./tmp/${base}_q$4_R2.fastq
#
# Then index the bac sequence(s)

echo "Indexing Vector and contaminent sequences"

/home/bin/bwa/last/x64/bwa index $3
mkdir -p indexVector_Cont
mv $3.* indexVector_Cont
cp $3 indexVector_Cont


# Then, map all paired reads on the bac sequence with bwa mem, use samtools to convert to bam, sort by reads name and reconvert to sam file

mkdir -p mapVectorCont

echo "Filtering of Vector and contaminents in reads"

echo ./indexVector_Cont/$(basename $3)
/home/bin/bwa/last/x64/bwa mem -M -t 1 ./indexVector_Cont/$(basename $3)  ./tmp/${base}_q$4_R1.fastq > ./mapVectorCont/tmp_1.sam
samtools view -bS ./mapVectorCont/tmp_1.sam > ./mapVectorCont/tmp_1.bam
samtools sort -n ./mapVectorCont/tmp_1.bam ./mapVectorCont/tmp_1.sorted
samtools view ./mapVectorCont/tmp_1.sorted.bam >  ./mapVectorCont/${base}.paired_1.sam 
rm tmp_1.*
/home/bin/bwa/last/x64/bwa mem -M -t 1 ./indexVector_Cont/$(basename $3)  ./tmp/${base}_q$4_R2.fastq > ./mapVectorCont/tmp_2.sam
samtools view -bS ./mapVectorCont/tmp_2.sam > ./mapVectorCont/tmp_2.bam
samtools sort -n ./mapVectorCont/tmp_2.bam ./mapVectorCont/tmp_2.sorted
samtools view ./mapVectorCont/tmp_2.sorted.bam >  ./mapVectorCont/${base}.paired_2.sam 
rm tmp_2.*
#
## Same with frag reads produced by flash, but without the useless sorting.
#
#/home/bin/bwa/last/x64/bwa mem -M -t 8 ./mapVectorCont/$(basename $3) ./tmp/${base}_q$4.extendedFrags.fastq > ./mapVectorCont/${base}.frag.sam
#
# Now use bwa results to remove vector sequences and output new reads with my own script

./Scripts/extract_BAC_reads.py ./mapVectorCont/${base}.paired_1.sam  
./Scripts/extract_BAC_reads.py ./mapVectorCont/${base}.paired_2.sam  


# Finally launch spades assembler on clean fastq. Results are stored in a dedicated directory. 
# Most probable scaffold can be found by filtering the scaffolds.fasta file by length +  coverage


/home/bin/SPAdes/3.5.0/x64/bin/spades.py -t 1 -k 21,33,55,77,99,127 --careful -1 ./mapVectorCont/${base}.paired_1.vector_filtered.fastq -2 ./mapVectorCont/${base}.paired_2.vector_filtered.fastq -o ${base}_assembly

# Now try to extract most probable scaffolds from results




