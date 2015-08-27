# bac2ass
Create assembly from BAC short read sequences

## Purpose: 
One step pipeline to clean, filter and assemble short reads into a contiguous bac sequence. 
   
## Requirement: 
- SGE cluster
- python 2.xx
- seqtk(https://github.com/lh3/seqtk) (cleaning bases and adapters)
- bwa mem(https://github.com/lh3/bwa) (to map reads on vector contaminents)
- SPAdes(http://bioinf.spbau.ru/spades) (assembler)


## TODO (yes, you !)

One needs to edit back2ass.sh in order to edit path of the different dependencies. 
So far, the script assumes that you are a CBGP cluster user at INRA, so path are good for you. 
Otherwise, well, change them :)


## Input: 
- Two paired fastq file (can be gz compressed) 
- A fasta file with the vector and possible contaminent sequences

## Usage:
```bash
qsub back2ass.sh Forward.fastq Reverse.fastq vector.fasta qual
```
see intput to make a guess on fastq and fasta origin. 
qual is a floating number use as a threshold when cleaning fastq files for bases accuracy (ex 0.01). see seqtk manual for additional informations.  

## Output:
- Various intermediary files in a tmp directory
- mapping results in a dedicated directory
- assembly results in a dedicated directory
