# back2ass
Create assembly from BAC short read sequences

## Purpose: 
One step pipeline to clean, filter and assemble short reads into a contiguous bac sequence. 
   
Requirement: 
- SGE cluster
- python 2.xx
- [seqtk](https://github.com/lh3/seqtk) (cleaning bases and adapters)
- [bwa mem](https://github.com/lh3/bwa) (to map reads on vector contaminents)
- [SPAdes](http://bioinf.spbau.ru/spades) (assembler)

## Input: 
Two pair end fastq file (can be gz compressed) 
A fasta file with the vector and possible contaminent sequences

## Usage:
```bash
qsub back2ass.sh Forward.fastq Reverse.fastq vector.fasta quality
```

## Output:
- Various intermediary files in a tmp directory
- mapping results in a dedicated directory
- assembly results in a dedicated directory
