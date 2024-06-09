#!/bin/bash
### This is a script for making the files needed for submitting dimers with a fixed walltime request (18 hours typically) and space request (1 A100 GPU, 8 CPUs @ 6 GB/CPU).
### Can also modify the time request based on the length of the .fasta files as needed in the file alphafold_multimer.sh
while read line ; do
    sed "s/\"name\"/"\"${line}\""/g" <alphafold_multimer.sh >alphafold_multimer_${line}.sh
    chmod 777 alphafold_multimer_${line}.sh
    sbatch alphafold_multimer_${line}.sh
done <~/scratch/batch_multimer
