#!/bin/bash
### This is a script for making the files needed for submitting dimers with a fixe$
### Can also modify the time request based on the length of the .fasta files 
while read line ; do
    sed "s/\"name\"/"\"${line}\""/g" <alphafold_multimer_with_clean_up.sh >alphafo$
    chmod 777 alphafold_multimer_${line}.sh
    sbatch alphafold_multimer_${line}.sh
done <~/scratch/batch_multimer