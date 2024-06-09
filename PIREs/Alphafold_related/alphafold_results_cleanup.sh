#!/bin/bash
## This script is to be run locally from scratch directory (~/scratch) using 'sh alphafold_results_cleanup.sh' after the results of Alphafold multimer have been completed for a set of files.
## After this is run, just the minimal set of files needed for analyses can be downloaded. This also enables larger utilization of the scratch space before the downloading of data becomes necessary.
## This script uses the paired fasta file names provided in a txt file called 'batch_multimer'

cd results
while read line ; do
    cd ${line}
    rm ranked_*
    rm unrelaxed_*
    relaxed_root_name=$(find . -name "relaxed_*" | tail -c +11 | head -c -5)
    mv result_${relaxed_root_name}.pkl _result_${relaxed_root_name}.pkl
    rm result_*
    mv _result_${relaxed_root_name}.pkl result_${relaxed_root_name}.pkl
    cd ..
done <~/scratch/batch_multimer
