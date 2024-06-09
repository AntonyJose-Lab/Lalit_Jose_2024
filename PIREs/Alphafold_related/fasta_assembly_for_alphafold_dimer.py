#!/usr/bin/env python3

# This program is for taking in two multi-fasta files for the all-by-all pairwise comparisons and generating the multiple putative dimer fasta files for alphafold_multimer.sh. 
# Proteins in list A are tested for potential interactions with proteins in list B.  
# Antony Jose, May 7, 2024.

# Import modules
import os

# Set paths
experiment_date='2024_5_4'
root='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Preparation/data'
path_to_fasta_lists =root+'/dimers_to_test/' 
path_to_list_A=path_to_fasta_lists+str(experiment_date)+'_A_list'
path_to_list_B=path_to_fasta_lists+str(experiment_date)+'_B_list'

# Since both files have strictly 2-line fasta formatting, simply iterating through and combining lines gives the needed files. 

with open(path_to_list_A) as A_list:
    for line in A_list:
        if line.startswith('>'):
            A_name=line
            A_name_only=line.strip('>')
            A_count=+1
        else:
            A_seq=line
            A_count = 0      
            with open(path_to_list_B) as B_list:
                for line in B_list:
                    if line.startswith('>'):
                        B_name=line
                        B_name_only=line.strip('>')
                        B_count=+1
                    else:
                        B_seq=line
                        B_count = 0
                        if not os.path.exists(root+'/'+ 'batch_' +A_name_only.strip('\n')):
                            os.mkdir(root+'/'+ 'batch_' +A_name_only.strip('\n'))
                        file = open(root+'/'+ 'batch_' +A_name_only.strip('\n') + '/'+ A_name_only.strip('\n') + B_name_only.strip('\n') +'.fasta', 'a')
                        file.write(A_name)
                        file.write(A_seq) # + '\n')
                        file.write(B_name)
                        file.write(B_seq) # + '\n')
                        file.close()
