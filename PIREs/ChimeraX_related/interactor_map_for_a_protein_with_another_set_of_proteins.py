#!/opt/anaconda3/bin/env python3

# This program is for generating a plot of all the amino acids in a PIRE protein (e.g., C55C3.3,C09G5.7,TIMM-17B.2,FBXB-97,W09B7.2,W09B7.1,C08F11.7,RNH-1.3,C38D9.2,F15D4.5,C18D4.6,HIL-4,Y47H10A.5,E01G4.5,Y17D7B.4,K02E2.6) that are constrained by 25 selected regulators of RNA silencing in C. elegans (RDE-4,ADR-2,ERI-1,RDE-1,ERGO-1,PRG-1,HRDE-1,NRDE-3,CSR-1,ALG-2,HRDE-2,RDE-8,RDE-3,MUT-16,MUT-7,RDE-10,DEPS-1,PGL-1,PID-2,SET-25,ZNFX-1,DCR-1,EGO-1,MET-2,NRDE-2). The total number of residues constrained by each interactor are also plotted for each PIRE. The identity of the constrained residues of PIRE for each interaction are written out to a file. 
# Antony Jose, May 14, 2024. 

# Import modules.
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

stats_date='2024_5_9' # This is the date when a batch was summarized.
experiment='predicted_influencers_of_RNA-regulated_expression'

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/pire_interactor_maps', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/pire_interaction_residue_numbers', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/pire_interaction_residue_identities', exist_ok=True)

path_to_figures_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/'+str(experiment)
path_to_tables_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/'+str(experiment)

## Names of all regulators being considered as a list.
af2=['RDE-4','ADR-2','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','MUT-16','PID-2','MUT-7','RDE-3','RDE-8','RDE-10','SET-25']
af3=['DCR-1','ZNFX-1','EGO-1','NRDE-2','MET-2']
all_reg=af2+af3
## Path to fasta files.
root='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Preparation/data'
path_to_fasta_lists =root+'/dimers_to_test/' 

## Path to chimera interaction details
path_to_chimerax_interaction_details='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/2024_5_9_chimerax_interaction_details/'

# Since 2024_5_4_A_list_finished has strictly 2-line fasta formatting, simply iterating through gives the needed sequence and names. 
with open(path_to_fasta_lists+'2024_5_4_A_list_finished') as A_list:
    for line in A_list:
        if line.startswith('>'):
            A_name=line
            A_name_only=line.strip('>').strip('\n')
            A_count=+1
        else:
            A_seq=line.strip('\n')
            A_count = 0      

            ## Create zero matrix.
            matrix=pd.DataFrame(0, index=range((len(A_seq))), columns=range((len(all_reg))))
            aa_num, num_aa, i=[], [], 1
            for aa in A_seq:
                num_aa=i
                aa_num.append(num_aa)
                i +=1
            matrix.columns=all_reg
            matrix.index=aa_num

            #### Get needed information from the interaction_details file.  
            constraint_list=[]
            for reg in af2:
                path_to_residues= str(path_to_chimerax_interaction_details) + str(A_name_only) + '_' + str(reg) + '_pae_5_distance_6_interaction_details'
                constrained_residues=0
                A_interactors_final=[]
                if os.path.getsize(path_to_residues) > 0:
                    residues_df = pd.read_csv(path_to_residues, sep='\s', header=None, lineterminator='\n')
                    A_interactors=residues_df[0].unique()
                    A_interactors_formatting1=re.sub('\' \'/A:|\'\n \'/A:', ',', str(A_interactors))
                    A_interactors_formatting2=re.sub('\[\'/A:', '', A_interactors_formatting1)
                    A_interactors_final=re.sub('\'\]', '', A_interactors_formatting2).split(',')
                    for residue in A_interactors_final:
                        matrix.loc[int(residue), reg]=1
                    constrained_residues=len(A_interactors_final)
                else:
                    constrained_residues=len(A_interactors_final)
                constraint_list.append(
                    {
                        'regulator': reg,
                        'residues': constrained_residues,
                    }
                )
                ### Write out constrained residues in a readable format as a csv...prt 1
                parse=re.sub('\'| |\[|\]', '', str(A_interactors_final))
                file1=open(path_to_tables_output + '/pire_interaction_residue_identities/'+ str(A_name_only) + '_interacting_residues.tsv','a')
                file1.write(str(A_name_only) +'\t'+ str(reg)+'\t'+ str(parse) + '\n')

            for reg in af3:
                constrained_residues=0 # re-initialize to count for AlphaFold3-run regulators
                A_interactors_final=[] # re-initialize to count for AlphaFold3-run regulators
                path_to_residues= str(path_to_chimerax_interaction_details) + str(reg) + '_' + str(A_name_only) + '_pae_5_distance_6_interaction_details'
                if os.path.getsize(path_to_residues) > 0:
                    residues_df = pd.read_csv(path_to_residues, sep='\s', header=None, lineterminator='\n')
                    A_interactors=residues_df[1].unique()
                    A_interactors_formatting1=re.sub('/B:|\'\n \'/B:', '', str(A_interactors))
                    A_interactors_formatting2=re.sub('\'\]|\[\'', '', A_interactors_formatting1)
                    A_interactors_final=re.sub('\' \'', ',', A_interactors_formatting2).split(',')
                    for residue in A_interactors_final:
                        matrix.loc[int(residue), reg]=1
                    constrained_residues=len(A_interactors_final)
                else:
                    constrained_residues=len(A_interactors_final)
                constraint_list.append(
                    {
                        'regulator': reg,
                        'residues': constrained_residues,
                    }
                )
                ### Write out constrained residues in a readable format as a csv...prt 2
                parse=re.sub('\'| |\[|\]', '', str(A_interactors_final))
                file1.write(str(A_name_only) +'\t'+ str(reg)+'\t'+ str(parse) + '\n')
            interaction_table=pd.DataFrame(constraint_list)

            ### Plot a summary of interactions based on numbers of residues constrained.
            RegulatorSeq=['DCR-1','RDE-4','ADR-2','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','ZNFX-1','PID-2','MUT-16','RDE-10','RDE-8','MUT-7','RDE-3','EGO-1','NRDE-2','SET-25','MET-2']
            regulator_sequence={rs:ix for ix,rs in enumerate(RegulatorSeq)}
            interaction_table=interaction_table.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

            sns.set_theme(rc={'figure.figsize': (3, 7), 'axes.grid': False},style= 'whitegrid', )
            plt.figure()
            fig0 = sns.barplot(
                data=interaction_table,
                x="residues",
                y="regulator",
                palette=["black"]
            )
            plt.axvline(x=10, color="grey")
            plt.axvline(x=20, color="blue")
            plt.xlabel(str(A_name_only)+ ' (number aa residues of ' + str(len(A_seq)) + ' total)')
            plt.xlim(0,75)
            plt.ylabel("RNA regulators")
            #plt.xticks(rotation=90)
            path_to_final_figure=path_to_figures_output + '/pire_interaction_residue_numbers/' + str(A_name_only)+'.svg'
            plt.savefig(path_to_final_figure, dpi=300, format='svg')

            # Reshape dataframe for plotting
            melted_matrix = (matrix.rename_axis('residue')
                    .reset_index()
                    .melt(id_vars=['residue'], var_name='regulator')
                )
            #Reorder based on logical sequence of RNA regulators according to current literature.
            RegulatorSeq=['DCR-1','RDE-4','ADR-2','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','MUT-16','RDE-10','ZNFX-1','PID-2','RDE-8','MUT-7','RDE-3','EGO-1','NRDE-2','SET-25','MET-2']
            regulator_sequence={rs:ix for ix,rs in enumerate(RegulatorSeq)}
            melted_matrix=melted_matrix.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

            plt.figure()
            plt.figure(figsize=(3, 7))
            fig1 = sns.scatterplot(
                data=melted_matrix,
                x="residue",
                y="regulator",
                size="value",
                sizes=(40,40),
                hue="value",
                marker="s",
                palette=["white","black"],
                linewidth=0,
                legend=False,
            )
            plt.xlabel(str(A_name_only)+ ' (amino acid residue)')
            plt.xlim(0,len(A_seq))
            plt.ylabel("RNA regulators")
            #plt.xticks(rotation=90)
            path_to_final_figure=path_to_figures_output + '/pire_interactor_maps/' + str(A_name_only)+'.svg'
            plt.savefig(path_to_final_figure, dpi=300, format='svg')
    file1.close()


