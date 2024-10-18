#!/opt/anaconda3/bin/env python3

# This program is for analyzing interactions between two sets of proteins using AlphaFold2 that were run on an HPCC (Zaratan @ UMD). 
# Run this program using python3 (v. 3.11.0) after running 'predicted_dimer_chimerax.py' from within chimerax (v. 1.7.1).
# Various plots and tables are generated in four parts as described in the comments below.  
# Antony Jose, Sept 14, 2024. 

# Import modules.
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import json
import math

## Setting variables
experiment='PID-2_IP_Ketting_Lab' # Name of the 'experiment' being performed.
start_dates=['2024_8_13']#, This is the date that the runs were started. 
download_dates=['2024_8_14'] # This is the date that the results were downloaded from the HPCC (Zaratan).
stats_date='2024_8_14' # This is the date when a batch was summarized by running 'predicted_dimer_chimerax.py'.

Prey=['IFD-1','IFB-2','PRDH-1','VHA-4','PID-5','PID-4','KIN-19','PAR-5','T07C4.3','APP-1','PMP-5']
Bait=['PID-2'] 

ranking_threshold = 0.4 # The 0.8*iptm + 0.2*ptm score that will be used as a cutoff for looking at candidate interactions. 
interacting_residues_threshold = 100 # The product of residues constrained in bait protein with those in prey proteins.

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables', exist_ok=True)
path_to_figures_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/'+str(experiment)+'_ranking_'+str(ranking_threshold)+'_interacting_'+str(interacting_residues_threshold)
path_to_tables_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/'+str(experiment)+'_ranking_'+str(ranking_threshold)+'_interacting_'+str(interacting_residues_threshold)
os.makedirs(path_to_figures_output, exist_ok=True)
os.makedirs(path_to_tables_output, exist_ok=True)
path_to_fasta_lists ='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Preparation/data/dimers_to_test/' 

### PART 1: Generating plots and tables that show interaction areas.

df1 = pd.read_csv("/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/alphafold2/summaries/"+str(stats_date)+"_alphafold2_summary_stats", header=None)
df1.columns=['A','B','PAE','distance','area'] ## Summary of run on local HPCC (Zaratan).

## Reorder based on logical sequence of RNA regulators according to current literature.
Prey_order=Prey
regulator_sequence={rs:ix for ix,rs in enumerate(Prey_order)}
df1=df1.sort_values(by='B', key=lambda rs: rs.map(regulator_sequence))

## Calculating and plotting absolute interaction area.

df_all_maxPae_5=df1[df1['PAE'] == 5]
df_all_maxPae_5.to_csv(path_to_tables_output + '/maxPae_5_abs_interaction_area.csv') 

sns.set_theme(rc={'figure.figsize': (7, 10), 'axes.grid': False},style= 'whitegrid', )

min5size = min(df_all_maxPae_5['area'])
max5size = max(df_all_maxPae_5['area'])
scaling=0.3

fig_height=0.5*len(df_all_maxPae_5.columns)
plt.figure(figsize=(3, fig_height))

fig1 = sns.scatterplot(
    data=df_all_maxPae_5,
    x="A",
    y="B",
    size="area",
    hue="PAE",
    palette=["blue"],
    legend=True,
    sizes=(min5size*scaling, max5size*scaling)
)
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("Bait")
plt.ylabel("Prey")
plt.xticks(rotation=90)
fig1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
path_to_final_figure=path_to_figures_output + '/'+ str(experiment)+ '_abs.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')

## Calculating and plotting normalized interaction area.

plt.figure()

A_names=[]
B_names=[]
A_seqs=[]
B_seqs=[]
for date in start_dates:
    with open(path_to_fasta_lists+str(date)+'_A_list') as A_list:
        for line in A_list:
            if line.startswith('>'):
                A_name=line.strip('>').strip('\n')
                A_names.append(A_name)
            else:
                A_seq=len(line.strip('\n'))
                A_seqs.append(A_seq)
    with open(path_to_fasta_lists+str(date)+'_B_list') as B_list:
        for line in B_list:
            if line.startswith('>'):
                B_name=line.strip('>').strip('\n')
                B_names.append(B_name)
            else:
                B_seq=len(line.strip('\n'))
                B_seqs.append(B_seq)

df1_size = pd.DataFrame(list(zip(A_names, A_seqs)),
               columns =['A', 'length'])
df1_size.drop_duplicates(keep='first', inplace=True)
df2_size = pd.DataFrame(list(zip(B_names, B_seqs)),
               columns =['B', 'length'])

df_all_df1=pd.merge(df1, df1_size, on='A')
df_with_length=pd.merge(df_all_df1, df2_size, on='B')
df_with_length.rename(columns={"length_x": "length_A", "length_y": "length_B"}, inplace=True)
df_with_length['NormArea']=df_with_length['area']/(df_with_length['length_A']*df_with_length['length_B'])

# Reorder based on logical sequence of RNA regulators according to current literature.
df_with_length=df_with_length.sort_values(by='B', key=lambda rs: rs.map(regulator_sequence))

df_all_maxPae_5=df_with_length[df_with_length['PAE'] == 5]
df_all_maxPae_5.to_csv(path_to_tables_output + '/maxPae_5_norm_interaction_area.csv') 

sns.set_theme(rc={'figure.figsize': (7, 10), 'axes.grid': False, 'legend.loc': 'best'},style= 'whitegrid')

min5size = min(df_all_maxPae_5['NormArea'])
max5size = max(df_all_maxPae_5['NormArea'])

scaling=30000

fig_height=0.5*len(df_all_maxPae_5.columns)
plt.figure(figsize=(3, fig_height))

fig2 = sns.scatterplot(
    data=df_all_maxPae_5,
    x="A",
    y="B",
    size="NormArea",
    hue="PAE",
    palette=["blue"],
    legend=True,
    sizes=(min5size*scaling, max5size*scaling)
)
sns.move_legend(fig2, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("Bait")
plt.ylabel("Prey")
plt.xticks(rotation=90)
path_to_final_figure=path_to_figures_output  + '/'+ str(experiment)+ '_norm.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')

### PART 2: Generating plots and data tables describing the interactions for all regardless of a positive interactions that pass threshold criteria. 

# Directory for saving files and setting working directory.
os.makedirs(path_to_figures_output+'/interaction_maps', exist_ok=True)
os.makedirs(path_to_figures_output+'/interaction_residue_numbers', exist_ok=True)
os.makedirs(path_to_tables_output+'/interaction_residue_identities', exist_ok=True)

## Path to chimera interaction details
path_to_chimerax_interaction_details_collected='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/'+str(stats_date)+'_chimerax_interaction_details/'
os.makedirs(path_to_chimerax_interaction_details_collected, exist_ok=True)
for date in download_dates:
    path_to_chimerax_interaction_details_current='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/alphafold2/'+str(date)+'/*details' 
    copy_cmd='cp '+str(path_to_chimerax_interaction_details_current)+' '+str(path_to_chimerax_interaction_details_collected)
    os.popen(copy_cmd)

for date in start_dates:
    # Since the fasta file list has strictly 2-line fasta formatting, simply iterating through gives the needed sequence and names. 
    with open(path_to_fasta_lists+str(date)+'_A_list') as A_list:
        for line in A_list:
            if line.startswith('>'):
                A_name=line
                A_name_only=line.strip('>').strip('\n')
                A_count=+1
            else:
                A_seq=line.strip('\n')
                A_count = 0      

                ## Create zero matrix.
                matrix=pd.DataFrame(0, index=range((len(A_seq))), columns=range((len(Prey))))
                aa_num, num_aa, i=[], [], 1
                for aa in A_seq:
                    num_aa=i
                    aa_num.append(num_aa)
                    i +=1
                matrix.columns=Prey
                matrix.index=aa_num

                ## Get needed information from the interaction_details file.  
                constraint_list=[]
                file1=open(path_to_tables_output + '/interaction_residue_identities/'+ str(A_name_only) + '_interacting_residues_all.tsv','a')
                file1.write('Bait\tPrey\tBait_Residues\tPrey_Residues\n')
                for reg in Prey:
                    path_to_residues= str(path_to_chimerax_interaction_details_collected) + str(A_name_only) + '_' + str(reg) + '_pae_5_distance_6_interaction_details'
                    constrained_residues=0
                    A_interactors_final=[]
                    B_interactors_final=[]
                    if os.path.getsize(path_to_residues) > 0:
                        residues_df = pd.read_csv(path_to_residues, sep='\s+', header=None, engine='python', lineterminator='\n')
                        A_interactors=residues_df[0].unique()
                        A_interactors_formatting1=re.sub('\' \'/A:|\'\n \'/A:', ',', str(A_interactors))
                        A_interactors_formatting2=re.sub('\[\'/A:', '', A_interactors_formatting1)
                        A_interactors_final=re.sub('\'\]', '', A_interactors_formatting2).split(',')

                        B_interactors=residues_df[1].unique()
                        B_interactors_formatting1=re.sub('\' \'/B:|\'\n \'/B:', ',', str(B_interactors))
                        B_interactors_formatting2=re.sub('\[\'/B:', '', B_interactors_formatting1)
                        B_interactors_final=re.sub('\'\]', '', B_interactors_formatting2).split(',')

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

                    ## Write out constrained residues in a readable format as a csv.
                    parseA=re.sub('\'| |\[|\]', '', str(A_interactors_final))
                    parseB=re.sub('\'| |\[|\]', '', str(B_interactors_final))
                    file1=open(path_to_tables_output + '/interaction_residue_identities/'+ str(A_name_only) + '_interacting_residues_all.tsv','a')
                    file1.write(str(A_name_only) +'\t'+ str(reg)+'\t'+ str(parseA) + '\t'+ str(parseB) + '\n')
                interaction_table=pd.DataFrame(constraint_list)

                ## Plot a summary of interactions based on numbers of residues constrained.

                regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
                interaction_table=interaction_table.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

                sns.set_theme(rc={'figure.figsize': (3, 7), 'axes.grid': False},style= 'whitegrid', )
                plt.figure()
                fig_height=8*len(interaction_table.columns)
                plt.figure(figsize=(3, fig_height))
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
                plt.ylabel("Prey")
                path_to_final_figure=path_to_figures_output + '/interaction_residue_numbers/' + str(A_name_only)+'_interacting_residues_all.svg'
                plt.savefig(path_to_final_figure, dpi=300, format='svg')

                # Reshape dataframe for plotting
                melted_matrix = (matrix.rename_axis('residue')
                        .reset_index()
                        .melt(id_vars=['residue'], var_name='regulator')
                    )
                # Reorder based on logical sequence of RNA regulators according to current literature.
                regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
                melted_matrix=melted_matrix.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))
                plt.figure()
                fig_height=2*len(melted_matrix.columns)
                plt.figure(figsize=(3, fig_height))
                if max(melted_matrix["value"])>0:
                    palette_var=["white","black"]
                else:
                    palette_var=["white"]

                fig1 = sns.scatterplot(
                    data=melted_matrix,
                    x="residue",
                    y="regulator",
                    size="value",
                    sizes=(40,40),
                    hue="value",
                    marker="s",
                    palette=palette_var,
                    linewidth=0,
                    legend=False,
                )
                plt.xlabel(str(A_name_only)+ ' (amino acid residue)')
                plt.xlim(0,len(A_seq))
                plt.ylabel("Prey")
                #plt.xticks(rotation=90)
                path_to_final_figure=path_to_figures_output + '/interaction_maps/' + str(A_name_only)+'_interaction_maps_all.svg'
                plt.savefig(path_to_final_figure, dpi=300, format='svg')
        file1.close()

### PART 3: Generating plots and data tables describing the interactions for positive interactions alone that pass threshold criteria. 

# Directory for saving files and setting working directory.

for date in start_dates:
    # Since fasta file list has strictly 2-line fasta formatting, simply iterating through gives the needed sequence and names. 
    with open(path_to_fasta_lists+str(date)+'_A_list') as A_list:
        for line in A_list:
            if line.startswith('>'):
                A_name=line
                A_name_only=line.strip('>').strip('\n')
                A_count=+1
            else:
                A_seq=line.strip('\n')
                A_count = 0      

                ## Create zero matrix.
                matrix=pd.DataFrame(0, index=range((len(A_seq))), columns=range((len(Prey))))
                aa_num, num_aa, i=[], [], 1
                for aa in A_seq:
                    num_aa=i
                    aa_num.append(num_aa)
                    i +=1
                matrix.columns=Prey
                matrix.index=aa_num

                ## Get needed information from the interaction_details file.  
                constraint_list=[]
                Prey_order=[]
                file1=open(path_to_tables_output + '/interaction_residue_identities/'+ str(A_name_only) + '_interacting_residues_positives.tsv','a')
                file1.write('Bait\tPrey\tBait_Residues\tPrey_Residues\n')
                for reg in Prey:
                    path_to_residues= str(path_to_chimerax_interaction_details_collected) + str(A_name_only) + '_' + str(reg) + '_pae_5_distance_6_interaction_details'
                    constrained_residues=0
                    A_interactors_final=[]
                    if os.path.getsize(path_to_residues) > 0:
                        residues_df = pd.read_csv(path_to_residues, sep='\s+', header=None, engine='python', lineterminator='\n')
                        A_interactors=residues_df[0].unique()
                        A_interactors_formatting1=re.sub('\' \'/A:|\'\n \'/A:', ',', str(A_interactors))
                        A_interactors_formatting2=re.sub('\[\'/A:', '', A_interactors_formatting1)
                        A_interactors_final=re.sub('\'\]', '', A_interactors_formatting2).split(',')

                        B_interactors=residues_df[1].unique()
                        B_interactors_formatting1=re.sub('\' \'/B:|\'\n \'/B:', ',', str(B_interactors))
                        B_interactors_formatting2=re.sub('\[\'/B:', '', B_interactors_formatting1)
                        B_interactors_final=re.sub('\'\]', '', B_interactors_formatting2).split(',')

                        for residue in A_interactors_final:
                            matrix.loc[int(residue), reg]=1
                        constrained_residues=len(A_interactors_final)
                        Prey_order.append(str(reg))
                        constraint_list.append(
                            {
                                'regulator': reg,
                                'residues': constrained_residues,
                            }
                        )
                        ## Write out constrained residues in a readable format as a csv...prt 1
                        parseA=re.sub('\'| |\[|\]', '', str(A_interactors_final))
                        parseB=re.sub('\'| |\[|\]', '', str(B_interactors_final))
                        file1=open(path_to_tables_output + '/interaction_residue_identities/'+ str(A_name_only) + '_interacting_residues_positives.tsv','a')
                        file1.write(str(A_name_only) +'\t'+ str(reg)+'\t'+ str(parseA) +'\t'+ str(parseB) +'\n')
                interaction_table=pd.DataFrame(constraint_list)

                ## Plot a summary of interactions based on numbers of residues constrained.

                if len(interaction_table)!=0:
                    regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
                    interaction_table=interaction_table.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

                    sns.set_theme(rc={'figure.figsize': (3, 7), 'axes.grid': False},style= 'whitegrid' )
                    plt.figure()
                    fig_height=0.1*len(matrix.columns)
                    plt.figure(figsize=(3, fig_height))
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
                    plt.ylabel("Prey")
                    #plt.xticks(rotation=90)
                    path_to_final_figure=path_to_figures_output + '/interaction_residue_numbers/' + str(A_name_only)+'_interacting_residues_positives.svg'
                    plt.savefig(path_to_final_figure, dpi=300, format='svg')

                    # Drop rows that have only zeroes and reshape dataframe for plotting
                    matrix=matrix.loc[:, (matrix!=0).any(axis=0)]
                    melted_matrix = (matrix.rename_axis('residue')
                            .reset_index()
                            .melt(id_vars=['residue'], var_name='regulator')
                        )

                    regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
                    melted_matrix=melted_matrix.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

                    plt.figure()
                    fig_height=0.2*len(matrix.columns)
                    plt.figure(figsize=(3, fig_height))
                    if max(melted_matrix["value"])>0:
                        palette_var=["white","black"]
                    else:
                        palette_var=["white"]

                    fig1 = sns.scatterplot(
                        data=melted_matrix,
                        x="residue",
                        y="regulator",
                        size="value",
                        sizes=(40,40),
                        hue="value",
                        marker="s",
                        palette=palette_var,
                        linewidth=0,
                        legend=False,
                    )
                    plt.xlabel(str(A_name_only)+ ' (amino acid residue)')
                    plt.xlim(0,len(A_seq))
                    plt.ylabel("Prey")
                    #plt.xticks(rotation=90)
                    path_to_final_figure=path_to_figures_output + '/interaction_maps/' + str(A_name_only)+'_interaction_maps_positives.svg'
                    plt.savefig(path_to_final_figure, dpi=300, format='svg')
        file1.close()

### PART 4: Generating plots that include ranking distributions and final thresholded interactions. 

## Path to chimera interaction details
path_to_rankings_json='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/'+str(stats_date)+'_rankings_json/'
os.makedirs(path_to_rankings_json, exist_ok=True)

## Collect ranking files together.
for date in download_dates:
    for bait in Bait:
        for prey in Prey:
            path_to_rankings_json_current='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Downloads/'+str(experiment)+'/'+str(date)+'/'+str(bait)+str(prey)+'/'+'ranking_debug.json'
            copy_cmd='cp '+str(path_to_rankings_json_current)+' '+str(path_to_rankings_json)+str(bait)+str(prey)+'_ranking_debug.json'
            os.popen(copy_cmd)

## Make dataframe of rankings that are larger than threshold.
df_rankings = pd.DataFrame(index=Prey, columns=Bait)
for bait in Bait:
    for reg in Prey:
        path_to_source='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/'+str(stats_date)+'_rankings_json/'+str(bait)+str(reg)+'_ranking_debug.json'
        if os.path.exists(path_to_source):
            with open(path_to_source) as f:
                ranking = max(json.load(f)["iptm+ptm"].values())
                if ranking > ranking_threshold:
                    df_rankings.loc[reg,bait] = ranking
                else:
                    df_rankings.loc[reg,bait] = 0

## Make dataframe of normalized interaction area.
path_to_data='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/'+str(experiment)+'_ranking_'+str(ranking_threshold)+'_interacting_'+str(interacting_residues_threshold)+'/maxPae_5_norm_interaction_area.csv'
df_area=pd.DataFrame(index=Prey, columns=Bait)
if os.path.exists(path_to_data):
    df = pd.read_csv(path_to_data)
    df_int=pd.concat([df['A'],df['B'],df['NormArea']], axis=1)
    for row in df_int.itertuples():
        df_area.loc[row.B, row.A]=row.NormArea

## Make dataframe of number of constrained residues.
path_to_data='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/'+str(experiment)+'_ranking_'+str(ranking_threshold)+'_interacting_'+str(interacting_residues_threshold)+'/interaction_residue_identities/'
df_residues=pd.DataFrame(columns=['Bait', 'Prey', 'Bait_Residues', 'Prey_Residues'])

for bait in Bait:
    temp=pd.DataFrame(pd.read_csv(str(path_to_data)+str(bait)+'_interacting_residues_all.tsv', sep='\t', skiprows=1, names=['Bait', 'Prey', 'Bait_Residues', 'Prey_Residues']))
    df_residues=pd.concat([df_residues, temp])
    df_residues=df_residues.fillna('')
df_residues['Bait_Residues'] = df_residues.Bait_Residues.map(lambda x: [i.strip() for i in x.split(",")])
df_residues['Prey_Residues'] = df_residues.Prey_Residues.map(lambda x: [i.strip() for i in x.split(",")])
df_residues['count']=df_residues.Bait_Residues.apply(len)*df_residues.Prey_Residues.apply(len) # product of numbers of interacting residues.
df_residues=df_residues[['Bait','Prey','count']]
df_residues=df_residues.loc[df_residues['Bait'] !='Bait']
df_residues.to_csv(path_to_tables_output +'/'+ 'Bait_Prey_Interacting_Residues_Products.csv')
interactors=[] # Collecting together the interactors.
for row in df_residues.itertuples():
    if row.count < interacting_residues_threshold:
        df_rankings.loc[row.Prey, row.Bait]=0 
    else:
        interactors.append(str(row.Bait)+str(row.Prey))
df_rankings.to_csv(path_to_tables_output +'/'+ 'Bait_Prey_Rankings_Thresholded.csv')

## Plot using infomration about areas, rankings, and count of residues constrained.
sns.set_theme(rc={'figure.figsize': (6, 7), 'axes.grid': False}, style= 'whitegrid')

# Reshape dataframe for plotting
melted_matrix1 = (df_rankings.rename_axis('Prey')
        .reset_index()
        .melt(id_vars=['Prey'], var_name='Bait')
    )
melted_matrix1.to_csv(path_to_tables_output +'/'+ 'Melted_Matrix1.csv')
# Reorder based on logical sequence of RNA regulators according to current literature.
regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
melted_matrix1=melted_matrix1.sort_values(by='Prey', key=lambda rs: rs.map(regulator_sequence))

melted_matrix2 = (df_area.rename_axis('Prey')
        .reset_index()
        .melt(id_vars=['Prey'], var_name='Bait')
    )
# Reorder based on logical sequence of RNA regulators according to current literature.
regulator_sequence={rs:ix for ix,rs in enumerate(Prey)}
melted_matrix2=melted_matrix2.sort_values(by='Prey', key=lambda rs: rs.map(regulator_sequence))


melted_matrix1["value2"]=melted_matrix2["value"]
min_size = min(melted_matrix1["value2"])
max_size = max(melted_matrix1["value2"])

scaling=30000
plt.figure()

fig_height=1.5*len(melted_matrix1.columns)
plt.figure(figsize=(3, fig_height))
fig1 = sns.scatterplot(
    data=melted_matrix1,
    y="Prey",
    x="Bait",
    size="value2",
    sizes=(min_size*scaling, max_size*scaling),
    hue="value",
    palette=sns.color_palette('Blues', as_cmap=True),
    linewidth=0,
    legend=True,
)
plt.xticks(rotation=90)
plt.xlabel('Bait')
plt.ylabel('Prey')
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
path_to_final_figure=path_to_figures_output + '/'+ str(experiment)+ '_area_and_rank.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')

# Create plots of ranking scores for all 25 models when one model that satisfies interaction criteria is identified. 
os.makedirs(path_to_figures_output+'/ranking_distributions', exist_ok=True)
plt.figure()
for interactor in interactors:
    path_to_source='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/'+str(stats_date)+'_rankings_json/'+str(interactor)+'_ranking_debug.json'
    path_to_ranking_distributions=path_to_figures_output + '/ranking_distributions/'+ str(interactor)+'_model_ranks.svg'
    if os.path.exists(path_to_source):
        with open(path_to_source) as f:
            ranking_dict = json.load(f)["iptm+ptm"]
            max_rank = max(ranking_dict.values())
            if max_rank > ranking_threshold:
                plt.bar(list(ranking_dict.keys()), ranking_dict.values(), color='grey')
                plt.ylabel('ranking score')
                plt.ylim(0,1)
                plt.xticks(rotation=90)
                plt.savefig(path_to_ranking_distributions, dpi=300, format='svg')
                plt.figure()

## Copy the program with the run parameters as a record.
path_to_program='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/code/predicted_dimer_python.py'
path_to_run_program='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/code/Successful_Runs_and_Parameters/predicted_dimer_python_'+str(experiment)+'_ranking_'+str(ranking_threshold)+'_interacting_'+str(interacting_residues_threshold)+'.py'
copy_cmd='cp '+str(path_to_program)+' '+str(path_to_run_program)
os.popen(copy_cmd)