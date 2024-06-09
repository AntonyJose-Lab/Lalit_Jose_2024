#!/opt/anaconda3/bin/env python3

# This program is for generating a bubble plot of the proteins encoded by URGEs (C55C3.3,C09G5.7,TIMM-17B.2,FBXB-97,W09B7.2,W09B7.1,C08F11.7,RNH-1.3,C38D9.2,F15D4.5,C18D4.6,HIL-4,Y47H10A.5,E01G4.5,Y17D7B.4,K02E2.6. Pseudogenes among the top 25 were ignored) and their interactions with 25 selected regulators of RNA silencing in C. elegans (RDE-4,ADR-2,ERI-1,RDE-1,ERGO-1,PRG-1,HRDE-1,NRDE-3,CSR-1,ALG-2,HRDE-2,RDE-8,RDE-3,MUT-16,MUT-7,RDE-10,DEPS-1,PGL-1,PID-2,SET-25,ZNFX-1,DCR-1,EGO-1,MET-2,NRDE-2). It is essentially generating a single plot where the extent of each interaction is represented by the size of a circle and shaded in three different scaled values that indicate confidence as measured by combinations of PAE and inter-residue distance by AlphaFold.
# Antony Jose, May 10, 2024. 

# Import modules.
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables', exist_ok=True)

path_to_figures_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/pire_interactor_maps_af2_vs_af3_af3_runs'
path_to_tables_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/pire_interaction_residue_identities_af2_vs_af3_af3_runs'

df_all = pd.read_csv("/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/af2_vs_af3/summaries/2024_5_18_alphafold3_summary_stats", header=None)
df_all.columns=['A','B','PAE','distance','area'] ## Run on the Alphafold 3 server.

#Reorder based on logical sequence of RNA regulators according to current literature.

RegulatorSeq=['RDE-4','ADR-2','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','MUT-16','RDE-10','ZNFX-1','PID-2','RDE-8','MUT-7','RDE-3','SET-25']
regulator_sequence={rs:ix for ix,rs in enumerate(RegulatorSeq)}
df_all=df_all.sort_values(by='B', key=lambda rs: rs.map(regulator_sequence))

## Calculating and plotting for absolute interaction area.

df_all_maxPae_5=df_all[df_all['PAE'] == 5]
df_all_maxPae_5.to_csv(path_to_tables_output + '/pires_at_maxPae_5_abs_interaction_area.csv') 
df_all_maxPae_20=df_all[df_all['PAE'] == 20]
df_all_maxPae_20.to_csv(path_to_tables_output + '/pires_at_maxPae_20_abs_interaction_area.csv') 
df_all_maxPae_30=df_all[df_all['PAE'] == 30]
df_all_maxPae_30.to_csv(path_to_tables_output + '/pires_at_maxPae_30_abs_interaction_area.csv') 

sns.set_theme(rc={'figure.figsize': (3, 8), 'axes.grid': False},style= 'whitegrid')

min5size = min(df_all_maxPae_5['area'])
max5size = max(df_all_maxPae_5['area'])
min20size = min(df_all_maxPae_20['area'])
max20size = max(df_all_maxPae_20['area'])
min30size = min(df_all_maxPae_30['area'])
max30size = max(df_all_maxPae_30['area'])
scaling=0.3

fig1 = sns.scatterplot(
    data=df_all_maxPae_30,
    x="A",
    y="B",
    size="area",
    hue="PAE",
    palette=["grey"],
    legend=True,
    sizes=(min30size*scaling, max30size*scaling)
)
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("understudied regulated proteins")
plt.ylabel("RNA regulators")
plt.xticks(rotation=90)

fig1 = sns.scatterplot(
    data=df_all_maxPae_20,
    x="A",
    y="B",
    size="area",
    hue="PAE",
    palette=["orange"],
    legend=True,
    sizes=(min20size*scaling, max20size*scaling)
)
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
fig1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig1 = sns.scatterplot(
    data=df_all_maxPae_5,
    x="A",
    y="B",
    size="area",
    hue="PAE",
    palette=["blue"],
    legend=True,
    sizes=(min20size*scaling, max20size*scaling)
)
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
path_to_final_figure=path_to_figures_output + '/pires_af2_vs_af3_af3_run_abs.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')

plt.figure()

## Calculating and plotting for normalized interaction area.
# import csv files with lengths of proteins.
df1_size=pd.read_csv("/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Preparation/data/dimers_to_test/2024_5_4_A_list_sizes")
df2_size=pd.read_csv("/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Preparation/data/dimers_to_test/2024_5_4_B_list_sizes")

df_all_df1=pd.merge(df_all, df1_size, on='A')
df_with_length=pd.merge(df_all_df1, df2_size, on='B')
df_with_length.rename(columns={"length_x": "length_A", "length_y": "length_B"}, inplace=True)
df_with_length['NormArea']=df_with_length['area']/(df_with_length['length_A']*df_with_length['length_B'])

#Reorder based on logical sequence of RNA regulators according to current literature.
df_with_length=df_with_length.sort_values(by='B', key=lambda rs: rs.map(regulator_sequence))

df_all_maxPae_5=df_with_length[df_with_length['PAE'] == 5]
df_all_maxPae_5.to_csv(path_to_tables_output + '/pires_at_maxPae_5_norm_interaction_area.csv') 
df_all_maxPae_20=df_with_length[df_with_length['PAE'] == 20]
df_all_maxPae_20.to_csv(path_to_tables_output + '/pires_at_maxPae_20_norm_interaction_area.csv') 
df_all_maxPae_30=df_with_length[df_with_length['PAE'] == 30]
df_all_maxPae_30.to_csv(path_to_tables_output + '/pires_at_maxPae_30_norm_interaction_area.csv') 

sns.set_theme(rc={'figure.figsize': (3, 8), 'axes.grid': False, 'legend.loc': 'best'},style= 'whitegrid')

min5size = min(df_all_maxPae_5['NormArea'])
max5size = max(df_all_maxPae_5['NormArea'])
min20size = min(df_all_maxPae_20['NormArea'])
max20size = max(df_all_maxPae_20['NormArea'])
min30size = min(df_all_maxPae_30['NormArea'])
max30size = max(df_all_maxPae_30['NormArea'])
scaling=30000

fig2 = sns.scatterplot(
    data=df_all_maxPae_30,
    x="A",
    y="B",
    size="NormArea",
    hue="PAE",
    palette=["grey"],
    legend=True,
    sizes=(min30size*scaling, max30size*scaling)
)
sns.move_legend(fig2, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("understudied regulated proteins")
plt.ylabel("RNA regulators")
plt.xticks(rotation=90)
fig2 = sns.scatterplot(
    data=df_all_maxPae_20,
    x="A",
    y="B",
    size="NormArea",
    hue="PAE",
    palette=["orange"],
    legend=True,
    sizes=(min20size*scaling, max20size*scaling)
)
sns.move_legend(fig2, "upper left", bbox_to_anchor=(1, 1))
fig2 = sns.scatterplot(
    data=df_all_maxPae_5,
    x="A",
    y="B",
    size="NormArea",
    hue="PAE",
    palette=["blue"],
    legend=True,
    sizes=(min20size*scaling, max20size*scaling)
)
sns.move_legend(fig2, "upper left", bbox_to_anchor=(1, 1))
path_to_final_figure=path_to_figures_output + '/pires_af2_vs_af3_af3_run_norm.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')
