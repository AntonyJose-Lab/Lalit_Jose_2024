#!/opt/anaconda3/bin/env python3

import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## Setting variables
ranking_threshold = 0.6
interacting_residues_threshold = 20
## Names of all regulators being considered as a list.
af2=['RDE-4','ADR-2','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','MUT-16','PID-2','MUT-7','RDE-3','RDE-8','RDE-10','SET-25']
af3=['DCR-1','ZNFX-1','EGO-1','NRDE-2','MET-2']
all_reg=af2+af3
URGEs=['C08F11.7','C09G5.7','C18D4.6','C38D9.2','C55C3.3','E01G4.5','F15D4.5','FBXB-97','HIL-4','K02E2.6','RNH-1.3','TIMM-17B.2','W09B7.1','W09B7.2','Y17D7B.4','Y47H10A.5']

os.makedirs('/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/2024_5_9_rankings_json/', exist_ok=True)
path_to_destination='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/reformatted_results/2024_5_9_rankings_json/'
path_to_figures_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/figures/'

for urge in URGEs:
    ### Move and rename the ranking files generated using AF2 run on Zaratan.
    for reg in af2:
        path_to_source='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Downloads/URGEs_to_RNA_regulators/*/'
        copy_cmd='cp '+str(path_to_source)+str(urge)+str(reg)+'/'+'ranking_debug.json'+' '+str(path_to_destination)+str(urge)+str(reg)+'_rankings.json'
        os.popen(copy_cmd)
    ### Move and rename the file with ranking generated using AF3 webserver.
    for reg in af3:
        path_to_source='/Users/antonyjose/Desktop/JoseLab/Computing/AlphaFold3_WebServer/URGEs_to_RNA_regulators/'
        copy_cmd='cp '+str(path_to_source)+str(reg)+str(urge)+'/'+'*_summary_confidences_0.json'+' '+str(path_to_destination)+str(urge)+str(reg)+'_rankings.json'
        os.popen(copy_cmd)

### While all files have been moved and renamed, the .json files for AF2 runs are different from those for AF3 runs. 

### Load up the .json files and extract the maximal values into a dataframe for plotting.
df_rankings = pd.DataFrame(index=all_reg, columns=URGEs)
for urge in URGEs:
    ### Note difference in code for parsing AF2 vs AF3 runs.
    for reg in af2:
        path_to_rankings=str(path_to_destination)+str(urge)+str(reg)+'_rankings.json'
        with open(path_to_rankings) as f:
            ranking = max(json.load(f)["iptm+ptm"].values())
            if ranking > ranking_threshold:
                df_rankings.loc[reg,urge] = ranking
            else:
                df_rankings.loc[reg,urge] = 0


    for reg in af3:
        path_to_rankings=str(path_to_destination)+str(urge)+str(reg)+'_rankings.json'
        with open(path_to_rankings) as f:
            ranking = json.load(f)["ranking_score"]
            if ranking > ranking_threshold:
                df_rankings.loc[reg,urge] = ranking
            else:
                df_rankings.loc[reg,urge] = 0

### Make dataframe of normalized interaction area.
path_to_data='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/pires_at_maxPae_5_norm_interaction_area.csv'
df_area=pd.DataFrame(index=all_reg, columns=URGEs)
df = pd.read_csv(path_to_data)
df_int=pd.concat([df['A'],df['B'],df['NormArea']], axis=1)
for row in df_int.itertuples():
    df_area.loc[row.B, row.A]=row.NormArea

### Make dataframe of number of constrained residues
path_to_data='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/tables/pire_interaction_residue_identities/'
df_residues=pd.DataFrame(columns=['URGE', 'regulator', 'residue'], index=range(0))

for urge in URGEs:
    temp=pd.DataFrame(pd.read_csv(str(path_to_data)+str(urge)+'_interacting_residues.tsv', sep='\t',  header=None))
    temp.columns=['URGE', 'regulator', 'residue']
    df_residues=pd.concat([df_residues, temp])
    df_residues=df_residues.fillna('')
df_residues['residue'] = df_residues.residue.map(lambda x: [i.strip() for i in x.split(",")])
df_residues['count']=df_residues.residue.apply(len)
df_residues=df_residues.drop('residue', axis=1)
df_residues=df_residues[['URGE','regulator','count']]

for row in df_residues.itertuples():
    if row.count < interacting_residues_threshold:
        df_rankings.loc[row.regulator, row.URGE]=0

### Plot using infomration about areas, rankings, and count of residues constrained.
sns.set_theme(rc={'figure.figsize': (6, 7), 'axes.grid': False}, style= 'whitegrid')

# Reshape dataframe for plotting
melted_matrix1 = (df_rankings.rename_axis('regulator')
        .reset_index()
        .melt(id_vars=['regulator'], var_name='URGE')
    )
#Reorder based on logical sequence of RNA regulators according to current literature.
RegulatorSeq=['RDE-4','ADR-2','DCR-1','ERI-1','RDE-1','ERGO-1','PRG-1','ALG-2','CSR-1','HRDE-1','NRDE-3','HRDE-2','PGL-1','DEPS-1','MUT-16','RDE-10','ZNFX-1','PID-2','RDE-8','MUT-7','RDE-3','EGO-1','NRDE-2','SET-25','MET-2']
regulator_sequence={rs:ix for ix,rs in enumerate(RegulatorSeq)}
melted_matrix1=melted_matrix1.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))

melted_matrix2 = (df_area.rename_axis('regulator')
        .reset_index()
        .melt(id_vars=['regulator'], var_name='URGE')
    )
#Reorder based on logical sequence of RNA regulators according to current literature.
regulator_sequence={rs:ix for ix,rs in enumerate(RegulatorSeq)}
melted_matrix2=melted_matrix2.sort_values(by='regulator', key=lambda rs: rs.map(regulator_sequence))


melted_matrix1["value2"]=melted_matrix2["value"]
min_size = min(melted_matrix1["value2"])
max_size = max(melted_matrix1["value2"])
scaling=30000
plt.figure()
fig1 = sns.scatterplot(
    data=melted_matrix1,
    y="regulator",
    x="URGE",
    size="value2",
    sizes=(min_size*scaling, max_size*scaling),
    hue="value",
    palette=sns.color_palette('Blues', as_cmap=True),
    linewidth=0,
    legend=True,
)
plt.xticks(rotation=90)
plt.xlabel('protein encoded by understudied regulated gene')
plt.ylabel('RNA regulators')
sns.move_legend(fig1, "upper left", bbox_to_anchor=(1, 1))
path_to_final_figure=path_to_figures_output + 'pire_interactor_area_and_rank.svg'
plt.savefig(path_to_final_figure, dpi=300, format='svg')
