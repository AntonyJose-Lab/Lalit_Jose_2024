#!/usr/bin/python
# encoding=utf8

# This program is for annotating batches of Alphafold3 multimer (DIMER only) runs from from the server (https://golgi.sandbox.google.com/).
# For each dimer, it will filter interactors based on MaxPae and distance. 
# The outputs include all the residues that are under mutual constraint for the interacting partners, and the surface area of interaction.
# A movie with appropriate color highlighting is also generated.
# Antony Jose, May 9, 2024. 

# Import modules.

import os
import re
import pandas as pd
from bs4 import BeautifulSoup
from chimerax.core.commands import run

## Set an array of conditions for which calculations would be performed (maxPae and distance variables)
maxPae_distance=[[5,6],[20,6],[30,2]]
# The distance option limits the number of pseudobonds by only drawing them between pairs of residues with any inter-residue distance ≤ d Å (default 3.0). These pseudobonds are drawn between α-carbons regardless of which atoms were within the distance cutoff. 
# The maxPae option can be used to further limit the pseudobonds to only those representing PAE values ≤ max-error (no default value; if the option is omitted, there is no restriction by PAE).
stats_date='2024_5_18' # This is the date when a batch was summarized.
### Get all the necessary variables from a .csv file that is formatted as follows. [include top row of column labels]
#   A,B,output_date
#   PROTEIN-1,PROTEIN-2,2024_4_24
#   ...
df = pd.read_csv("/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/data/batch_dimer_predictions_alphafold3_google_af2_vs_af3.csv")
for index, row in df.iterrows():
    for m,d in maxPae_distance:
        # variables 
        variables_list=list(row)
        A=variables_list[0]
        B=variables_list[1]
        max_pae=m
        distance=d
        output_date=variables_list[2]
        sequences_name=A+B

        # use the format to identify the name of the model .cif file and the associated scores in the .json file.
        results_directory='/Users/antonyjose/Desktop/JoseLab/Computing/af2_vs_af3/af3/' + str(sequences_name)
        os.chdir(results_directory)

        json_name=os.popen("find . -name \"*full_data_0.json\" | tail -c +3").read()
        print(json_name)
        cif_name=os.popen("find . -name \"*_0.cif\"").read()

        # paths to files
        root_to_input='/Users/antonyjose/Desktop/JoseLab/Computing/af2_vs_af3/af3/'
        root_to_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/analyses/af2_vs_af3/'
        path_to_optional_programs =root_to_output+'code/' ## example: chimera code for flipping the structures when making the movie.
        path_to_cif=root_to_input+str(sequences_name) 
        path_to_json=root_to_input+str(sequences_name)

        # create results directory as needed
        if not os.path.exists(root_to_output + str(output_date)):
            os.makedirs(root_to_output + str(output_date))
        if not os.path.exists(root_to_output+'summaries'):
            os.makedirs(root_to_output+'summaries')

        path_to_residues=root_to_output + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction_details'
        path_to_chimera_log=root_to_output + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_log'
        path_to_summary=root_to_output + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction'
        path_to_all_stats=root_to_output + 'summaries/' + str(stats_date) + '_alphafold3_summary_stats'
        path_to_movie=root_to_output+ str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '.mp4'
        path_to_movie2=root_to_output+str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction_surface.mp4'

        # Set chimeraX commands with paths to files to be defined within run(session, ...).
        open_cif='open ' + path_to_cif + '/'+ cif_name
        open_pae='alphafold pae #1 file ' + path_to_json + '/'+ json_name
        filtered_interactions='alphafold contacts #1/A toResidues #1/B distance ' + str(distance) + ' maxPae ' + str(max_pae) + ' outputFile ' + path_to_residues
        color_A_green = 'color #1/A green'
        color_B_magenta = 'color #1/B magenta'
        bg_color = 'set bgColor white'
        run(session, open_cif)
        run(session, color_A_green)
        run(session, color_B_magenta)
        run(session, bg_color)
        run(session, open_pae)
        run(session, filtered_interactions)
          
        if os.path.getsize(path_to_residues) > 0: ## Process only when there are some residues interacting.

            residues_df = pd.read_csv(path_to_residues, sep='\s', header=None, lineterminator='\n')

            ### Getting all the interacting residues and formatting them for running commands.
            A_interactors=residues_df[0].unique()
            A_interactors_formatting1=re.sub('\' \'/A:|\'\n \'/A:', ',', str(A_interactors))
            A_interactors_formatting2=re.sub('\[\'/A:', '#1/A:', A_interactors_formatting1)
            A_interactors_final=re.sub('\'\]', '', A_interactors_formatting2)
            B_interactors=residues_df[1].unique()
            B_interactors_formatting1=re.sub('\' \'/B:|\'\n \'/B:', ',', str(B_interactors))
            B_interactors_formatting2=re.sub('\[\'/B:', '#1/B:', B_interactors_formatting1)
            B_interactors_final=re.sub('\'\]', '', B_interactors_formatting2)

            ### Prepare the models by annotating and highlighting the model before saving a movie
            add_label_A='2dlabels create ' + A + ' text ' + A + ' color green size 30 xpos 0.05 ypos 0.9'
            if A==B:
                add_label_B='2dlabels create ' + B + '_b' + ' text ' + B + '_b' + ' color magenta size 30 xpos 0.05 ypos 0.85'
            else:
                add_label_B='2dlabels create ' + B + ' text ' + B + ' color magenta size 30 xpos 0.05 ypos 0.85'
            movie_save='movie record; turn y 2 180; wait; movie encode ' + path_to_movie

            ### Optional flipping of the model along the x, y, or z-axis, which might be useful/needed for presentations by convention. Uncomment lines below as needed.
            # open_flip='open ' + path_to_optional_programs + 'flip_structure_along_y.py'
            # run(session, open_flip)
            # flip_it='flip #1'
            # view = 'view #1' 
            # run(session, flip_it)
            # run(session, view)
            ### End of optional flipping code.  

            run(session, add_label_A)
            run(session, add_label_B)
            run(session, movie_save)

            ## Measure the buried surface area.
            buried_area='measure buriedarea '+ A_interactors_final + ' with_atoms2 ' + B_interactors_final + ' probeRadius 0.5'
            run(session, buried_area)
            log_save='log save ' + path_to_chimera_log
            run(session, log_save)

            # Parsing the log file to extract buried surface area.
            with open(path_to_chimera_log, 'r', encoding='utf8') as file:
                html_content = file.read()

            clean_text = '\n'.join(BeautifulSoup(html_content, "html.parser").stripped_strings)
            log_list = re.split(r'\n', clean_text)
            buried_area_line = [item for item in log_list if "Buried" in item]
            interaction_area=re.split('=', str(buried_area_line[0]))[1]

            file1 = open(path_to_summary, 'a')
            file1.write(str(interaction_area) +'\n')
            file1.write(A_interactors_final +'\n')
            file1.write(str(B_interactors_final))
            file1.close()

            file2 = open(path_to_all_stats, 'a')
            results_line=str(A)+','+str(B)+','+str(max_pae)+','+str(distance)+','+str(interaction_area)
            file2.write(str(results_line) + '\n')
            file2.close()

            # Make movie of the interacting residues alone.
            hide_all='hide all cartoons'
            run(session, hide_all)
            interacting_residues='show ' + str(A_interactors_final) + ' ' + str(B_interactors_final) 
            run(session, interacting_residues)
            movie_save2='movie record; turn y 2 180; wait; movie encode ' + path_to_movie2
            run(session, movie_save2)
            delete_labels='2dlabels delete all'
            run(session, delete_labels) ## This is necessary to make new labels for the next pair of interactors.

        if os.path.getsize(path_to_residues) == 0:
            interaction_area='0'
            file2 = open(path_to_all_stats, 'a')
            results_line=str(A)+','+str(B)+','+str(max_pae)+','+str(distance)+','+str(interaction_area)
            file2.write(str(results_line) + '\n')
            file2.close()  
        
        ## A final save of all the steps done using ChimeraX, which will serve as a log for the code that was run.
        log_save='log save ' + path_to_chimera_log
        run(session, log_save)
        log_clear='log clear'
        run(session, log_clear)
        close_model= 'close #1'
        run(session, close_model)
