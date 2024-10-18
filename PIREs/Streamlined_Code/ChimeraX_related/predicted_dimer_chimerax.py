#!/usr/bin/python
# encoding=utf8

# This program is for processing batches of Alphafold multimer (DIMER only) runs from the command line within chimerax.
# It will filter each dimer based on MaxPae and distance. 
# Although many combinations can be tested, a maxPae of 5 and a distance of 6 seem reasonable for identifying interactors. 
# With more experience comparing predictions with experimental results, these criteria can be refined further.
# The outputs include all the residues that are under mutual constraint for the interacting partners, and the surface area of interactions.
# A movie with color highlighting of each protein is also generated.
# Antony Jose, Sept 14, 2024. 

# Import modules.
import os
import re
import pandas as pd
from bs4 import BeautifulSoup
from chimerax.core.commands import run

## Set the PAE allowed and the distance considered for various thresholded calculations.
maxPae_distance=[[5,6]]
# The distance option limits the number of pseudobonds by only drawing them between pairs of residues with any inter-residue distance ≤ d Å. These pseudobonds are drawn between α-carbons regardless of which atoms were within the distance cutoff. 
# The maxPae option can be used to further limit the pseudobonds to only those representing PAE values ≤ max-error.

experiment='PID-2_IP_Ketting_Lab' # Name of the 'experiment' being performed.
download_date='2024_8_14' # This is the date that the results were downloaded from the HPCC (Zaratan)
stats_date='2024_8_14' # This is the date when a batch was summarized by running this program: 'predicted_dimer_chimerax.py'.

## List of the two sets of proteins being compared. 
Prey=['IFD-1','IFB-2','PRDH-1','VHA-4','PID-5','PID-4','KIN-19','PAR-5','T07C4.3','APP-1','PMP-5'] ## Provide in the order that the interactions, if any, would be plotted.
Bait=['PID-2'] 

for m,n in maxPae_distance:
    max_pae,distance = m,n
    for bait in Bait:
        for prey in Prey:
            A=bait
            B=prey
            input_date=download_date
            output_date=download_date
            sequences_name=A+B

            # use the fact that only the highest ranking model was subject to AMBER relaxation to identify the name of the model and the associated scores in the .pkl file.
            results_directory="/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/Downloads/"+str(experiment)+"/"+str(input_date)+"/" + str(sequences_name)
            os.chdir(results_directory)
            final_model_name=os.popen("find . -name \"relaxed_*\" | tail -c +11 | head -c +26").read()
            pdb_path='/'+str(sequences_name)+'/relaxed_'+ str(final_model_name)+'.pdb' 
            pkl_path='/'+str(sequences_name)+'/result_'+str(final_model_name)+'.pkl' 

            # paths to files
            root_to_input='/Users/antonyjose/Desktop/JoseLab/Computing/Zaratan/'
            root_to_output='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/'
            path_to_optional_programs =root_to_output+'code/' ## example: chimera code for flipping the structures when making the movie.

            path_to_pdb=root_to_input+'Downloads/'+str(experiment)+'/'+str(input_date)+str(pdb_path)
            path_to_pae=root_to_input+'Downloads/'+str(experiment)+'/'+str(input_date)+str(pkl_path)

            # create results directories as needed
            if not os.path.exists(root_to_output+'analyses/alphafold2/' + str(output_date)):
                os.makedirs(root_to_output+'analyses/alphafold2/' + str(output_date))
            if not os.path.exists(root_to_output+'analyses/alphafold2/summaries'):
                os.makedirs(root_to_output+'analyses/alphafold2/summaries')

            path_to_residues=root_to_output+'analyses/alphafold2/' + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction_details'
            path_to_chimera_log=root_to_output+'analyses/alphafold2/' + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_log'
            path_to_summary=root_to_output+'analyses/alphafold2/' + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction'
            path_to_all_stats=root_to_output + 'analyses/alphafold2/summaries/' + str(stats_date) + '_alphafold2_summary_stats'
            path_to_movie=root_to_output+'analyses/alphafold2/' + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '.mp4'
            path_to_movie2=root_to_output+'analyses/alphafold2/' + str(output_date) + '/' + str(A) + '_' + str(B) + '_pae_'+ str(max_pae) + '_distance_' + str(distance) + '_interaction_surface.mp4'

            # Set chimeraX commands with paths to files to be defined within run(session, ...).
            open_pdb='open ' + path_to_pdb
            open_pae='alphafold pae #1 file ' + path_to_pae
            filtered_interactions='alphafold contacts #1/A toResidues #1/B distance ' + str(distance) + ' maxPae ' + str(max_pae) + ' outputFile ' + path_to_residues
            color_A_green = 'color #1/A green'
            color_B_magenta = 'color #1/B magenta'
            bg_color = 'set bgColor white'
            run(session, open_pdb)
            run(session, color_A_green)
            run(session, color_B_magenta)
            run(session, bg_color)
            run(session, open_pae)
            run(session, filtered_interactions)
            
            if os.path.getsize(path_to_residues) > 0: ## Process only when there are some residues interacting.

                residues_df = pd.read_csv(path_to_residues, sep='\s', header=None, lineterminator='\n')

                ### Get all the interacting residues and format them for running commands.
                A_interactors=residues_df[0].unique()
                A_interactors_formatting1=re.sub('\' \'/A:|\'\n \'/A:', ',', str(A_interactors))
                A_interactors_formatting2=re.sub('\[\'/A:', '#1/A:', A_interactors_formatting1)
                A_interactors_final=re.sub('\'\]', '', A_interactors_formatting2)
                B_interactors=residues_df[1].unique()
                B_interactors_formatting1=re.sub('\' \'/B:|\'\n \'/B:', ',', str(B_interactors))
                B_interactors_formatting2=re.sub('\[\'/B:', '#1/B:', B_interactors_formatting1)
                B_interactors_final=re.sub('\'\]', '', B_interactors_formatting2)

                ### Prepare the models by annotating and highlighting the model before saving a movie.
                add_label_A='2dlabels create ' + A + ' text ' + A + ' color green size 30 xpos 0.05 ypos 0.9'
                if A==B:
                    add_label_B='2dlabels create ' + B + '_b' + ' text ' + B + '_b' + ' color magenta size 30 xpos 0.05 ypos 0.85'
                else:
                    add_label_B='2dlabels create ' + B + ' text ' + B + ' color magenta size 30 xpos 0.05 ypos 0.85'
                movie_save='movie record; turn y 2 180; wait; movie encode ' + path_to_movie

                # ### Optional flipping of the model along the x, y, or z-axis, which might be useful/needed for presentations by convention. Uncomment lines below as needed.
                # open_flip='open ' + path_to_optional_programs + 'flip_structure_along_y.py'
                # run(session, open_flip)
                # flip_it='flip #1'
                # view = 'view #1' 
                # run(session, flip_it)
                # run(session, view)
                # ### End of optional flipping code.  

                run(session, add_label_A)
                run(session, add_label_B)
                run(session, movie_save)

                ## Measure the buried surface area.
                buried_area='measure buriedarea '+ A_interactors_final + ' with_atoms2 ' + B_interactors_final + ' probeRadius 0.5'
                run(session, buried_area)

                ## Save the log file to extract measurements.
                log_save='log save ' + path_to_chimera_log
                run(session, log_save)

                ## Parse the log file to extract buried surface area. Can also use this approach to extract any information logged on ChimeraX.
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

                ## Make movie of the interacting residues alone.
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

## Copy the program with the run parameters as a record.
path_to_program='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/code/predicted_dimer_chimerax.py'
path_to_run_program='/Users/antonyjose/Desktop/JoseLab/Computing/ChimeraX/code/Successful_Runs_and_Parameters/predicted_dimer_chimerax_'+str(experiment)+'.py'
copy_cmd='cp '+str(path_to_program)+' '+str(path_to_run_program)
os.popen(copy_cmd)