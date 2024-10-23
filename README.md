# Selecting genes for analysis using historically contingent progress: from RNA changes to protein-protein interactions

Farhaan Lalit<sup>1</sup> and Antony M. Jose<sup>1,*</sup>

<sup>1</sup>Department of Cell Biology and Molecular Genetics, University of Maryland, College Park, USA.

<sup>*</sup>Corresponding author email: amjose@umd.edu

This paper analyzes data from the field of RNA silencing in <i>C. elegans</i> to identify understudied regulated genes (URGEs) and uses AlphaFold to analyze the proteins encoded by them to reveal predicted influencers of RNA-regulated expression (PIREs). 

<b>URGEs</b>: Data tables were reformatted manually and/or using custom scripts, and filtered to generate lists of significant genes. The top ‘g’ genes that occur in the greatest numbers of tables were culled as the most frequently identified genes. A measure for the extent of regulation of each gene (<i>r<sub>g</sub></i>) was used to aid prioritization for detailed study. Co-occurrence patterns of genes in different tables were captured using the Jaccard distance (<i>d<sub>J</sub></i>) or as a symmetric measure of normalized mutual information, defined as Historical Mutual Information (HMI). The <i>d<sub>J</sub></i> values were used to generate a dendrogram using the average linkage method. HMI was used to group genes into clusters according to the Girvan-Newman algorithm and different sets of genes were highlighted. 

<b>PIREs</b>: Proteins encoded by the top 25 URGEs were tested for predicted interactions with 25 key regulators of RNA silencing using AlphaFold 2 and/or AlphaFold 3. The resulting models were analyzed and filtered using various criteria to define PIRE proteins. Similar analyses performed on the proteins encoded by some of the top 100 genes, RDE-3 and its interactors, STAU-1 and its interactors, PID-2 and its co-immunoprecipitating proteins were included when revising the manuscript in response to reviews. Additional analyses performed using more streamlined code included plotting of distributions of model rank scores.     

The code used, key datasets, and supporting files are available here. The streamlined code developed during revisions is presented in a separate section below. Please change paths as necessary.

## Project Structure

- `URGEs/`  
    - `code_and_resources/`: Directory containing all the code used for analyses and visualizations.
        - `data/`: Directory containing subdirectories to be used with the code relevant to standardizing files. These subdirectories are empty. Users should put their cleaned CSV files into the `cleaned/` folder to start with. After running `0_sanitization.py`, files will automatically get added to the `sanitized/` folder. Similarly, after running `0_filter_pvals.py`, files will automatically get added to the `filtered/` folder. 
            - `cleaned/`: Directory containing all the datasets in CSV format after cleaning.
            - `filtered/`: Directory containing all the datasets after being filtered for p-value < 0.05, where applicable.
            - `sanitized/`: Directory containing all the data files after gene name sanitization.
        - `0_rg_simulation_box-whisker.ipynb`: A simulation looking at the relationship between S<sub>i</sub>, T<sub>i</sub>, and g.
        - `0_dataset_wormexp_overlap.ipynb`: Code used to determine how many datasets chosen for this study that were published until 2017 (Extended_Data_Table_until_2017) overlap with the 461 included in WormExp 2.0 as on 27 Jul 2017.
        - `0_filter_pvals.py`: Code used to filter datasets for p-value < 0.05.
        - `0_sanitization.py`: Code used to sanitize gene names.
        - `1_TableOccupancy.py`: Code used to create a TableOccupancy file that counts how many datasets each gene occurs in.
        - `2_fig1D_r25_references.ipynb`: Code used to generate Figure 1D. 
        - `3_fig1F_fig1G_r25_related.ipynb`: Code used to generate Figure 1F and Figure 1G.
        - `4_heatmaps_normalized_100_full.ipynb`: Code used to generate heatmaps showing presence of genes across datasets for the top 100 genes.
        - `5_sklearn_nmi.ipynb`: Code used to cluster genes and create a distance matrix based on normalized mutual information between genes.
        - `6_HMI_explorer.py`: App for generating the graphical user interface for exploring genes clustered based on historical mutual information. The app can be run locally by following these steps:
          1. Download the program to your local machine.
          2. Set download folder as the working directory.
          3. In the terminal, run `python 6_HMI_explorer.py`.
          4. Open a browser window and navigate to the URL given in the terminal.
        - `7_fig8.ipynb`: Code used to generate Figure 8.
        - `8_table_S5_r25_with_100_total.ipynb`
        - `8_table_S5_r25_with_1000_total.ipynb`
        - `2016_Ni_et_al_Epigenetics_and_Chromatin_hrde1_dep_genes_15C_common_ONLY`: A list of HRDE-1-dependent genes.
        - `TableOccupancy_full.csv`: CSV file listing how many datasets each gene occurs in.
        - `df_na_full.csv`: A CSV file that shows gene presence and absence in datasets.
        - `communities.csv `: A file containing clustered communities after running the app.
        - `community_graph.pickle`: File containing graph data for generation of Figure 8.
        - `dist_matrix_100_full.csv`: Distance matrix used to generate an adjacency matrix for creating the network.
        - `Extended_Data_Table_until_2017.xlsx`: An Excel file containing a list of the studies until 2017 that were used to compare with those in WormExp.
        - `simplemine_results_TableOccupancy_full_27Sep2024.xlsx`: Simplemine results for all genes in occupancy table downloaded on September 27, 2024.
        - `Wormbase_Gene_Sanitizer_Database_DownloadedOn_10-30-2023.txt`: Wormbase Gene Sanitizer Database downloaded on October 30, 2023.
        - `WormExp_2_0_info_20170727_Unique_Studies.xlsx`: An Excel file containing URLs for the 461 datasets included in WormExp 2.0 as on 27 Jan 2023.
        - `r_100_simplemine_results_TableOccupancy_full_27Sep2024.xlsx`: Annotated table highlighting particular subsets of genes.
    - `data_full`: Directory containing all the formatted datasets used in the study as .csv files.

- `PIREs/`
    - `Alphafold_related`: Directory containing all the scripts used to run Alphafold and process the results on HPCC.
        - `fasta_assembly_for_alphafold_dimer.py`: This script assembles pairs of fasta files from two lists of multi-fasta files for testing interactions between them using AlphaFold.
        - `run_alphafold.py`: This script runs alphafold v. 2.3.2. It was set up and run as instructed in https://github.com/google-deepmind/alphafold.
        - `alphafold_multimer.sh`: This script calls `run_alphafold.py` for running on the HPCC (Zaratan at UMD). 
        - `alphafold_multimer_batch_submission.sh`: This script calls `alphafold_multimer.sh` and loops through the list of protein pairs in `batch_multimer` (see `Support_File_examples` below) to submit multiple jobs at once to the HPCC.
        - `alphafold_results_cleanup.sh`: This script deletes all models except the top ranking model to save storage space.
    - `ChimeraX_related`: Directory containing all the scripts used for processing results using ChimeraX and for analyzing the results.
        - `alphafold2_dimer_batch_computed_on_zaratan.py`: This script analyzes the models predicted by AlphaFold 2. It extracts interaction areas between the proteins for a variety of criteria and generates annotated movies as needed using ChimeraX 1.7.1.
        - `alphafold3_dimer_batch_computed_on_google.py`: This script is very similar to `alphafold2_dimer_batch_computed_on_zaratan.py` and is used for analyzing models predicted by the AlphaFold 3 server. 
        - `alphafold3_dimer_batch_computed_on_google_comparing_af2_af3.py`: This script is very similar to `alphafold2_dimer_batch_computed_on_zaratan.py` and is used for comparing structures predicted by AlphaFold 2 with those predicted by AlphaFold 3.
        - `predicted_influencer_of_RNA_regulated_expression_d2.py`: This script is used to analyze the summary data obtained for all interactions (e.g., as in `2024_5_9_alphafold2_summary_stats`) and other outputs of `alphafold2_dimer_batch_computed_on_zaratan.py`  to generate plots, and to write out the residues that are constrained as per the criteria for interaction used (e.g., maxPae <5 and distance <6)
        - `predicted_influencer_of_RNA_regulated_expression_d2_af2_vs_af3_af3_run.py`: This script is very similar to `predicted_influencer_of_RNA_regulated_expression_d2.py` and is used for comparing structures run on AlphaFold 2 with results from AlphaFold 3.
        - `interactor_map_for_a_protein_with_another_set_of_proteins.py`: This script is used to generate plots that show the interacting residues along a 'bait' protein, generate plots for the total numbers of residues constrained by each interaction based on the criteria used in `predicted_influencer_of_RNA_regulated_expression_d2.py`, and write out tables with interacting residues of the 'bait' protein for each interactor.
        - `interactor_map_for_a_protein_with_another_set_of_proteins_comparing_af2_vs_af3_rerun_on_af3.py`: This script is very similar to `interactor_map_for_a_protein_with_another_set_of_proteins.py` and is used for comparing structures run on AlphaFold 2 with results from AlphaFold 3.
        - `final_interactors_filtered_by_model_rankings.py`: This script is used to integrate all the data, incorporate the model ranking scores, impose a cutoff (e.g., 0.6 ranking score) and generate a final figure where interactions below the cutoff that 'constrain' at least 20 amino acid residues in the bait are in grey. 
    - `Support_file_examples/`: Directory containing example files needed for running some of the scripts above.
        - `2024_5_4_A_list`: This is a multi-fasta file that provides all the 'bait' proteins followed by their sequence in successive lines. For example, '>A1\nMETADSCF...\n>A2\nMYPEDWQN...'. 
        - `2024_5_4_B_list`: This is a multi-fasta file that provides all the 'prey' proteins followed by their sequence in successive lines. For example, '>B1\nMETADSCF...\n>B2\nMYPEDWQN...'. 
        - `2024_5_4_B_list_finished`: This is an intermediate multi-fasta file that provides all the 'bait' proteins followed by their sequence in successive lines only for those that have already been run on the HPCC. This is useful for intermediate checks before all the models of all proteins have been computed, which can take many days depending on available compute time and queues on the HPCC.
        - `2024_5_4_A_list_sizes`: This file provides a list of 'bait' proteins and their sizes. Each line lists 'A,length'.
        - `2024_5_4_B_list_sizes`: This file provides a list of 'prey' proteins and their sizes. Each line lists 'B,length'.
        - `2024_5_9_alphafold2_summary_stats`: This file provides the summary statistics of the models collected together as a batch with a particular start date using  `alphafold2_dimer_batch_computed_on_zaratan.py` or using `alphafold3_dimer_batch_computed_on_google.py`. Each line lists 'A,B,maximum PAE,distance,interaction area'
        - `batch_dimer_predictions_alphafold2_zaratan.csv`: This file provides details of the runs that were performed. Each line lists 'A,B,input_date,output_date,pair_name'.
        - `batch_multimer`: This file lists the pairs of proteins (A and B) being tested for interactions. Each line lists a pair_name as 'AB'.
    - `Streamlined_Code/`
        - `AlphaFold_related`: AlphaFold runs are performed essentially as above using these two scripts:
          - `alphafold_multimer_with_clean_up.sh`: This is a script for running a single dimer by calling `run_alphafold.py`.
          - `alphafold_multimer_batch_submission_with_clean_up.sh`: This is a script for making the files needed for submitting dimers with a fixed time request.
        - `ChimeraX_related`: Directory containing scripts needed for processing results using ChimeraX and for analyzing the results. 
          - `predicted_dimer_chimerax.py`: The outputs of this script include all the residues that are under mutual constraint for the interacting partners, the surface area of interactions and movies highlighting the interactions. It is run from within ChimeraX.
          - `predicted_dimer_python.py`: This script computes on the interactions identified using ChimeraX by `predicted_dimer_chimerax.py`, and plots figures and prints tables in four parts. It is run from the terminal. 
        - For these streamlined scripts, the `2024_5_4_A_list` , `2024_5_4_B_list` , and `batch_multimer` from the `Support_file_examples` above are sufficient.
