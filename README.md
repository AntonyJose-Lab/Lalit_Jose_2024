# Selecting genes for analysis using historically contingent progress: from RNA changes to protein-protein interactions

Farhaan Lalit<sup>1</sup> and Antony M. Jose<sup>1,*</sup>

<sup>1</sup>Department of Cell Biology and Molecular Genetics, University of Maryland, College Park, USA.

<sup>*</sup>Corresponding author email: amjose@umd.edu

This paper analyzes data from the field of RNA silencing in <i>C. elegans</i> to identify understudied regulated genes (URGEs) and uses AlphaFold to analyze the proteins encoded by them to reveal predicted influencers of RNA-regulated expression (PIREs). 

<b>URGEs</b>: Data tables were reformatted manually and/or using custom scripts, and filtered to generate lists of significant genes. The top ‘g’ genes that occur in the greatest numbers of tables were culled as the most frequently identified genes. A measure for the extent of regulation of each gene (<i>r<sub>g</sub></i>) was used to aid prioritization for detailed study. Co-occurrence patterns of genes in different tables were captured using the Jaccard distance (<i>d<sub>J</sub></i>) or as a symmetric measure of normalized mutual information, defined as Historical Mutual Information (HMI). The <i>d<sub>J</sub></i> values were used to generate a dendrogram using the average linkage method. HMI was used to group genes into clusters according to the Girvan-Newman algorithm and different sets of genes were highlighted. 

<b>PIREs</b>: Proteins encoded by the top 25 URGEs (16 proteins) were tested for predicted interactions with 25 key regulators of RNA silencing using a mix of AlphaFold 2 and/or AlphaFold 3. The resulting models were analyzed and filtered using various criteria to define PIRE proteins.    

The code used, key datasets, and supporting files are available here. Please change paths as necessary.

## Installation

The `requirements.txt` file lists the packages and dependencies used for the 'URGEs' portion of this project. To install all external packages, run `pip install -r requirements.txt` in your environment.

## Project Structure

- `URGEs/`
    - `app/`: Directory containing code and required files for the graphical user interface (GUI) for exploring regulated yet understudied genes. The app can be run locally by following these steps:
        1. Download the `app/` directory and its contents to your local machine.
        2. Set the `app/` folder as the working directory.
        3. In the terminal, run `python app.py`.
        4. Open a browser window and navigate to the URL given in the terminal after completing step 3.
        - `app.py`: Source code for the GUI.
        - `communities.csv `: A file containing clustered communities after running the app.
        - `dist_matrix.csv`: Distance matrix used to generate an adjacency matrix for creating the graph.
    - `code_and_resources/`: Directory containing all the code used for analyses and visualizations.
        - `data/`: Directory containing subdirectories to be used with the code relevant to standardizing files. These subdirectories are empty. Users should put their cleaned CSV files into the `cleaned/` folder to start with. After running `sanitization.py`, files will automatically get added to the `sanitized/` folder. Similarly, after running `filter_pvals.py`, files will automatically get added to the `filtered/` folder. 
            - `cleaned/`: Directory containing all the datasets in CSV format after cleaning.
            - `filtered/`: Directory containing all the datasets after being filtered for p-value < 0.05, where applicable.
            - `sanitized/`: Directory containing all the data files after gene name sanitization.
        - `community_graph.pickle`: File containing graph data for generation of Figure 2.
        - `dataset_wormexp_overlap.ipynb`: Code used to determine how many datasets chosen for this study overlap with the 461 included in WormExp 2.0 as on 27 Jan 2023.
        - `df_na_full.csv`: A CSV file that shows gene presence and absence in datasets.
        - `Extended_Data_Table.xlsx`: An Excel file containing names of all studies and their corresponding URLs used in this study.
        - `fig1D.ipynb`: Code used to generate Figure 1D. 
        - `fig1F_fig1G.ipynb`: Code used to generate Figure 1F and Figure 1G.
        - `fig5.ipynb`: Code used to generate Figure 5.
        - `filter_pvals.py`: Code used to filter datasets for p-value < 0.05.
        - `sanitization.py`: Code used to sanitize gene names.
        - `simplemine_results_100_11Feb2024.xlsx`: Simplemine results for top 100 genes, downloaded on February 11, 2024.
        - `simplemine_results_1000_31Jan2024.xlsx`: Simplemine results for top 1000 genes, downloaded on January 31, 2024.
        - `sklearn_nmi.ipynb`: Code used to cluster genes and create a distance matrix based on normalized mutual information between genes.
        - `TableOccupancy_full.csv`: CSV file listing how many datasets each gene occurs in.
        - `TableOccupancy.py`: Code used to create a TableOccupancy file listng how many datasets each gene occurs in.
        - `Wormbase_Gene_Sanitizer_Database_DownloadedOn_10-30-2023.txt`: Wormbase Gene Sanitizer Database downloaded on October 30, 2023.
        - `WormExp_2_0_info_20170727_Unique_Studies.xlsx`: An Excel file containing URLs for the 461 datasets included in WormExp 2.0 as on 27 Jan 2023.
        - `rg_simulation_box-whisker.ipynb`: A simulation looking at the relationship between S<sub>i</sub>, T<sub>i</sub>, and g.

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
    - `Support_File_examples/`: Directory containing example files needed for running some of the scripts above.
        - `2024_5_4_A_list`: This is a multi-fasta file that provides all the 'bait' proteins followed by their sequence in successive lines. For example, '>A1\nMETADSCF...\n>A2\nMYPEDWQN...'. 
        - `2024_5_4_B_list`: This is a multi-fasta file that provides all the 'prey' proteins followed by their sequence in successive lines. For example, '>B1\nMETADSCF...\n>B2\nMYPEDWQN...'. 
        - `2024_5_4_B_list_finished`: This is an intermediate multi-fasta file that provides all the 'bait' proteins followed by their sequence in successive lines only for those that have already been run on the HPCC. This is useful for intermediate checks before all the models of all proteins have been computed, which can take many days depending on available compute time and queues on the HPCC.
        - `2024_5_4_A_list_sizes`: This file provides a list of 'bait' proteins and their sizes. Each line lists 'A,length'.
        - `2024_5_4_B_list_sizes`: This file provides a list of 'prey' proteins and their sizes. Each line lists 'B,length'.
        - `2024_5_9_alphafold2_summary_stats`: This file provides the summary statistics of the models collected together as a batch with a particular start date using  `alphafold2_dimer_batch_computed_on_zaratan.py` or using `alphafold3_dimer_batch_computed_on_google.py`. Each line lists 'A,B,maximum PAE,distance,interaction area'
        - `batch_dimer_predictions_alphafold2_zaratan.csv`: This file provides details of the runs that were performed. Each line lists 'A,B,input_date,output_date,pair_name'.
        - `batch_multimer`: This file lists the pairs of proteins (A and B) being tested for interactions. Each line lists a pair_name as 'AB'.
