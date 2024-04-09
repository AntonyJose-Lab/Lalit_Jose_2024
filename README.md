# Selecting genes for analysis using historically contingent progress

Farhaan Lalit<sup>1</sup> and Antony M. Jose<sup>1,*</sup>

<sup>1</sup>Department of Cell Biology and Molecular Genetics, University of Maryland, College Park, USA.

<sup>*</sup>Corresponding author email: amjose@umd.edu

This paper analyzes data from the field of RNA silencing in <i>C. elegans</i>. Data tables were reformatted manually and/or using custom scripts, and filtered to generate lists of significant genes. The top ‘g’ genes that occur in the greatest numbers of tables were culled as the most frequently identified genes. A measure for the extent of regulation of each gene (<i>r<sub>g</sub></i>) was used to aid prioritization for detailed study. Co-occurrence patterns of genes in different tables were captured using the Jaccard distance (<i>d<sub>J</sub></i>) or as a symmetric measure of normalized mutual information, defined as Historical Mutual Information (HMI). The <i>d<sub>J</sub></i> values were used to generate a dendrogram using the average linkage method. HMI was used to group genes into clusters according to the Girvan-Newman algorithm and different sets of genes were highlighted. The code used and the key datasets are available here.

## Installation

The `requirements.txt` file lists the packages and dependencies used for this project. To install all external packages, run `pip install -r requirements.txt` in your environment.

## Project Structure

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
    - `fig2.ipynb`: Code used to generate Figure 2.
    - `filter_pvals.py`: Code used to filter datasets for p-value < 0.05.
    - `sanitization.py`: Code used to sanitize gene names.
    - `simplemine_results_100_11Feb2024.xlsx`: Simplemine results for top 100 genes, downloaded on February 11, 2024.
    - `simplemine_results_1000_31Jan2024.xlsx`: Simplemine results for top 1000 genes, downloaded on January 31, 2024.
    - `sklearn_nmi.ipynb`: Code used to cluster genes and create a distance matrix based on normalized mutual information between genes.
    - `TableOccupancy_full.csv`: CSV file listing how many datasets each gene occurs in.
    - `TableOccupancy.py`: Code used to create a TableOccupancy file listng how many datasets each gene occurs in.
    - `Wormbase_Gene_Sanitizer_Database_DownloadedOn_10-30-2023.txt`: Wormbase Gene Sanitizer Database downloaded on October 30, 2023.
    - `WormExp_2_0_info_20170727_Unique_Studies.xlsx`: An Excel file containing URLs for the 461 datasets included in WormExp 2.0 as on 27 Jan 2023.
    - `rg_simulation_box-whisker.ipynb`: A simulation looking at the relationship between Si, Ti, and g.
