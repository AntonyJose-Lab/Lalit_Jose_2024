"""Code for sanitization of files"""

# Importing required modules
import pandas as pd
import os

# Importing the wormbase database
geneDB = pd.read_csv("Wormbase_Gene_Sanitizer_Database_DownloadedOn_10-30-2023.txt", comment="#", sep="\t")

# Defining a list of column names to search for
termList = ["Gene", "Genes", "genes", "gene", "Name", "name", "Gene Name", "Gene name", "gene name", "gene Name",
            "GeneName", "genename", "Genename", "Gene ID", "GeneID", "geneID", "gene_id", "gene ID", "symbol", "Symbol",
            "GENE", "geneName", "Gene locus", "Endo-siRNA locus", "Gene Locus", "piRNA", "Cosmid ID", "Cosmid", "cosmid", "Unnamed: 0",
            "wormbase gene id", "Sequence ID", "...1", "GENEID", "wormbase_gene", "Accession",
            "Sup. Table 1: 239 Examined epigenetic genes", "WormBaseID", "WBGene Name", "row_names", "Closest gene",
            "EAG", "Gene ID", "WBGene", "Wormbase ID:", "WORMBASE_GENE_ID (ce6)", "WORMBASE_GENE_ID ( ce10)",
            "Gene_Symbol", "CSR-1b enriched 22G-RNAs (Nguyen & Phillips, 2021)", "Gene.IDs", "WB.ID", "pseudogene",
            "CSR-1b enriched 22G-RNAs (Nguyen & Phillips, 2021)", "Transcript", "Gene Name 1", "This study and Kaminsky et al.",
            "This study and Drabikowski et al.", "CSR-1 targets exhibiting ectopic 22G-RNAs Related to Figure 4.", "Gene…1",
            "my_label", "P Granule Bound Transcripts", "Heritable Downregulated Genes and P Granule Bound Transcripts", 
            "Heritable Upregulated Genes and P Granule Bound Transcripts", "gene_name", "Protein", "Stable gene ID",
            "CSR-1 targets exhibiting ectopic 22G-RNAs Related to Figure 4.", "spermatogenesis genes"]

errorList = []  # Creating an empty list to append the names of all problematic files to
nonExistentColList = []    # Creating an empty list to consider the names of files with differing column names

fileList = os.listdir("data/cleaned") # Importing all files from the current directory

os.chdir("data/cleaned")  # Changing the current directory

fileCount = 0
l = len(fileList)

# Looping through each file
for file in fileList:
    fileCount += 1
    try:
        commentList = []    # Creating an empty list to contain all the comments
        os.chdir("data/cleaned")
        fhand = open(file)  # Creating a file handle and importing the csv as a text file
        
        # Checking if the line is a comment by checking if the first character is "#". If it is, it is appended to the comment list for that file
        for line in fhand:
            if line.startswith("#"):
                commentList.append(line)
        df0 = pd.read_csv(file, comment="#")    # Reading the csv as a dataframe
        df0 = df0.dropna(how='all')

        geneCount = 0
        inputList = []  # Creating an empty list to append all entries to the "Input" column
        suggestList = []    # Creating an empty list to append all entries to the "Suggested Match" column
        historyList = []    # Creating an empty list to append all entries to the "History" column

        for term in termList:   # Iterating through each element in the termList
            if term in df0.columns:     # Checking if the column exists in the dataframe
                for gene in (df0[term]):    # If the column exists, each gene in the column is checked.
                    if gene in (geneDB["Input"].values):    # If the gene exists in the "Input" column of the Wormbase Database, the corresponding suggested match and history are stored.
                        geneCount += 1
                        rownum = geneDB[geneDB["Input"] == gene].index[0]
                        inputValue = geneDB.at[rownum, "Input"]
                        suggestValue = geneDB.at[rownum, "Suggested Match"]
                        historyValue = geneDB.at[rownum, "History"]

                    else:   # If the gene does not exist in the "Input" column of the Wormbase Databse, the values passed is a blank string.
                        inputValue = ""
                        suggestValue = ""
                        historyValue = ""

                    # Appending the values from all three columns to their respective lists.
                    inputList.append(inputValue)
                    suggestList.append(suggestValue)
                    historyList.append(historyValue)

                    # Printing the following values just to know how many files have been converted and the progress of the program
                    print(file)
                    print(f"Num genes: {geneCount} / {len(df0)}")
                    print(f"Num files: {fileCount} / {l}")
                    print()
                    
                break
            else:
                nonExistentColList.append(file)
            
        # Converting all the lists into dataframes
        dfInput = pd.DataFrame(inputList)
        dfSuggest = pd.DataFrame(suggestList)
        dfHistory = pd.DataFrame(historyList)

        # Checking if any of the inputs were Empty Strings. If so, the file gets appended to the error list to be checked later.
        if "" in inputList:
            errorList.append(file)

        # Concatenating all the dataframes
        newdf = pd.concat([dfInput, dfSuggest], axis=1)
        newdf = pd.concat([newdf, dfHistory], axis=1)
        newdf = pd.concat([df0, newdf], axis=1)
        for term in termList:
            # Renaming the columns. By default, new columns get set to "0".
            if term in df0.columns:
                column_names = newdf.columns.values
                column_names[(len(column_names)) - 1] = "History"
                column_names[(len(column_names)) - 2] = "Suggested Match"
                column_names[(len(column_names)) - 3] = "Input"
                newdf.columns = column_names

        # Creating a new file name
        sanitized = "_sanitized.csv"
        fileName = file
        fileName1 = fileName.removesuffix(".csv")
        newFileName = fileName1 + sanitized

        # Changing the working directory and writing the comments and dataframe to the file.
        os.chdir("data/sanitized")
        f = open(newFileName, 'a')
        for comment in commentList:
            f.write(comment)
        newdf.to_csv(f, index=False)
        f.close()
    
    except Exception as e:  # Capture any exception that arises
        print(f"Error processing file {file}. Error message: {e}")
        errorList.append(f"{file} - Error: {e}")    # Add the file and error message to the error list

# Converting the error list to a dataframe and saving it.
errordf = pd.DataFrame(errorList, columns=["col1"])
coldf = pd.DataFrame(nonExistentColList, columns=["col1"])
errordf.to_csv("ERRORS.csv")
coldf.to_csv("ColNameIssues.csv")
