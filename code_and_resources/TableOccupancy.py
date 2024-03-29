# Table Occupancy single occurrences - Filtered

# Importing required modules
import pandas as pd
import os

fileList = os.listdir("../data_full")    # Importing all files from the directory
fileList.remove(".DS_Store") if ".DS_Store" in fileList else None # Removing a system generated file

os.chdir("../data_full")     # Changing the working directory

genedict = {}   # Initializing an empty dictionary
count = 0   # Initializing a counter

geneList = []

# Iterating through each file
for file in fileList:
    count = count + 1
    print(count)
    print(file)
    df = pd.read_csv(file, comment="#")  # Reading each file as a dataframe
    geneList = []
    for gene in df["Suggested Match"]:  # Iterating through each gene in the Suggested Match column
        genestr = str(gene)     # Converting the contents to the value to a string
        if "," in genestr:  # Splitting values which have multiple genes present. Denoted by a "," between the two genes.
            print(genestr)
            gene1 = genestr[0:14]   # The first gene
            print(gene1)
            gene2 = genestr[16:30]      # The second gene
            print(gene2)
            genestr = gene1
            if genestr not in geneList:
                geneList.append(genestr)
                if genestr != "nan":
                    if genestr in genedict:
                        genedict[genestr] = genedict[
                                                genestr] + 1  # If the gene already exists in the dictionary, it's value is incremented
                    else:
                        genedict[
                            genestr] = 1  # If the gene does not already exist in the dictionary, it is added to the dictionary
            genestr = gene2     # The second gene in the cell
            if genestr not in geneList:
                geneList.append(genestr)
                if genestr != "nan":
                    if genestr in genedict:
                        genedict[genestr] = genedict[
                                                genestr] + 1  # If the gene already exists in the dictionary, it's value is incremented
                    else:
                        genedict[
                            genestr] = 1  # If the gene does not already exist in the dictionary, it is added to the dictionary


        if genestr not in geneList:
            geneList.append(genestr)
            if genestr != "nan":
                if genestr in genedict:
                    genedict[genestr] = genedict[genestr] + 1     # If the gene already exists in the dictionary, it's value is incremented
                else:
                    genedict[genestr] = 1  # If the gene does not already exist in the dictionary, it is added to the dictionary
 
for gene in genedict:
    
    # changing F07B7.1 to W09B7.1
    if gene == "WBGene00017185":
        genedict["WBGene00021106"] = genedict["WBGene00021106"] + genedict["WBGene00017185"]
    
    # changing F07B7.2 to W09B7.2
    if gene == "WBGene00017186":
        genedict["WBGene00021107"] = genedict["WBGene00021107"] + genedict["WBGene00017186"]
    
del genedict["WBGene00017185"] 
del genedict["WBGene00017186"]
                    
print(count)

df = pd.DataFrame([genedict])   # Converting the dictionary to a dataframe
df = df.transpose()     # Transposing the dataframe
df = df.dropna()    # Dropping null values from the dataframe

    
os.chdir("../code_and_resources")
df.to_csv("TableOccupancy_full.csv", index=True)     # Saving the Dataframe as a CSV.
