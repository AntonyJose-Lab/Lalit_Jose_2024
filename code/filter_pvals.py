"""Code for filtering of records for a p-value of < 0.05."""

import pandas as pd
import os

# Function to safely convert a column to float.
def safe_float_conversion(series):
    try:
        return series.astype(float)
    except:
        return series

# Setting the working directory and reading in the list of files.
os.chdir("data/sanitized")
fileList = os.listdir("data/sanitized")
fileList = [f for f in fileList if not f.startswith(".")]  # To avoid any hidden files

match1 = ["adj", "Qvalue", "q_value", "q-value", "FDR", "padj"]
match2 = ["Pvalue", "p value", "p-value"]
col_prob = list()

for file in fileList:
    df = pd.read_csv(file, comment="#")

    for col in df.columns:
        colname = col.lower()

        # Convert potential p-value column to float
        df[col] = safe_float_conversion(df[col])

        if df[col].dtype == 'float64':
            # First, checking for any match in match1
            if any(x in colname for x in match1):
                df_fil = df[df[col] <= 0.05]
                filename = file.removesuffix(".csv")
                filename = filename + "_filtered.csv"
                os.chdir("data/filtered")
                df_fil.to_csv(filename, index=False)
                os.chdir("data/filtered")
            # If no match in match1, checking for any match in match2
            elif any(x in colname for x in match2):
                df_fil = df[df[col] <= 0.05]
                filename = file.removesuffix(".csv")
                filename = filename + "_filtered.csv"
                os.chdir("data/filtered")
                df_fil.to_csv(filename, index=False)
                os.chdir("data/sanitized")
            else:
                col_prob.append(file)
                print(file)

df1 = pd.DataFrame([col_prob])
df1 = df1.transpose()
os.chdir("data/filtered")
