import pandas as pd
import sys

"""
Used to read in calculate relative abundances from qiime2. 
Suzette Palmer
06-15-2023

#To use: python relative_abundance.py <path to level-5.csv> <output file path/name with archea> <output file path/name no archea>
#Example path: "/home2/s180020/Desktop/Gut_BSI_16S_Data/OutputFile_noclinda7/level-5.csv"

Example command: python relative_abundance.py /home2/s180020/Desktop/Gut_BSI_16S_Data/OutputFile_noclinda7/level-5.csv /home2/s180020/Desktop/outfile_archea.csv /home2/s180020/Desktop/outfile_noarchea.csv
read in level-5.csv file from https://view.qiime2.org/

"""

def calrel_abundance(df):
    #calculate relative abundance and sort by phylum
    df2 = pd.DataFrame(df/df.sum(axis=0))
    val = pd.DataFrame(df2.index)
    val = val[0].str.split(';',expand=True)
    val.index = df2.index
    result = pd.concat([val, df2], axis=1, join='inner')
    sorted_df = result.sort_values(by=[1], ascending=True)
    return(sorted_df)


#this should be taken from the barplot.qzv output
#Path = sys.argv[0]
df = pd.DataFrame(pd.read_csv(str(sys.argv[1]), header =0, index_col='index'))
#transpose df and remove group name
df_T = df.T
df_T = df_T[:-1]
#remove archaea before calculations 
df_noa = df_T[df_T.index.str.startswith('d__Bacteria')]

df_archea = calrel_abundance(df_T)
df_noarchea = calrel_abundance(df_noa)

#export to csv
df_archea.to_csv(str(sys.argv[2]))
df_noarchea.to_csv(str(sys.argv[3]))
print("Success! Check your output file.")






