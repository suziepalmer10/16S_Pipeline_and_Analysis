import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])
new_df = df['sample-id'].str.split('_', expand=True)
df['sample-id1'] = new_df[0]
df['sample-id2'] = new_df[1]
df['sample-id3'] = new_df[2]
df_for = df[df["sample-id2"]=='R1']
df_rev = df[df["sample-id2"]=='R2']
df_merge = pd.merge(df_for, df_rev, on="sample-id1")
df_merge = pd.DataFrame(df_merge)
df_sub = df_merge[["sample-id1", " filepath_x", " filepath_y"]]
df_sub1 = df_sub.rename({"sample-id1": "sample-id",  " filepath_x":"forward-absolute-filepath", " filepath_y":"reverse-absolute-filepath"}, axis=1)
df_sub1.to_csv(sys.argv[2], index=False)  
