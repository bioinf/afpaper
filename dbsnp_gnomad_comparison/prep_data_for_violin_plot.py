import pandas as pd
import numpy as np

df = pd.read_csv('all_data_AF_vs_GNOMAD.tsv', sep='\t', header=None)
df['KNOWN_NOVEL'] = np.where(df[0]=='.', 'novel', 'known')
df = df.loc[~df[1].str.contains(','),]
df.to_csv('data_for_violin_plots.tsv', header=None, index=None, sep='\t')
