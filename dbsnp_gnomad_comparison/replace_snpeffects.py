import pandas as pd
import numpy as np

df = pd.read_csv('SnpEffects.tsv', sep='\t', header=None)

df.loc[df[1].str.contains('HIGH'), 1] = 'HIGH'
df.loc[df[1].str.contains('MODERATE'), 1] = 'MODERATE'
df.loc[df[1].str.contains('LOW'), 1] = 'LOW'
df['KNOWN_NOVEL'] = np.where(df[0]=='.', 'novel', 'known')
df.nunique()
df.to_csv('SnpEffects_correct.tsv', sep='\t', index=None, header=None)
