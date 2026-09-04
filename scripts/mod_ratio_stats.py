import pandas as pd
import numpy as np

# Load the data
df = pd.read_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/data.site_withsampleids.proba.csv")

df = df[(df["n_reads"] >= 10)]
# Compute number of modified and unmodified reads
df["reads_mod"] = (df["n_reads"] * df["mod_ratio"]).round().astype(int)
df["reads_unmod"] = df["n_reads"] - df["reads_mod"]

df["chrom:genomepos"] = df["chrom"].astype(str) + ":" + df["genomepos"].astype(str)


# Pivot for modified reads
df_mod = df.pivot(index="chrom:genomepos", columns="sample_name", values="reads_mod").fillna(0).astype(int)
df_unmod = df.pivot(index="chrom:genomepos", columns="sample_name", values="reads_unmod").fillna(0).astype(int)

# Add suffix to row names
df_mod.index = [f"{idx}-mod" for idx in df_mod.index]
df_unmod.index = [f"{idx}-unmod" for idx in df_unmod.index]

# Concatenate modified and unmodified rows one after the other
df_final = pd.concat([df_mod, df_unmod], axis=0)


#base_site = (df_final.index.str.replace("-mod", "", regex=False).str.replace("-unmod", "", regex=False))

#mod_order = df_final.index.str.endswith("-mod").map({True: 0, False: 1})

#df_final = df_final.iloc[np.lexsort((mod_order, base_site))]
# Ensure modified row comes immediately before unmodified row for each site
#df_final = df_final.sort_index(key=lambda x: [s.replace("-mod","").replace("-unmod","") for s in x])
df_final = df_final.sort_index(key=lambda x: list(zip([s.replace("-mod","").replace("-unmod","") for s in x],[0 if s.endswith("-mod") else 1 for s in x])))


#add up the mod and unmod values for the clones
#m1csc,m1dmso,m2csc,m2dmso,w1csc,w1dmso,w2csc,wt2dmso
df_final['wtdmso'] = df_final['wt1dmso'] + df_final['wt2dmso']
df_final['wtcsc'] = df_final['wt1csc'] + df_final['wt2csc']
df_final['mtdmso'] = df_final['m1dmso'] + df_final['m2csc']
df_final['mtcsc'] = df_final['m1csc'] + df_final['m2csc']


df_final = df_final.drop(columns=['wt2dmso'])
df_final = df_final.drop(columns=['wt1dmso'])
df_final = df_final.drop(columns=['m2dmso'])
df_final = df_final.drop(columns=['m1dmso'])
df_final = df_final.drop(columns=['m2csc'])
df_final = df_final.drop(columns=['m1csc'])
df_final = df_final.drop(columns=['wt2csc'])
df_final = df_final.drop(columns=['wt1csc'])


# Save to CSV
df_final.to_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/m6anet_output_clones_combined.csv", index=True, index_label="chrom:genomepos")
print("Saved reshaped data to output.csv") 
