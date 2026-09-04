import pandas as pd

df = pd.read_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/data.indiv_genomic_samplenames.csv")

sample_map = {
    "wt1csc_4": "wt1csc",
    "w1c2-3": "wt1csc",
    "w1c_2": "wt1csc",


    "w2c2": "wt2csc",

    "wt1dmso_3": "wt1dmso",
    "w1d2": "wt1dmso",

    "w2d2": "wt2dmso",

    "m1_csc": "m1csc",
    "m1c2":"m1csc",
    
    "m2c2": "m2csc",

    "m1dmso_7": "m1dmso",
    "m1d2": "m1dmso",

    "m2d2": "m2dmso",
}

# Rename the sample names
df["sample_name"] = df["sample_name"].replace(sample_map)
set_samplename = set(df["sample_name"])
print(set_samplename)

# Save the renamed file
df.to_csv("data.indiv_genomic_samplenames_renamed.csv", index=False)
