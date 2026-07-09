import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

df_mut_only = pd.read_csv('U2AF1_specific_mutation_status.tsv', sep = '\t')
df_all_samples = pd.read_csv('U2AF1_mutation_vs_cigarette_use.txt', sep='\t')

print(df_all_samples.columns.tolist())

df_all_samples["Sample Id"] = df_all_samples["Sample Id"].astype(str).str.strip()
df_all_samples = df_all_samples.rename(columns={'Sample Id' : 'Sample ID'})
df_mut_only["Sample ID"] = df_mut_only["Sample ID"].astype(str).str.strip()

merged = df_all_samples.merge(
    df_mut_only[["Sample ID", "Protein Change"]],
    on="Sample ID",
    how="left"
)

def classify_u2af1_mut_status(protein_change):
    if pd.isna(protein_change) or str(protein_change).strip() == "":
        return "none"

    pc = str(protein_change).upper()

    if re.search(r"S34F", pc):
        return "s34f"
    return "other"

merged["U2AF1 mutation status"] = merged["Protein Change"].apply(classify_u2af1_mut_status)


merged["U2AF1 mutation status"] = pd.Categorical(
    merged["U2AF1 mutation status"],
    categories=["none", "s34f"],
    ordered=True
)

#def collapse_smoke_status(smoke_status):
#    if pd.isna(smoke_status):
#        return np.nan
#    smoke_status = smoke_status.strip()
    # Combine Current user + Former user + Former unknown time
#    if smoke_status in ["Current user", "Former user (quit >1 year)", 
 #                       "Former user (quit < 1 year)", "Former user (unknown time)"]:
 #       return "Former/Current user"
  #  elif smoke_status in ["Never used"]:
   #     return "none"
   # else:
   #     return np.nan

smoke_col = "Non-Small Cell Lung Cancer: Cigarette Use at Time of Diagnosis"
#merged["former/Current User"] = merged[smoke_col].apply(collapse_smoke_status)
#merged = merged.dropna(subset=["former/Current User"])

smoke_table = pd.crosstab(
            merged["U2AF1 mutation status"],
                merged[smoke_col]
                )


smoke_table = smoke_table.reindex(["none", "s34f"])

print(smoke_table)




smoke_table.to_csv("u2af1_smoking_fisher_table.tsv", sep="\t")


smoke_col = "Non-Small Cell Lung Cancer: Cigarette Use at Time of Diagnosis"

plt.rcParams.update({'font.size': 15}) 

ct = smoke_table

ax = ct.plot(kind="bar", stacked=True, figsize=(7,5))

ax.set_xlabel("U2AF1 mutation status")
ax.set_ylabel("Number of samples")
ax.set_title("Smoking status by U2AF1 mutation group")
ax.legend(title="Smoker status", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.savefig("counts_u2af1_mut_smoking.png", format = "png", bbox_inches='tight')

ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

ax = ct_pct.plot(kind="bar", stacked=True, figsize=(7,5))

ax.set_xlabel("U2AF1 mutation status")
ax.set_ylabel("Percent of samples (%)")
#ax.set_title("Smoking status by U2AF1 mutation group (percent)")
ax.legend(title="Smoker status", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.savefig("percentage_u2af1_mut_smoking.pdf", format = "pdf", dpi = 1200, bbox_inches='tight')





