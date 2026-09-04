import pandas as pd
import numpy as np
from itertools import combinations
import re
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

data_matrix = pd.read_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/fishers_results/wtdmso_wtcsc.tsv", index_col=0, sep="\t")

sites = data_matrix.index.str.replace("-mod|-unmod", "", regex=True).unique()
site_pattern = re.compile(r"^[^:]+:\d+$")
n_before = len(sites)
sites = [s for s in sites if site_pattern.match(s)]
print(f"Kept {len(sites)} of {n_before} sites matching chrom:genomepos format")

conds = data_matrix.columns.tolist()

def process_site(site):
    mod_row = f"{site}-mod"
    unmod_row = f"{site}-unmod"
    
    if mod_row not in data_matrix.index or unmod_row not in data_matrix.index:
        return []
    
    mod_counts = data_matrix.loc[mod_row].values
    unmod_counts = data_matrix.loc[unmod_row].values
    
    if np.all(mod_counts == 0) and np.all(unmod_counts == 0):
        return []
    
    results = []
    for cond1, cond2 in combinations(conds, 2):
        table = np.array([
            [mod_counts[conds.index(cond1)], mod_counts[conds.index(cond2)]],
            [unmod_counts[conds.index(cond1)], unmod_counts[conds.index(cond2)]]
        ])
        
        if table.sum() == 0:
            continue
        
        odds_ratio, pvalue = fisher_exact(table)
        
        results.append({
            "site": site,
            "cond1": cond1,
            "cond2": cond2,
            "mod_cond1": table[0,0],
            "mod_cond2": table[0,1],
            "unmod_cond1": table[1,0],
            "unmod_cond2": table[1,1],
            "odds_ratio": odds_ratio,
            "pvalue": pvalue
        })
    
    return results

all_results = []
for site in sites:
    all_results.extend(process_site(site))

results_df = pd.DataFrame(all_results)
sig_results = results_df[results_df["pvalue"] <= 0.05]

#I wan to also find the delta mod ratio

#maybe
#results_df["p_adj"] = multipletests(results_df["pvalue"], method="fdr_bh")[1]

sig_results.to_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/fishers_results/wtdmso_wtcsc_pairwise_fisher_results.tsv", sep="\t", index=False)
results_df.to_csv("/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/fishers_results/wtdmso_wtcsc_pairwise_fisher_unfiltered.tsv", sep="\t", index = False)

print("Done! Results saved")

