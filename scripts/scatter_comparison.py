import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser()

parser.add_argument("--inFile", "-i", required = True, help = "input table to be plotted")
parser.add_argument("--site_file", "-s", required = True, help = "m6anet output data.site.proba")
parser.add_argument("--output", "-o", required = True, help = "output")
args = parser.parse_args()
inFile = args.inFile
site_file = args.site_file
output = args.output

df = pd.read_csv(inFile, sep = ",")
df_site = pd.read_csv(site_file, sep= ",")


sns.set(style="whitegrid")
plt.figure(figsize=(8,5))
sns.scatterplot(
    x=df["mod_ratio"],   
    y=df_site["mod_ratio"],    
    s=15,
    alpha=0.6
)

plt.plot([0,1], [0,1], ls="--", color="black")
plt.xlabel("modratio_calculated")
plt.ylabel("modratio_m6anet")

plt.savefig(args.output, dpi=1500)


