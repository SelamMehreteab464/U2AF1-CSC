import pandas as pd
import sys

def rename_with_condition(colname):
    if colname.startswith("wt") and "dmso" in colname:
        return f"wtdmso"
    elif colname.startswith("wt") and "csc" in colname:
        return f"wtcsc"
    elif colname.startswith("mt") and "dmso" in colname:
        return f"mtdmso"
    elif colname.startswith("mt") and "csc" in colname:
        return f"mtcsc"
    else:
        return colname





def drop_zero_sites(sub, main_sample="wtdmso"):
    sub['base_site'] = sub.index.str.replace("-mod", "").str.replace("-unmod", "")
    
    mod_rows = sub.index.str.endswith("-mod")

    keep_main_sites = sub[mod_rows].loc[sub[mod_rows][main_sample] > 0, 'base_site'].unique()
    sub = sub[sub['base_site'].isin(keep_main_sites)]
            
    other_samples = [c for c in sub.columns if c not in ['base_site', main_sample]]
    def site_valid(df):
        return all(df[other_samples].iloc[0] > 0)  # only check mod row (first row of group)
                                    
    sub = sub.groupby('base_site').filter(site_valid)
                                        
    sub = sub.drop(columns=['base_site'])
    return sub

def main(input_file):
    df = pd.read_csv(input_file)
    wd_cols = [c for c in df.columns if c.startswith("wt") and "dmso" in c]
    wc_cols = [c for c in df.columns if c.startswith("wt") and "csc" in c]
    md_cols = [c for c in df.columns if c.startswith("mt") and "dmso" in c]
    mc_cols = [c for c in df.columns if c.startswith("mt") and "csc" in c]

    id_cols = ['chrom:genomepos']


    combos = {
            "wtdmso_wtcsc.tsv": wd_cols + wc_cols,
            "wtdmso_mtdmso.tsv": wd_cols + md_cols,
            "wtdmso_mtcsc.tsv": wd_cols + mc_cols,
            }


    for outname, cols in combos.items():
        sub = df[id_cols + cols].copy()
        samples =  [rename_with_condition(c) for c in cols]
        sub.columns = id_cols + samples
        sub.set_index(id_cols, inplace=True)
        #I want to add a line here that drops a row from all tables if wtdmso is 0 for mod and unmod and if wtdmso is not 0,0 then if one of the other samples are 0 mod and 0 unmod then drop that site from that df
        #sub = sub[~((sub['wtdmso'] == 0))]
        sub = drop_zero_sites(sub, main_sample="wtdmso")
        # Drop rows if any of the other samples are 0
        #other_samples = [c for c in sub.columns if c != 'wtdmso']
        #sub = sub[~(sub[other_samples] == 0).any(axis=1)]
        #sub.insert(0,'', range(len(sub)))
        sub.to_csv(outname, sep = "\t", index = True)

    print("Finished writing 3 tables: wd_wc.tsv, wd_md.tsv, wd_mc.tsv")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_combos.py <input.tsv>")
        sys.exit(1)
    main(sys.argv[1])

