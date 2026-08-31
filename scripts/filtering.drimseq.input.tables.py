import pandas as pd
import sys

def rename_with_condition(colname):
    if "_S" in colname:
        colname = colname.replace("_S", "-S")
    base = colname.split("_")[0]
    if colname.startswith("W") and "D" in colname:
        return f"{base}_wtdmso_b1"
    elif colname.startswith("W") and "C" in colname:
        return f"{base}_wtcsc_b1"
    elif colname.startswith("M") and "D" in colname:
        return f"{base}_mtdmso_b1"
    elif colname.startswith ("M") and "C" in colname:
        return f"{base}_mtcsc_b1"
    else:
        return colname

def main(input_file):
    df = pd.read_csv(input_file, sep = "\t")

    wd_cols = [c for c in df.columns if c.startswith("W") and "D" in c]
    wc_cols = [c for c in df.columns if c.startswith("W") and "C" in c]
    md_cols = [c for c in df.columns if c.startswith("M") and "D" in c]
    mc_cols = [c for c in df.columns if c.startswith("M") and "C" in c]

    id_cols = ['gene_id', 'feature_id']


    combos = {
            "wtdmso_wtcsc.tsv": wd_cols + wc_cols,
            "wtdmso_mtdmso.tsv": wd_cols + md_cols,
            "wtdmso_mtcsc.tsv": wd_cols + mc_cols,
            "wtdmso_all.tsv": wd_cols + wc_cols + md_cols + mc_cols
            }


    for outname, cols in combos.items():
        sub = df[id_cols + cols].copy()
        sub.columns = id_cols + [rename_with_condition(c) for c in cols]
        sub.insert(0,'', range(len(sub)))
        sub.to_csv(outname, sep = "\t", index = False)

    print("Finished writing 4 tables: wd_wc.tsv, wd_md.tsv, wd_mc.tsv, wd_all.tsv")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_combos.py <input.tsv>")
        sys.exit(1)
    main(sys.argv[1])

