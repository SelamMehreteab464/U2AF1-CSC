import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

def extract_condition_from_fastq(fastq: str) -> dict:
    read_to_sample = {}
    with open(fastq, 'r') as read_input:
        while True:
            header = read_input.readline().strip()
            if not header:
                break
            read_input.readline()  
            read_input.readline()  
            read_input.readline()  

            if header.startswith('@'):
                sections = header[1:].split()
                # Your working string isolation logic
                read_name = sections[0] if sections else ""
                line_sections = next((t for t in sections if t.startswith("sampleid=")), None)
                if line_sections:
                    sample_name = line_sections.split("=")[1]
                    read_to_sample[read_name] = sample_name
                    
    print(f"Extracted {len(read_to_sample):,} read-to-sample mappings")
    return read_to_sample

def map_sample_to_condition(sample_name: str) -> str:
    if pd.isna(sample_name):
        return np.nan
    
    # Standardize to lowercase for robust string matching
    sample_lower = str(sample_name).lower()
    
    if 'w' in sample_lower and '1' in sample_lower and 'd' in sample_lower:
        return 'wt1dmso'
    elif 'w' in sample_lower and '2' in sample_lower and 'd' in sample_lower:
        return 'wt2dmso'
    elif 'w' in sample_lower and '1' in sample_lower and 'c' in sample_lower:
        return 'wt1csc'
    elif 'w' in sample_lower and '2' in sample_lower and 'c' in sample_lower:
        return 'wt2csc'
    elif 'm' in sample_lower and '1' in sample_lower and 'd' in sample_lower:
        return 'mt1dmso'
    elif 'm' in sample_lower and '2' in sample_lower and 'd' in sample_lower:
        return 'mt2dmso'
    elif 'm' in sample_lower and '1' in sample_lower and 'c' in sample_lower:
        return 'mt1csc'
    elif 'm' in sample_lower and '2' in sample_lower and 'c' in sample_lower:
        return 'mt2csc'
    else:
        return 'unknown'

    
def process_tsv_and_box_plot(tsv_input: str, tsv_output: str, fastq_file: str):
    read_to_sample = extract_condition_from_fastq(fastq_file)
    
    print(f"Reading {tsv_input}...")
    df = pd.read_csv(tsv_input, sep='\t')
    
    df['sample_name'] = df['readname'].map(read_to_sample)
    
    df['condition'] = df['sample_name'].apply(map_sample_to_condition)
    
    df.to_csv(tsv_output, sep='\t', index=False)
    print(f"Saved updated dataset to {tsv_output}")
    
    plot_df = df[
        (df['qc_tag'] == 'PASS') & 
        (df['condition'].isin(['wt1dmso', 'wt2dmso','wt1csc', 'wt2csc', 'mt1dmso', 'mt2dmso', 'mt1csc', 'mt2csc']))
    ].copy()


    if plot_df.empty:
        print("Warning: No matching 'PASS' reads found. Skipping plot generation.")
        return

    plot_df['clone'] = plot_df['condition'].str.extract(r'([12])')[0]

    plot_df['condition_group'] = plot_df['condition'].str.replace(r'[12]', '', regex=True)

    clone_order = ['1', '2']

    #condition_order = ['wt1dmso','wt2dmso', 'wt1csc', 'wt1csc','wt2csc','mt1dmso', 'mt2dmso', 'mt1csc', 'mt2csc']
    condition_order = ['wtdmso', 'wtcsc', 'mtdmso', 'mtcsc']


    #custom_palette = {
     #   'wt1dmso': 'gray',
      #  'wt2dmso': 'gray',
       # 'wt1csc': 'red',
        #'wt2csc': 'red',
        #'mt1dmso': 'blue',
        #'mt2dmso': 'blue',
        #'mt1csc': 'purple',
        #'mt2csc': 'purple'
    #}
    custom_palette = {'wtdmso':'gray', 'wtcsc':'red', 'mtdmso':'blue', 'mtcsc': 'purple'}

    plt.figure(figsize=(10, 7))

    sns.boxplot(
        data=plot_df,
        x='clone',
        y='polya_length',
        order=clone_order,
        hue_order=condition_order,
        palette=custom_palette,
        legend=True,
        showfliers=False)


    #pairs_to_test = [
     #   ("wt1dmso", "wt1csc"),
      #  ("wt2dmso", "wt2csc"),
       # ("wt1dmso", "mt1dmso"),
        #("wt2dmso", "mt2dmso"),
        #("wt1dmso", "mt1csc"),
        #("wt2dmso", "mt2csc")
    #]



    
    print("\n" + "="*50)
    print("STATISTICAL MANN-WHITNEY U TEST RESULTS:")
    print("="*50)
    
    annotator = Annotator(max, pairs_to_test, data=plot_df, x='clone', y='polya_length', order=condition_order)
    
    annotator.configure(
        test='Mann-Whitney', 
        text_format='star', 
        loc='inside', 
        verbose=True
    )
    annotator.apply_and_annotate()
    print("="*50 + "\n")

    plt.title('PolyA Length Distribution Comparison (QC PASS)', fontsize=14, weight='bold', pad=15)
    plt.xlabel('clone', fontsize=12, labelpad=10)
    plt.ylabel('PolyA Length', fontsize=12, labelpad=10)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plot_output = "polya_boxplot_wt_control_only.png"
    plt.savefig(plot_output, dpi=300)
    plt.show()
    print(f"Boxplot significance plot saved successfully as {plot_output}")


    # =========================================================================
    print("Generating Multi-Condition Histogram Line Plot...")
    plt.figure(figsize=(10, 6))

    # Using element="poly" and fill=False turns the histogram blocks into single tracing lines
    sns.histplot(
        data=plot_df,
        x='polya_length',
        hue='condition',
        hue_order=condition_order,
        palette=custom_palette,
        element='poly',
        fill=False,
        stat='probability',
        common_norm=False,
        bins=50,
        linewidth=2.5
    )

    plt.title('PolyA Length Frequency Distribution by Condition', fontsize=14, weight='bold', pad=15)
    plt.xlabel('PolyA Length (nt)', fontsize=12, labelpad=10)
    plt.ylabel('Norm Read Count', fontsize=12, labelpad=10)
    plt.xlim(0, 300)  # Capping the window focus on standard biological tail boundaries
    plt.grid(axis='both', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    hist_output = "polya_histogram.png"
    plt.savefig(hist_output, dpi=300)
    plt.show()
    print(f"Histogram line plot saved successfully as {hist_output}")


if __name__ == "__main__":
    input_tsv = "polya_results.tsv"
    output_tsv = "polya_results_samplename.tsv"
    fastq_path = "/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample/all_samples.fastq"
    
    process_tsv_and_box_plot(input_tsv, output_tsv, fastq_path)

