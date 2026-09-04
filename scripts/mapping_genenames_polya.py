import os
import re
import pandas as pd

def parse_gtf_gene_names(gtf_path: str) -> dict:
    """
    Parses a GTF file to build a mapping dictionary from 
    Ensembl Gene IDs (ENSG) to Gene Names (Symbols).
    """
    ensg_to_name = {}
    print(f"Parsing GTF file: {gtf_path}...")
    
    gene_id_regex = re.compile(r'gene_id "([^"]+)"')
    gene_name_regex = re.compile(r'gene_name "([^"]+)"')
    
    with open(gtf_path, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
                
            sections = line.strip().split('\t')
            if len(sections) >= 9 and sections[2] == 'gene':
                attributes = sections[8]
                
                id_match = gene_id_regex.search(attributes)
                name_match = gene_name_regex.search(attributes)
                
                if id_match and name_match:
                    ensg_id = id_match.group(1)
                    gene_name = name_match.group(1)
                    
                    clean_ensg = ensg_id.split('.')[0]
                    ensg_to_name[clean_ensg] = gene_name
                    
    print(f"Extracted {len(ensg_to_name):,} unique gene ID mappings from GTF.")
    return ensg_to_name

def extract_ensg_from_contig(contig_string: str) -> str:
    """
    Extracts the clean ENSG ID from the contig string.
    Example: '2-1_ENST00000327044.7_ENSG00000188976.11' -> 'ENSG00000188976'
    """
    if pd.isna(contig_string):
        return None
        
    # Search for the string segment starting with ENSG
    match = re.search(r'(ENSG\d+)', str(contig_string))
    if match:
        return match.group(1)
    return None

def add_gene_names_to_results(input_tsv: str, output_tsv: str, gtf_path: str):
    ensg_to_name = parse_gtf_gene_names(gtf_path)
    
    print(f"Loading data from {input_tsv}...")
    df = pd.read_csv(input_tsv, sep='\t')
    
    print("Extracting Ensembl Gene IDs from contig metadata...")
    df['ensg_id'] = df['contig'].apply(extract_ensg_from_contig)
    
    print("Mapping gene names...")
    df['gene_name'] = df['ensg_id'].map(ensg_to_name)
    
    df['gene_name'] = df['gene_name'].fillna(df['ensg_id'])
    
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Successfully processed file! Saved with gene symbols to: {output_tsv}")
    
    print("\nSample Preview:")
    print(df[['contig', 'ensg_id', 'gene_name']].dropna().head())

if __name__ == "__main__":
    input_file = "polya_results_samplename.tsv"
    output_file = "polya_results_with_genenames.tsv"
    gtf_file = "/private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf" 
    
    if os.path.exists(input_file):
        add_gene_names_to_results(input_file, output_file, gtf_file)
    else:
        print(f"Error: Could not locate input file '{input_file}' in your workspace directory.")

