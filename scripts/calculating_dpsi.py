import argparse
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-i', type=str, action='store', help= 'drim output file', required = True)
parser.add_argument('--outFile', '-o', type=str, action='store', help = 'table with dpsi vals', required = True)
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile


def parser_drim(inFile):
    genes = []
    with open(inFile, 'r') as f:
        header = f.readline().strip().split('\t')

        psi_cols = header[2:9]
        for line in f:
            cols = line.strip().split('\t')
            feature_id = cols[0].split(';')

            inclusion_exclusion = feature_id[0].split('_')[0] if len(feature_id[0].split('_')) > 0 else 'NA'
            contains_novel_or_only_known_junctions = feature_id[0].split('_')[1] if len(feature_id[0].split('_')) > 1 else 'NA'
            
            if inclusion_exclusion != 'inclusion':
                continue


            while len(feature_id) < 5:
                feature_id.append('NA')
            
            genes.append({'inclusion_exclusion': inclusion_exclusion, 'contains_novel_or_only_known_junctions':contains_novel_or_only_known_junctions,'as_event_type':feature_id[1],'gene_name': feature_id[2],'gene_coordinates':feature_id[4], 'strand':feature_id[3], 'psi_stat_values': cols[2:10], 'lr': cols[10], 'padj': cols[11]})

        return genes,psi_cols

def calculate_dpsi(genes):
    dpsi_val = []
    skipped_count = 0 
    for gene in genes:
        psi_vals = gene['psi_stat_values']
        try:
            wt_psi = float(psi_vals[0]) * 100
            mut_psi = float(psi_vals[4]) * 100
            dpsi = mut_psi - wt_psi
            gene['dpsi'] = dpsi
            dpsi_val.append(gene)
        except (ValueError, IndexError) as e:
            print(f"Error processing PSI values {psi_vals}: {e}")
            skipped_count += 1
            continue
                                                
    print(f"Calculated dPSI for {len(dpsi_val)} genes")
    return dpsi_val


def write_output(dpsi_genes, outFile, psi_cols):
    # Define the columns
    columns = [
        'inclusion_exclusion',
        'contains_novel_or_only_known_junctions',
        'as_event_type',
        'gene_name',
        'gene_coordinates',
        'strand'
    ] + psi_cols + ['lr', 'padj', 'dpsi']  

    with open(outFile, 'w') as f:
        f.write('\t'.join(columns) + '\n')
        for gene in dpsi_genes:
            row = [
                gene['inclusion_exclusion'],
                gene['contains_novel_or_only_known_junctions'],
                gene['as_event_type'],
                gene['gene_name'],
                gene['gene_coordinates'],
                gene['strand']
            ] + [str(x) for x in gene['psi_stat_values']] + [gene['lr'], gene['padj'], str(gene['dpsi'])]
            f.write('\t'.join(row) + '\n')

genes, psi_cols = parser_drim(inFile)
dpsi_genes = calculate_dpsi(genes)

if not dpsi_genes:
    print("No genes with valid PSI values! Exiting.")
else:
    print(f"Writing {len(dpsi_genes)} genes to output file.")
    write_output(dpsi_genes, outFile, psi_cols)

