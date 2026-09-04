import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import os
import csv
import logging
from typing import Dict, List, Any

# Constant values
READ_THRESHOLD = 10
MOD_PROB_THRESHOLD = 0.9

def parse_bed_file(bed_file: str) -> Dict[str, List[int]]:
    """Parses a BED file to extract cumulative exon edges for each transcript, excluding the first and last exons."""
    c, d = 0, 0
    isotoexonedges = {}
    mapping_to_gene = defaultdict(lambda: defaultdict(list))
    with open(bed_file, 'r') as bed:
        for line in bed:
            if line.startswith('#') or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            if len(fields) < 12:
                logging.warning(f"Skipping line due to insufficient columns: {line.strip()}")
                continue
            transcript_id = fields[3]
            chrom = fields[0]
            chrom_start = int(fields[1])
            chrom_end = int(fields[2])
            gene_name = fields[3]
            strand = fields[5]
            mapping_to_gene[chrom][chrom_start,chrom_end].append(gene_name)
            #print((type(mapping_to_gene))
            try:
                exon_sizes = list(map(int, fields[10].rstrip(',').split(',')))
                exon_starts = list(map(int, fields[11].rstrip(',').split(',')))
            except ValueError:
                logging.warning(f"Skipping line due to invalid exon sizes: {line.strip()}")
                continue

            d += 1
            if len(exon_sizes) > 1:
                c += 1
                cumulative_position = 0
                exon_edges = [0,] #exon edge list to start with 0.
                genomic_edges = []
    #            cumulative_genomic_edges = [chrom_start] 
     #           for i in range(len(exon_sizes)):
      #              cumulative_genomic_edges.append(cumulative_genomic_edges[-1] + exon_sizes[i])
                for i in range(len(exon_sizes)):
                    cumulative_position += exon_sizes[i]
                    exon_edges.append(cumulative_position)

                    #first i want all the genomic coordinates of exon starts

                    g_start = chrom_start + exon_starts[i]
                    
                    genomic_edges.append((g_start, exon_sizes[i]))

                    g_end = chrom_start + exon_starts[i] + exon_sizes[i]
                    genomic_edges.append((g_end, exon_sizes[i]))
    

                genomic_edges = genomic_edges[1:-1]
                
                #isotoexonedges.append({"chrom": chrom, "exon_edges" : exon_edges, "genomic_edges" : genomic_edges})
                #gene_name = mapping_to_gene
                isotoexonedges[gene_name] = genomic_edges
            

    #print(isotoexonedges)
    #print(mapping_to_gene)
    print(c/d, c, d)
    return isotoexonedges, mapping_to_gene

def process_sample_data(sample_file: str, mapping_to_gene: Dict[str, Dict[tuple, List[str]]], isotoexonedges: Dict[str, List[int]], read_threshold: int, mod_prob_threshold: float) -> List[List[float]]:
    """Processes m6A data to calculate modification probabilities relative to exon edges."""
    exonsizebins = [600,800,1000]
    all_bin_mod_values = [[[] for _ in range(2001)] for i in range(len(exonsizebins))]
    with open(sample_file, 'r') as m6a:
        csv_reader = csv.reader(m6a)
        next(csv_reader)  # Skip header row if present
        for row in csv_reader:
            chrom = row[0]
            genome_pos = int(row[1])
            read_count = int(row[3])
            mod_prob = float(row[4])
                #i really feel like the indentation thing is really throwing my game off here I will look at this more closely when I get home
            if read_count > read_threshold: #and gene_name in isotoexonedges
                mapping_to_mod = []
                if chrom in mapping_to_gene:
                    for (chrom_start,chrom_end), gene_list in mapping_to_gene[chrom].items():
                        if chrom_start <= genome_pos <= chrom_end:
                            mapping_to_mod.extend(gene_list)
                    #print(mapping_to_mod.keys())
                #print(mapping_to_gene.keys())
                    #if not mapping_to_mod:
                     #   continue

                for gene_name in mapping_to_mod:
                    mygenomicedges = isotoexonedges.get(gene_name, [])
                    distances = [genome_pos - edge for edge, size in mygenomicedges]
                    closest_distance = min(distances, key=abs)  # Find the closest distance
                    distindex = distances.index(closest_distance)
                    exonsize = mygenomicedges[distindex][1]
                    

                    ###new ccheck for whether pos is in first or last exon
                    #infirstlast = False
                    #if (distindex == 1 and closest_distance < 0) or (distindex == 0 and closest_distance >= 0): #in first exon
                     #   infirstlast = True
                    #if (distindex == len(distances) - 2 and closest_distance >= 0) or (distindex == len(distances) - 1):
                     #   infirstlast = True
                    
                    # if not infirstlast: #internal exons
                    #if infirstlast: #only first + last exons, remove internal
            #            if closest_distance >= 0:
             #           #exonedgeindexes = [distindex, distindex + 1]
              #              genomicindexes = [distindex, distindex + 1]
               #         elif closest_distance < 0:
                        #exonedgeindexes = [distindex -1, distindex]
                        #genomicindexes = [distindex -1, distindex]
                #            genomicindexes = [distindex - 1, distindex]
                        #print(mygenomicedges)
                        #exonsize = int(mygenomicedges[1])
                 
                    #print(exonsize)
                    #print(genomicindexes)
                    mybinindex = -1
                    for index in range(len(exonsizebins) - 1):
                        if exonsizebins[index] <= exonsize < exonsizebins[index + 1]:
                            mybinindex = index
                    if exonsize >= exonsizebins[-1]: mybinindex = len(exonsizebins) - 1
                    #if exonsize > 600: print(exonsize, mybinindex)

                    if -1000 <= closest_distance <= 1000 and mybinindex >= 0:
                        modposindex = closest_distance + 1000  # Map -1000 to 1000 -> index 0 to 2000
                        thresholded_value = 1 if mod_prob >= mod_prob_threshold else 0
                        all_bin_mod_values[mybinindex][modposindex].append(thresholded_value)

    #print(all_bin_mod_values[1][1000])
    return all_bin_mod_values#combined_relative_mod_values

def calculate_relative_means(all_bin_mod_values: List[List[List[float]]]) -> List[List[float]]:
    """Calculates mean modification probabilities for each relative position."""
    all_bin_rel_means = []
    for one_bin_values in all_bin_mod_values:
        relative_means = []
        for values in one_bin_values:
            if values:
                #relative_means.append(sum(values) / len(values))
                relative_means.append(len(values))
            else:
                relative_means.append(0)
        all_bin_rel_means.append(relative_means)
    return all_bin_rel_means

def bin_data(all_relative_means: List[List[float]], bin_size: int) -> List[List[float]]:
    """Bins the data into specified bin size for smoothing."""
    all_bin_data = []
    for relative_means in all_relative_means:
        binned_means = []
        for i in range(0, len(relative_means), bin_size):
            bin_values = relative_means[i:i+bin_size]
            binned_means.append(sum(bin_values) / len(bin_values))
        all_bin_data.append(binned_means)
    return all_bin_data

def plot_single_sample(dataset_name: str, sample_name: str,means: List[List[float]],bin_size: int = 50):
    fig, ax = plt.subplots(figsize=(8, 6))

    x_positions = np.arange(-1000, 1001, bin_size)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    bin_labels = ['600-800', '800-1000', '1000+']

    binned_means = bin_data(means, bin_size)

    for i, (vals, color) in enumerate(zip(binned_means, colors)):
        ax.plot(
            x_positions,
            vals,
            label=f"Exon size {bin_labels[i]}",
            linewidth=2,
            color=color,
            alpha=0.7
        )

    ax.axvline(0, color='red', linestyle='--')
    ax.set_xlabel("Distance to Closest Exon Edge", fontsize=14)
    #ax.set_ylabel("Number of modifiable positions", fontsize=14)
    ax.set_ylabel("fraction of pos modified", fontsize = 14)
    ax.set_title(f"{dataset_name} — {sample_name}", fontsize=16)
    ax.set_xlim(-500, 500)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()

    #output_file = f"{dataset_name}_{sample_name}_all_600upbins.png"
    output_file = f"{dataset_name}_{sample_name}_mean_all_600upbins.png"
    plt.savefig(output_file, dpi=600)
    plt.show()

    logging.info(f"Plot saved as {output_file}")


def analyze_and_plot_m6a_relative_positions(bed_file: str, sample_dirs: Dict[str, Dict[str, str]], read_threshold: int, mod_prob_threshold: float) -> None:

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Starting analysis and plotting")

    isotoexonedges, mapping_to_gene = parse_bed_file(bed_file)

    for dataset_name, sample_file in sample_dirs.items():
        logging.info(f"Processing {dataset_name}")

        mod_values = process_sample_data(
            sample_file,
            mapping_to_gene,
            isotoexonedges,
            read_threshold,
            mod_prob_threshold
        )

        means = calculate_relative_means(mod_values)

        plot_single_sample(
            dataset_name,
            dataset_name,
            means
        )


#bed_file = "/private/groups/brookslab/smehrete/RNA_modification/FLAIR/092525_cdna_FLAIR/100625_m6anet_combined_transcriptome.bed"

bed_file = "/private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.bed"


#sample_dirs = {
 #       "Combined_genomic": "/scratch/smehrete/10272025.m6anet.results/data.site_withsampleids.proba.csv"
  #      }

sample_dirs = {"07232026_raw_qc": "/scratch/smehrete/07212026_m6anet_results/data.site_withsampleids.proba.csv"}


analyze_and_plot_m6a_relative_positions(
    bed_file,
    sample_dirs,
    READ_THRESHOLD,
    MOD_PROB_THRESHOLD

)




