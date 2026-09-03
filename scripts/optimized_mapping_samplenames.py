import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
from typing import List, Tuple
import os


def map_readname(eventalign: str, indiv_file: str, chunksize: int = 5_000_000):
    df_indiv = pd.read_csv(indiv_file, sep='\t')
    df_indiv['read_index'] = df_indiv['read_index'].astype('int32')
    needed_keys = set(zip(df_indiv['transcript_id'],df_indiv['read_index']))
    print(f"Loaded indiv file: {len(df_indiv):,} rows, {len(needed_keys):,} unique keys")

    reader = pd.read_csv(
        eventalign, sep='\t', chunksize=chunksize,
        usecols=['contig', 'read_index', 'read_name'],
        dtype={'contig': 'category','read_index': 'int32', 'read_name': 'category'},
    )

    lookup = {}
    total_scanned = 0

    for i, chunk in enumerate(reader):
        total_scanned += len(chunk)
        chunk = chunk.rename(columns={'contig': 'transcript_id'})

        for tid, ridx, rname in zip(chunk['transcript_id'],chunk['read_index'], chunk['read_name']):
            key = (tid,ridx)
            if key in needed_keys and key not in lookup:
                lookup[key] = rname

        print(f"  chunk {i}: scanned {total_scanned:,} rows so far, "
              f"lookup has {len(lookup):,} entries ({len(lookup)/len(needed_keys):.1%} of needed keys)", end='\r')

    print(f"\nFinished streaming eventalign. Matched {len(lookup):,} / {len(needed_keys):,} keys.")
    return df_indiv, lookup


def extract_condition_from_fastq(fastq: str) -> dict:
    # write script to parse through fastq headers and grab the read name and sample id
    read_to_condition = {}
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
                read_name = sections[0]
                line_sections = next((t for t in sections if t.startswith("sampleid=")), None)
                if line_sections:
                    sample_name = line_sections.split("=")[1]
                    read_to_condition[read_name] = sample_name
    print(f"Extracted {len(read_to_condition):,} read-to-sample mappings")
    return read_to_condition


def compute_site_level_probs(df_indiv: pd.DataFrame, lookup: dict, read_to_condition: dict,
                              output_file: str, chunksize: int = 1_000_000):
    if os.path.exists(output_file):
        os.remove(output_file)

    n = len(df_indiv)
    wrote_header = False
    for start in range(0, n, chunksize):
        chunk = df_indiv.iloc[start:start + chunksize].copy()

        keys = list(zip(chunk['transcript_id'], chunk['read_index']))
        chunk['read_name'] = [lookup.get(k, np.nan) for k in keys]
        chunk['sample_name'] = chunk['read_name'].map(read_to_condition)

        chunk.to_csv(output_file, index=False, mode='a', header=not wrote_header)
        wrote_header = True
        print(f"  wrote rows {start:,}-{min(start + chunksize, n):,} of {n:,}", end='\r')

    print(f"\nSaved site-level probabilities to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Compute site proba from read proba.")
    parser.add_argument("--indiv_file", "-i", required=True, help="Input CSV file with per-read probabilities.")
    parser.add_argument("--eventalign_mapping", "-e", required=True, help="Input the eventalign.txt file to help map read names")
    parser.add_argument("--fastq", "-f", required=True, help="input fastq file")
    parser.add_argument("--output_file", "-o", required=True, help="Output CSV file for site-level probabilities.")
    parser.add_argument("--chunksize", type=int, default=5_000_000, help="Rows per chunk when streaming eventalign file")
    parser.add_argument("--write_chunksize", type=int, default=1_000_000, help="Rows per chunk when writing output")
    args = parser.parse_args()

    df_indiv, lookup = map_readname(args.eventalign_mapping, args.indiv_file, chunksize=args.chunksize)
    read_to_condition = extract_condition_from_fastq(args.fastq)

    compute_site_level_probs(
        df_indiv=df_indiv,
        lookup=lookup,
        read_to_condition=read_to_condition,
        output_file=args.output_file,
        chunksize=args.write_chunksize,
    )


if __name__ == "__main__":
    main()
