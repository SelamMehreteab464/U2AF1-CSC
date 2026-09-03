import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
from typing import List, Tuple
from Bio import SeqIO

def _calculate_site_proba(task: Tuple[np.ndarray, int, int]) -> float:
	probs, num_iterations, n_samples = task
	probs = np.random.choice(probs, num_iterations * n_samples, replace=True).reshape(num_iterations, n_samples)
	site_prob = (1 - np.prod(1 - probs, axis=1)).mean()
	return site_prob


def calculate_site_proba(read_probs: List[np.ndarray], num_iterations: int, n_samples: int, n_processes: int) -> List[float]:
	tasks = [(p, num_iterations, n_samples) for p in read_probs]
	with Pool(n_processes) as pool:
		site_probs = pool.map(_calculate_site_proba, tasks)
	return site_probs


def compute_site_level_probs(df_indiv: pd.DataFrame, output_file: str, num_iterations: int, n_samples: int,
                             n_processes: int, threshold: float):

    #print(df["transcript_id"].unique()[:20])
    #print(df["transcript_id"].dtype)
    grouped = df_indiv.groupby(["transcript_id", "transcript_position"], sort=False)["probability_modified"].apply(np.array)

    site_probs = calculate_site_proba(grouped.tolist(), num_iterations, n_samples, n_processes)

    mod_ratios = [np.mean(x >= threshold) for x in grouped]


    result = pd.DataFrame({
        "transcript_id": [g[0] for g in grouped.index],
        "transcript_position": [g[1] for g in grouped.index],
        "n_reads": [len(x) for x in grouped],
        "site_proba": site_probs,
        "mod_ratio": mod_ratios})

    result.to_csv(output_file, index=False)
    print(f"Saved site-level probabilities to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Compute site proba from read proba.")
    parser.add_argument("--indiv_file","-i", required=True, help="Input CSV file with per-read probabilities.")
    parser.add_argument("--output_file","-o", required=True, help="Output CSV file for site-level probabilities.")
    parser.add_argument("--num_iterations", type=int, default=1000, help="Number of iterations for sampling.")
    parser.add_argument("--n_samples", type=int, default=20, help="Number of samples per iteration.")
    parser.add_argument("--n_processes", type=int, default=4, help="Number of parallel processes.")
    parser.add_argument("--threshold", type=float, default=0.033379376, help="Threshold for modified read ratio.")
    args = parser.parse_args()

	# Step 1: load per-read probabilities

    df_indiv = pd.read_csv(args.indiv_file, sep='\t')

    compute_site_level_probs(
        df_indiv=df_indiv,
        output_file=args.output_file,
        num_iterations=args.num_iterations,
        n_samples=args.n_samples,
        n_processes=args.n_processes,
        threshold=args.threshold
    )


if __name__ == "__main__":
    main()

