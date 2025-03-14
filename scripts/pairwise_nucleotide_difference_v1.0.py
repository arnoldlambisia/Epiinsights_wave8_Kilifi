#A script that calculates pairwise nucleotide differences from a multiple sequence alignment
#Created by Arnold Lambisia and ChatGPT
#Email: alambisia@kemri-wellcome.org
#This script can be used under a Creative Commons (CC) license.


import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from itertools import combinations

def pairwise_diff(seq1, seq2):
    valid_bases = {"A", "C", "T", "G"}
    differences = 0
    valid_sites = 0
    
    for nt1, nt2 in zip(seq1, seq2):
        if nt1 in valid_bases and nt2 in valid_bases:
            valid_sites += 1
            if nt1 != nt2:
                differences += 1

    return differences, valid_sites

def interquartile_range(data):
    lower_q = np.percentile(data, 25)
    upper_q = np.percentile(data, 75)
    return lower_q, upper_q

def plot_histogram(data, output_file):
    plt.hist(data, bins=20, edgecolor='black', alpha=0.7)
    plt.xlabel("Pairwise Nucleotide Differences")
    plt.ylabel("Frequency")
    plt.title("Histogram of Pairwise Differences")
    plt.savefig(output_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Calculate median pairwise nucleotide difference from an alignment file.")
    parser.add_argument("alignment_file", help="Path to the alignment file (FASTA format).")
    parser.add_argument("-o", "--output", help="Output file name (default: pairwise_nucleotide_diff_results.txt)", default="pairwise_nucleotide_diff_results.txt")
    parser.add_argument("--data_output", help="Output file for raw pairwise differences (default: pairwise_differences.txt)", default="pairwise_differences.txt")
    
    args = parser.parse_args()
    
    alignment = AlignIO.read(args.alignment_file, "fasta")
    sequences = {record.id: str(record.seq) for record in alignment}
    
    pairwise_diffs = []
    valid_sites_list = []
    total_valid_sites = 0
    pairwise_data = []
    
    for i, ((name1, seq1), (name2, seq2)) in enumerate(combinations(sequences.items(), 2)):
        diff, valid_sites = pairwise_diff(seq1, seq2)
        if valid_sites > 0:
            pairwise_diffs.append(diff)
            valid_sites_list.append(valid_sites)
            total_valid_sites += valid_sites
            pairwise_data.append(f"{name1} - {name2}: Differences={diff}, Valid Sites={valid_sites}")
    
    median_diff = np.median(pairwise_diffs)
    median_divergence = np.median([d / v for d, v in zip(pairwise_diffs, valid_sites_list)]) if valid_sites_list else 0
    
    if len(pairwise_diffs) > 1:
        lower_q, upper_q = interquartile_range(pairwise_diffs)
    else:
        lower_q, upper_q = median_diff, median_diff
    
    with open(args.output, "w") as f:
        f.write(f"Alignment File: {args.alignment_file}\n")
        f.write(f"Median Pairwise Nucleotide Difference: {median_diff:.4f}\n")
        f.write(f"Interquartile Range (IQR): ({lower_q:.4f}, {upper_q:.4f})\n")
        f.write(f"Median Pairwise Divergence (normalized by valid sites): {median_divergence:.4f}\n")
        f.write(f"Total Valid Sites: {total_valid_sites}\n\n")
        f.write("Pairwise Differences:\n")
        f.write("\n".join(pairwise_data))
    
    with open(args.data_output, "w") as f:
        f.write("\n".join(pairwise_data))
    
    histogram_output = args.output.replace(".txt", "_histogram.png")
    plot_histogram(pairwise_diffs, histogram_output)
    
    print(f"Results saved to {args.output}")
    print(f"Raw pairwise differences saved to {args.data_output}")
    print(f"Histogram saved to {histogram_output}")

if __name__ == "__main__":
    main()
