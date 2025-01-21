import gzip
import numpy as np
from Bio import SeqIO
import subprocess
import os


# Function to read the FASTQ file
def read_fastq(file_path):
    with gzip.open(file_path, "rt") if file_path.endswith("gz") else open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield record


# Function to compute average quality score of the sequences in a FASTQ file
def compute_average_quality(file_path):
    total_quality = 0
    total_bases = 0
    for record in read_fastq(file_path):
        quality_scores = record.letter_annotations["phred_quality"]
        total_quality += sum(quality_scores)
        total_bases += len(quality_scores)

    average_quality = total_quality / total_bases if total_bases > 0 else 0
    return average_quality


# Function to compute base composition (A, T, G, C) in the sequences of a FASTQ file
def compute_base_composition(file_path):
    base_counts = {"A": 0, "T": 0, "G": 0, "C": 0}
    total_bases = 0
    for record in read_fastq(file_path):
        seq = record.seq
        for base in seq:
            if base in base_counts:
                base_counts[base] += 1
                total_bases += 1

    base_composition = {base: count / total_bases for base, count in
                        base_counts.items()} if total_bases > 0 else base_counts
    return base_composition


# Function to compute sequence length distribution
def compute_length_distribution(file_path):
    lengths = []
    for record in read_fastq(file_path):
        lengths.append(len(record.seq))
    return np.histogram(lengths, bins=50)


# Run FastQC tool using subprocess
def run_fastqc(file_path, output_dir="fastqc_reports"):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Run FastQC and save the results to the specified directory
    subprocess.run(["fastqc", file_path, "-o", output_dir])


# Main function to perform quality check
def perform_quality_check(file_path):
    print("Performing basic quality checks...\n")

    # Average Quality Score
    avg_quality = compute_average_quality(file_path)
    print(f"Average Quality Score: {avg_quality}")

    # Base Composition
    base_composition = compute_base_composition(file_path)
    print(f"Base Composition (A, T, G, C): {base_composition}")

    # Sequence Length Distribution
    lengths, bins = compute_length_distribution(file_path)
    print(f"Sequence Length Distribution (first 10 counts): {lengths[:10]}")

    # print("\nRunning FastQC for detailed analysis...")
    # run_fastqc("F:\DEEPTHI\SELF STUDY\fastqc_reports\")  # You can also specify an output directory for the reports


# Usage example: Replace with your actual FASTQ file path
fastq_file = "F:\DEEPTHI\SELF STUDY\PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq"  # Example FASTQ file path

perform_quality_check(fastq_file)
