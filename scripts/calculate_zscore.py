from Bio import AlignIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import sys
import math

# Read the original dataset
original_dataset = sys.argv[1]
with open(original_dataset, "r") as file:
    original_sequences = {title: seq for title, seq in SimpleFastaParser(file)}

# Extract alignment information
original_taxa = list(original_sequences.keys())
original_alignment_length = len(next(iter(original_sequences.values())))

# Calculate mean AA diversity per site for the original dataset
original_total_diversity = 0
for i in range(original_alignment_length):
    site_states = set(
        seq[i] for seq in original_sequences.values() if seq[i] not in {"-", "X", "?"}
    )
    original_total_diversity += len(site_states)

original_mean_diversity = original_total_diversity / original_alignment_length


files = [file for file in os.listdir(".") if file.endswith(".fa_gaps")]
simulated_diversities = []

for infile in files:
    print(f"Processing {infile}")
    with open(infile, "r") as file:
        simulated_sequences = {title: seq for title, seq in SimpleFastaParser(file)}

    alignment_length = len(next(iter(simulated_sequences.values())))
    total_diversity = 0

    for i in range(alignment_length):
        site_states = set(
            seq[i] for seq in simulated_sequences.values() if seq[i] not in {"-", "X", "?"}
        )
        total_diversity += len(site_states)

    mean_diversity = total_diversity / alignment_length
    simulated_diversities.append(mean_diversity)


# Calculate statistics for simulated datasets
average_simulated_diversity = sum(simulated_diversities) / len(simulated_diversities)
variance = sum(
    (diversity - average_simulated_diversity) ** 2 for diversity in simulated_diversities
) / (len(simulated_diversities) - 1)
standard_deviation = math.sqrt(variance)

z_score = -(original_mean_diversity - average_simulated_diversity) / standard_deviation

# Write output
with open("diversity.pbr_gaps", "w") as out:
    out.write(
        "Test of model adequacy / pattern heterogeneity / across sites compositional heterogeneity\n"
    )
    out.write(f"Diversity original data: {original_mean_diversity}\n")
    out.write(f"Average diversity simulated data: {average_simulated_diversity}\n")
    out.write(f"SD simulated data: {standard_deviation}\n")
    out.write(f"Z-score: {z_score}\n")

with open("diversity_scores_bootstrapped_data.txt_gaps", "w") as out1:
    out1.write("file \t diversity\n")
    for file, diversity in zip(files, simulated_diversities):
        out1.write(f"{file}\t{diversity}\n")
