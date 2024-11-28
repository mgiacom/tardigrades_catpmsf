import sys

# Read data from the first file (alignment with gaps)
alignment_with_gaps = []
with open(sys.argv[1], "r") as file:
    for line in file:
        if line.startswith(">"):
            header = line.strip()
            alignment_with_gaps.append({"header": header, "sequence": ""})
        else:
            alignment_with_gaps[-1]["sequence"] += line.strip()

# Read data from the second file (alignment without gaps)
alignment_without_gaps = []
with open(sys.argv[2], "r") as file:
    for line in file:
        if line.startswith(">"):
            header = line.strip()
            alignment_without_gaps.append({"header": header, "sequence": ""})
        else:
            alignment_without_gaps[-1]["sequence"] += line.strip()

# Step 1: Identify Gap Positions
gap_positions = []
for aligned_entry in alignment_with_gaps:
    gaps = [i for i, aa in enumerate(aligned_entry["sequence"]) if aa == "-"]
    gap_positions.append(gaps)

# Step 2: Transfer Gaps
alignment_with_gaps_transferred = []
for i, aligned_entry in enumerate(alignment_without_gaps):
    new_aligned_seq = list(aligned_entry["sequence"])
    for gap_pos in gap_positions[i]:
        new_aligned_seq[gap_pos] = '-'  # Replace amino acid with gap
    alignment_with_gaps_transferred.append("".join(new_aligned_seq))

# Write aligned sequences with transferred gaps into a new FASTA file
with open(sys.argv[3], "w") as file:
    for i, seq_entry in enumerate(alignment_with_gaps_transferred):
        file.write(f"{alignment_with_gaps[i]['header']}\n")
        file.write(f"{seq_entry}\n")
