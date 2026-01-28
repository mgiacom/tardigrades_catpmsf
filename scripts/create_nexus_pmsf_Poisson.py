import csv
from collections import defaultdict 
import sys

# File inputs
my_ssf = sys.argv[1]          # site-specific frequencies file 
my_rate = sys.argv[2]         # site-specific rate file (-wsr flag in iqtree)

my_dict = defaultdict(list)

# Load frequencies into the dictionary
with open(my_ssf, 'r') as frequencies:
    reader = csv.reader(frequencies, delimiter=' ')
    for row in reader:
        my_dict[row[0]].append(row[1:])  # Store frequencies as a list

# Load rates into the dictionary
with open(my_rate, 'r') as rates:
    reader = csv.reader(rates, delimiter='\t')
    for row in reader:
        site = row[0]
        rate = float(row[1])
        gamma = float(row[3])

        # Add both rate and gamma to
        my_dict[site].append(rate)
        my_dict[site].append(gamma)

# Begin the output .nex
print('#nexus\nbegin sets;')

# Print charset site declaration
for k in my_dict:
    print(f"    charset site_{k} = {k};")

# Print the charpartition block
print('    charpartition mine = ')

# Iterate through the dictionary and construct the Model+SSF+G+tree_length string for each site
for idx, (k, v) in enumerate(my_dict.items()):
    freq = '/'.join(v[0])        # Frequency values (list of strings)
    tree_len = str(v[1])  # site rate to be scaled with tree length
    gamma = str(v[2])            # Gamma value

    # If it's the last entry, put ';', otherwise ','
    if idx == len(my_dict) - 1:
        print(f"        Poisson+F{{{freq}}}+G{{{gamma}}}:site_{k}{{{tree_len}}};")   #Change the model here below
    else:
        print(f"        Poisson+F{{{freq}}}+G{{{gamma}}}:site_{k}{{{tree_len}}},")

print('end;')
