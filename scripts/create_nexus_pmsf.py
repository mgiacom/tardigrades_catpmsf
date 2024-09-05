import csv
from collections import defaultdict 
import sys

my_ssf = sys.argv[1]          #site-specific frequencies file 
my_rate = sys.argv[2]         #site-specific rate file (-wsr flag in iqtree)
my_tree_length = sys.argv[3]  #tree length according to the .iqtree output

my_dict = defaultdict(list)

# Load frequencies into the dictionary
with open(my_ssf, 'r') as frequencies:
    reader = csv.reader(frequencies, delimiter=' ')
    for row in reader:
        my_dict[row[0]].append(row[1:])

# Load rates into the dictionary
with open(my_rate, 'r') as rates:
    reader = csv.reader(rates, delimiter='\t')
    for row in reader:
        my_dict[row[0]].append(float(row[1]))

with open(my_rate, 'r') as rates:
    reader = csv.reader(rates, delimiter='\t')
    for row in reader:
        my_dict[row[0]].append(float(row[3]))

print('#nexus'+'\n'+'begin sets;')

# Print charset site declaration
for k, v in my_dict.items():
    print(f"    charset site_{k} = {k};")

print('    charpartition mine = ')

# Print the Poisson+SSF strings
keys = list(my_dict.keys())
for idx, k in enumerate(keys):
    v = my_dict[k]
    freq = '/'.join(v[0])
    gamma = str(v[2])
    tree_len = str(v[1] / my_tree_length)

    # Check if it's the last element
    if idx == len(keys) - 1:
        print(f"        Poisson+F{{{freq}}}+G{{{gamma}}}:site_{k}{{{tree_len}}};")
    else:
        print(f"        Poisson+F{{{freq}}}+G{{{gamma}}}:site_{k}{{{tree_len}}},")

print('end;')
