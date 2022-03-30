import os, sys
from io import open

data_path = "sample.out"
read_flag = False
n = 0
with open(data_path) as of:
    for line in of.readlines():
        if "Standard orientation" in line:
            read_flag = True
            atomic_data = []
            n = 0
        elif "-----------------" in line and n > 4:
            read_flag = False
        else:
            n += 1

        if n>4 and read_flag:
            atomic_data.append(line)

for atom in atomic_data:
    print(atom)
