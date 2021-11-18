import numpy as np
import pandas as pd

# Name of file
name = "LLZO"
# Size of supercell
supersize = 2
new_name = name + "_" + str(supersize)

# Open file and get lines
structure_file = open(name+".vasp", 'r')
lines = structure_file.readlines()
structure_file.close()
scale = lines[1]
lattice = lines[2:5]
composition = lines[5:7]
method = lines[7]
atomic_positions = lines[8:]

# Set new scale and copy original lattice
new_scale = float(scale) * supersize
new_lattice = lattice

# Get multiply atom quantity by supercell size ** 3, then sort by Element alphabetically
new_composition = ["", ""]
new_composition[0] = np.array(composition[0].rstrip('\n').split(" "))
new_composition[1] = np.array(composition[1].split(" "), dtype="int64")
new_composition[1] *= supersize**3
new_composition = np.array(new_composition)
new_composition = pd.DataFrame(new_composition.transpose(), columns=['Element', 'Count'])
new_composition = new_composition.sort_values(by=['Element'])
new_composition = new_composition.to_numpy().transpose()

# Copy method
new_method = method

# Split lines into an array of atomic positions
new_atomic_positions = np.empty((0, 3))
split = [x.split(" ") for x in atomic_positions]
split = np.array(split)
# Add offset matrix to fill supercell and concatenate data to the last
for i in range(supersize):
    for j in range(supersize):
        for k in range(supersize):
            x = np.array([[(1/supersize)*i, 0, 0] for x in range(len(atomic_positions))])
            y = np.array([[0, (1/supersize)*j, 0] for x in range(len(atomic_positions))])
            z = np.array([[0, 0, (1/supersize)*k] for x in range(len(atomic_positions))])
            sum = x+y+z
            sum += (1/supersize) * split[:, 0:3].astype('float64')
            new_atomic_positions = np.concatenate( (new_atomic_positions, sum), axis=0)
# Reshape Element names to concatenate on right hand side
temp = np.array([np.reshape(split[:,3], split[:,3].size) for x in range(supersize**3)])
temp = np.reshape(temp, temp.size)
temp = np.reshape(temp, (temp.size, 1))
new_atomic_positions = np.concatenate( (new_atomic_positions, temp), axis=1)
# Sort by element alphabetically
new_atomic_positions = pd.DataFrame(new_atomic_positions, columns=['X', 'Y', 'Z', 'Element'])
new_atomic_positions = new_atomic_positions.sort_values(by=['Element'])
new_atomic_positions = new_atomic_positions.to_numpy()

# Write to new file
with open(new_name+".vasp", "w") as f:
    f.write(new_name+"\n")
    f.write(str(new_scale)+"\n")
    for line in new_lattice:
        f.write(line)
    f.write(" ".join(new_composition[0]) + "\n")
    f.write(" ".join(map(str, new_composition[1])) + "\n")
    f.write(new_method)
    for line in new_atomic_positions:
        f.write(" ".join(line))
