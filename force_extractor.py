'''
This piece of code takes out the forces in the xyz file and writes them in another file.
'''
import numpy as np

# Specify the path to your .extxyz file
extxyz_file_path = '{the directory of the input file}'

# Initialize list to store data
result = []

# Read and process the file
with open(extxyz_file_path, 'r') as file:
    lines = file.readlines()

# Iterate through lines to process data
i = 0
while i < len(lines):
    # Number of atoms (this line is generally used to identify the start of a new block)
    num_atoms = int(lines[i].strip())
    i += 1

    # Skip the comment line
    i += 1

    # Extract the atom data
    for _ in range(num_atoms):
        values = lines[i].split()
        if len(values) >= 4:  # Example: atom type, x, y, z, possibly more
            # Extract coordinates and convert to floats
            coordinates = [float(x) for x in values[5:8]]  # Adjust indexing if more columns are present
            result.append(coordinates)
        i += 1

    # Handle additional data if present
    # Here you would add code to handle extra data if needed
    # For example, if you have energy values or forces, you can extract them similarly

# Convert the result list to a NumPy array
result_array = np.array(result)

# Write the NumPy array to a CSV file
output_file_path = '{the directory of the output file}'
np.savetxt(output_file_path, result_array, delimiter=",", fmt="%.6f")
print("completed saving!")
