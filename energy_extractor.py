import pandas as pd
import re

# Specify the path to your file
file_path = '/{the directory of the file from which we are going to extract energy}'

# Read the data
with open(file_path, 'r') as file:
    lines = file.readlines()

# Create a DataFrame to store the extracted data
data = {
    'energy': []
}

# Regular expression pattern to extract MACE_energy values
pattern = re.compile(r'energy=(-?\d+\.\d+)')

# Process each line to find MACE_energy
i = 0
while i < len(lines):
    num = int(lines[i])
    match = pattern.search(lines[i + 1])
    if match:
        # Append the extracted value to the list
        data['energy'].append(float(match.group(1)))
    else:
        # Append NaN if MACE_energy is not found in the line
        data['energy'].append(None)
    i += num + 2
    

# Create a DataFrame from the extracted data
df = pd.DataFrame(data)

# Display the DataFrame
print(df)

# Optionally, save the DataFrame to a new CSV file
output_file_path = '/{the directory in to which we are going to save the output file}'

df.to_csv(output_file_path,header = False, index=False)
