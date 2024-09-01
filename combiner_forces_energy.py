'''
This piece of code allows to take the forces and the energies from each file and add them up and write those summed up forces and energies in 
a third file in which positions of the atoms will remain constant.
'''
def write_extxyz(data, filename):
    with open(filename, 'w') as file:
        for entry in data:
            file.write(f"{entry['num_atoms']}\n")
            # Write the header with the new total energy
            new_header = []
            for part in entry['header'].split():
                if part.startswith("energy="):
                    new_header.append(f"energy={entry['total_energy']}")
                else:
                    new_header.append(part)
            file.write(f"{' '.join(new_header)}\n")
            for (atom, pos), force, energy in zip(entry['atoms'], entry['forces'], entry['energies']):
                # Format the position and force values to align them into columns
                pos_str = "  ".join(f"{p:14.8f}" for p in pos)
                force_str = "  ".join(f"{f:14.8f}" for f in force)
                energy_str = f"{energy:14.8f}"
                file.write(f"{atom:<2}   {pos_str}   {energy_str}   {force_str}\n")

def parse_extxyz(content):
    data = []
    lines = content.split('\n')
    i = 0
    while i < len(lines):
        if lines[i].strip():  # Skip empty lines
            num_atoms = int(lines[i].strip())
            header = lines[i + 1].strip()
            atoms = []
            forces = []
            energies = []
            for j in range(num_atoms):
                parts = lines[i + 2 + j].split()
                atom = parts[0]
                pos = list(map(float, parts[1:4]))
                energy = float(parts[4])
                force = list(map(float, parts[5:8]))
                atoms.append((atom, pos))
                forces.append(force)
                energies.append(energy)
            data.append({
                'num_atoms': num_atoms,
                'header': header,
                'atoms': atoms,
                'forces': forces,
                'energies': energies,
            })
            i += num_atoms + 2
        else:
            i += 1
    return data

def merge_datasets(data1, data2):
    merged_data = []
    for entry1, entry2 in zip(data1, data2):
        if entry1['num_atoms'] != entry2['num_atoms']:
            raise ValueError("The number of atoms in the corresponding entries do not match.")
        
        merged_entry = {
            'num_atoms': entry1['num_atoms'],
            'header': entry1['header'],  # assuming headers are identical except for energy
            'atoms': entry1['atoms'],
            'forces': [],
            'energies': [],
            'total_energy': 0.0,
        }
        
        for force1, force2 in zip(entry1['forces'], entry2['forces']):
            merged_force = [f1 + f2 for f1, f2 in zip(force1, force2)]
            merged_entry['forces'].append(merged_force)
        
        merged_entry['energies'] = [e1 + e2 for e1, e2 in zip(entry1['energies'], entry2['energies'])]
        merged_entry['total_energy'] = sum(merged_entry['energies'])  # Sum of the total energies

        merged_data.append(merged_entry)
    
    return merged_data

# Load the datasets
with open('/pscratch/sd/t/tsedalya/second_tutorial/generated_data/lj_data.xyz', 'r') as file:
    content1 = file.read()

with open('/pscratch/sd/t/tsedalya/second_tutorial/generated_data/normal_data.xyz', 'r') as file:
    content2 = file.read()

# Parse the datasets
data1 = parse_extxyz(content1)
data2 = parse_extxyz(content2)

# Merge the datasets
merged_data = merge_datasets(data1, data2)

# Write the merged dataset to a new file
write_extxyz(merged_data, '/{insert the name of the directory that you want to save the file to}')

print("Merged dataset created and saved as '{name_file.xyz}'.")
