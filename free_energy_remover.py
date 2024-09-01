'''
This piece of  code removes free_energy from the xyz file
'''
def remove_free_energy_from_dataset(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if 'free_energy' in line:
                # Remove the free_energy field
                parts = line.split()
                filtered_parts = [part for part in parts if not part.startswith('free_energy')]
                new_line = ' '.join(filtered_parts) + '\n'
                outfile.write(new_line)
            else:
                outfile.write(line)


input_file = '/{the directory of the input file}'
output_file = '/{the directory of the output file}'
remove_free_energy_from_dataset(input_file, output_file)
print("removed successfully!")
