'''
This piece of code takes the the lammps test file and evaluates it. The output is an xyz file that evaluates as per the trained model
'''
import subprocess

command = [
    "mace_eval_configs",
    "--configs={the directory of the lammps test file}",
    "--model=/{the diretory of the trained model}",
    "--output=/{the directory of the mace output xyz file}"
]

# Run the command
result = subprocess.run(command, capture_output=True, text=True)

# Print the output
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)
print("evaluation completed!")