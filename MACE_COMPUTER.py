'''
This piece of code uses MACE calculator to calculate kinetic energy, potential energy, radial distribution functiona and mean-squared displacement.
'''
###################################### The one below calculates the kinetic energy and the potential energy using MACE
from ase import Atoms
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from mace.calculators.mace import MACECalculator
from ase.units import kB
from ase.geometry.rdf import get_rdf 
from ase import units
import csv
import numpy as np

# Set up your atomic system
atoms = read('/{the directory of the dataset file }', "0", format="extxyz")

# Initialize the MACE calculator
calculator = MACECalculator(model_paths='/{the MACE model path}', device = "cuda")
atoms.calc = calculator
atoms.set_calculator(atoms.calc)

# Initialize velocities at a given temperature 
MaxwellBoltzmannDistribution(atoms, temperature_K=0)

# Define the molecular dynamics method
dyn = VelocityVerlet(atoms, dt=1.0 * units.fs)
pe_list = []
ke_list = []
def print_kinetic_energy(atoms):
    ekin = atoms.get_kinetic_energy()  # in eV
    epot = atoms.get_potential_energy()
    pe_list.append(epot)
    ke_list.append(ekin)
# Attach the print function to the MD engine
dyn.attach(print_kinetic_energy, interval=1, atoms=atoms)

# Run the dynamics
dyn.run(steps=1000)
csv_file_path = '/{the potential energy output file}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in pe_list:
        writer.writerow([row])
csv_file_path2 = '/{the kinetic energy output file}'
with open(csv_file_path2, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in ke_list:
        writer.writerow([row])
        
print("Save to CSV 1")

######################################################################### The one below calculates the RDF
atoms = read("/{the directory of the dataset file }", "0", format="extxyz")
# Create a MACE calculator instance
mace_calculator = MACECalculator(
    model_paths='/{the MACE model path}',  # Path to your MACE model file
    device="cuda"  # Adjust device if necessary
)
calc = mace_calculator
atoms.calc = calc
# Attach the calculator to the atoms object
atoms.set_calculator(atoms.calc)
MaxwellBoltzmannDistribution(atoms, temperature_K=1000)
# Set up the MD simulation
timestep = 1.0  # Time step in fs
steps = 1000  # Number of MD steps
# Compute RDF after MD simulation
rmax = 5  # Maximum distance for RDF calculation
nbins = 100  # Number of bins for RDF calculation

dyn = VelocityVerlet(atoms, dt=timestep* units.fs)
list1 = list()

def rdf_calculation(atoms):
    rdf_result = get_rdf(atoms, rmax=rmax, nbins=nbins)
    member=np.array(rdf_result[0]).tolist()
    list1.append(member)
    
dyn.attach(rdf_calculation, interval=1, atoms=atoms)

dyn.run(steps=steps) 

csv_file_path = '/{the rdf output file}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    first_row = np.arange(0, 100)
    writer.writerow(first_row)
    for row in list1:
        writer.writerow(row)
print(f"Data successfully written to the csv 2")

######################################################  The one below calculates the mean-squared displacement.
# Define your atomic configuration
atoms = read("{the directory of the dataset file}", "0", format="extxyz")
# Create a MACE calculator instance
mace_calculator = MACECalculator(
    model_paths='/{the MACE model path}',  # Path to your MACE model file
    device= "cuda"
)
calc = mace_calculator
atoms.calc = calc
# Attach the calculator to the atoms object
atoms.set_calculator(atoms.calc)
temperature = 0  # Temperature in Kelvin
dt = 1
ntimesteps = 100 # Time step in fs
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
dyn = VelocityVerlet(atoms, ntimesteps)



msd_values = []
# Store initial positions for MSD calculation
initial_positions = atoms.get_positions()
msd_list = []
def calculate_msd(atoms):
    current_positions = atoms.get_scaled_positions()

    # Ensure the shapes of initial_positions and current_positions are correct
    if initial_positions.shape != current_positions.shape:
        raise ValueError("Initial and current positions must have the same shape")

    # Calculate the displacement in fractional coordinates
    displacement = current_positions - initial_positions

    # Apply PBC by wrapping displacements into the range [-0.5, 0.5)
    displacement = displacement - np.round(displacement)

    # Convert displacement back to Cartesian coordinates
    displacement_cartesian = np.dot(displacement, atoms.cell[:3, :3])

    # Ensure displacement_cartesian is 2D
    if displacement_cartesian.ndim == 1:
        displacement_cartesian = displacement_cartesian[:, np.newaxis]

    # Compute the squared displacement
    squared_displacement = np.sum(displacement_cartesian ** 2, axis=1)

    # Average over all atoms for this timestep
    msd = np.mean(squared_displacement)

    # Append to MSD values
    msd_values.append(msd)
   

# Initialize arrays to store MSD values
dyn.attach(calculate_msd, interval=1, atoms = atoms)
dyn.run(ntimesteps)
msd_list = np.array(msd_values).tolist()
print(msd_list)   
csv_file_path = '/{the msd output file}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in msd_list:
        print(row)
        writer.writerow([row])
print(f"Data successfully written to the csv 3")
