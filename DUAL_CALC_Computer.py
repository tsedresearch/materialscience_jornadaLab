'''
This piece of code calculates the kinetic energy, potential energy, the radial distribution function and the mean-squared displacement using two different MACE calculators that take two different models.    
'''
################################# The piece of code below calculates the KE and PE
from ase import Atoms
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from mace.calculators.mace import MACECalculator
from ase.calculators.calculator import Calculator, all_changes
from ase.units import kB
from ase.geometry.rdf import get_rdf
import csv
from ase import units
import numpy as np
class DualMACECalculator(Calculator):
    implemented_properties = ["energy", "forces"]
    def __init__(self, model_path1, model_path2, device="cuda"):
        # Initialize the first calculator
        Calculator.__init__(self)
        self.calc1 = MACECalculator(
            model_paths=model_path1,
            device=device,
            energy_units_to_eV=1.0,
            length_units_to_A=1.0,
            default_dtype="",
            charges_key="Qs",
            model_type="MACE"
        )

        # Initialize the second calculator
        self.calc2 = MACECalculator(
            model_paths=model_path2,
            device=device,
            energy_units_to_eV=1.0,
            length_units_to_A=1.0,
            default_dtype="",
            charges_key="Qs",
            model_type="MACE"
        )
        self.results = {}

    def calculate(self, atoms, properties=None, system_changes=all_changes):
        # Use the first calculator
        atoms.set_calculator(self.calc1)
        pe_energy1 = atoms.get_potential_energy()
        ke_energy1 = atoms.get_kinetic_energy()
        force1 = atoms.get_forces()
        # Use the second calculator
        atoms.set_calculator(self.calc2)
        pe_energy2 = atoms.get_potential_energy()
        force2 = atoms.get_forces()
        ke_energy2 = atoms.get_kinetic_energy()
        self.results["energy"] = pe_energy1 + pe_energy2
        self.results["forces"] = force1 + force2
        self.results["KE_energy"] =  ke_energy1 + ke_energy2
        return self.results

        
atoms = read('/{the xyz file directory that contains the atom info}', "0", format="extxyz")

model_path1 = "/{the directory of the trained model one}"
model_path2 = "/{the direcotry of the trained model two}" 
device = "cuda"  # Ensure we are using the correct device
 
dual_calc = DualMACECalculator(model_path1, model_path2, device)

atoms.dual_calc = dual_calc

# Initialize the MACE calculator
atoms.set_calculator(atoms.dual_calc)

# Initialize velocities at a given temperature 
MaxwellBoltzmannDistribution(atoms, temperature_K=1000)

# Define the molecular dynamics method
dyn = VelocityVerlet(atoms, dt=1.0 * units.fs)

pe_list = []
ke_list = []
    # Function to print the kinetic energy at each step
def print_energies(atoms):
    ekin = atoms.get_kinetic_energy()  # in eV
    epot = atoms.get_potential_energy()
    pe_list.append(epot)
    ke_list.append(ekin)

# Attach the print function to the MD engine
dyn.attach(print_energies, interval=1, atoms=atoms)

dyn.run(steps=1000)
csv_file_path = '/{the output}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in pe_list:
        writer.writerow([row])
csv_file_path2 = '/pscratch/sd/t/tsedalya/second_tutorial/generated_data/computed_properties/1000K/ke/dual_calc_ke_1000.csv'
with open(csv_file_path2, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in ke_list:
        writer.writerow([row])
        
################################### The piece of code below calculates the RDF
atoms = read("/pscratch/sd/t/tsedalya/second_tutorial/generated_data/LiF_data_latest2.xyz", "0", format="extxyz")

model_path1 = "/{The file directory of MACE model one}"
model_path2 = "/{The file directory of MACE model two}" 
device = "cuda"  # Ensure we are using the correct device
 
dual_calc = DualMACECalculator(model_path1, model_path2, device)

atoms.dual_calc = dual_calc

# Initialize the MACE calculator
atoms.set_calculator(atoms.dual_calc)
MaxwellBoltzmannDistribution(atoms, temperature_K=1000)

timestep = 1.0  # Time step in fs
steps = 1000  # Number of MD steps
# Compute RDF after MD simulation
rmax = 5  # Maximum distance for RDF calculation
nbins = 100  # Number of bins for RDF calculation

dyn = VelocityVerlet(atoms, dt=1.0 * units.fs)
list1 = list()
def rdf_calculation(atoms):
    rdf_result = get_rdf(atoms, rmax=rmax, nbins=nbins)
    member=np.array(rdf_result[0]).tolist()
    list1.append(member)
    
dyn.attach(rdf_calculation, interval=1, atoms=atoms)

dyn.run(steps=steps) 

csv_file_path = '/{the output directory that saves the calculated rdf}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    first_row = np.arange(0, 100)
    writer.writerow(first_row)
    for row in list1:
        writer.writerow(row)
print(f"Data successfully written to the csv")


######################################## The code below calculates the mean-squared displacement (MSD)
atoms = read("/{the file directory that has}", "0", format="extxyz")

model_path1 = "/{the directory of the trained model one}"
model_path2 = "/{the directory of the trained model two}" 
device = "cuda"  # Ensure we are using the correct device
 
dual_calc = DualMACECalculator(model_path1, model_path2, device)

atoms.dual_calc = dual_calc

# Initialize the MACE calculator
atoms.set_calculator(atoms.dual_calc)
MaxwellBoltzmannDistribution(atoms, temperature_K=0)
# Set up the MD simulation
timestep = 1.0  # Time step in fs
steps = 1000  # Number of MD steps
# Compute RDF after MD simulation
dyn = VelocityVerlet(atoms, dt=1.0 * units.fs)

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
dyn.run(steps)
msd_list = np.array(msd_values).tolist()
print(msd_list)   
csv_file_path = '/pscratch/sd/t/tsedalya/second_tutorial/generated_data/computed_properties/0K/msd/dual_calc_msd_0.csv'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in msd_list:
        print(row)
        writer.writerow([row])
print(f"Data successfully written to the csv")
