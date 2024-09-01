'''
This piece of code calculates the hybrid Embedded-ion (eim) potential energy and forces and Lennard-jonnes(LJ) potential energy and forces. It calculates the kinetic energy, potential energy, radial distribution function and mean-squared displacement.
'''
import ase
from ase.io import read, write
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.calculator import all_changes
from ase.atoms import Atoms
from ase.data import atomic_numbers, atomic_masses
import numpy as np


from numpy.typing import NDArray
from numpy import int32

from time import perf_counter

from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from ase.optimize import FIRE
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.geometry.rdf import get_rdf
import csv

class LiFCalculator(LAMMPSlib):
    implemented_properties = ["energy", "energies", "forces"]

    def __init__(self, eim_file: str, data_file: str = "{the data file directory which has the position of atoms}", *args, **kwargs):
        self.units = "metal"  # Setting units for the simulation

        lmp_commands = [
            "units metal",
            "atom_style atomic",
            "boundary p p p",
            "dimension 3",
            f"read_data {data_file}",
        ]

        atom_types = {
            "Li": 1,  # Lithium is type 1 in LAMMPS
            "F": 2    # Fluorine is type 2 in LAMMPS
        }
        #check out the lammps documentation to change the parameters as needed in the lammps command.
        lmp_commands = [ 
            "pair_style hybrid/overlay eim lj/cut 5.0",
            f"pair_coeff * * eim Li F {eim_file} Li F",
            "pair_coeff 1 1 lj/cut 1.0 1.0",
            "pair_coeff 1 2 lj/cut 1.0 1.0",
            "pair_coeff 2 2 lj/cut 1.0 1.0",
            "neighbor 2.0 bin",
        ]

        super().__init__(
            lmpcmds=lmp_commands,
            atom_types=atom_types,
            keep_alive=True,
            *args, **kwargs
        )

    def calculate(self, atoms=None, properties=["energy", "forces"], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        return self.results

if __name__ == "__main__":
    at = read("/{the directory which has the dataset}", "0", format="extxyz")

    calc = LiFCalculator("/global/cfs/cdirs/m3606/akash/ml_resources_jornada_group/lammps/potentials/ffield.eim")
    at.calc = calc
        
    dyn = FIRE(at, trajectory="/{the directory in which the relaxed atom trajectory will be stored.}")
    dyn.run(fmax=5e-4)
    print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
    T = 1000  # MD Temperature in K
    dt = 1  # Timestep in fs
    ntimesteps = 1000
    MaxwellBoltzmannDistribution(at, temperature_K=T)
    dyn = VelocityVerlet(at, dt * units.fs, trajectory="/{the directory in which the velocityVerlet trajectory is saved}")
    list1 = list()
    rmax = 5  # Maximum distance for RDF calculation
    nbins = 100  # Number of bins for RDF calculation

    def rdf_calculation(atoms):
        rdf_result = get_rdf(atoms, rmax=rmax, nbins=nbins)
        member=np.array(rdf_result[0]).tolist()
        list1.append(member)
    
    dyn.attach(rdf_calculation, interval=1, atoms=at)

    dyn.run(steps=ntimesteps) 

    csv_file_path = '/{the output file directory}'
    with open(csv_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        first_row = np.arange(0, 100)
        writer.writerow(first_row)
        for row in list1:
            writer.writerow(row)
    
    print("Finished VerletVelocity training!!! 1")
    
    ########################################### The code below calculates the potential energy and kinetic energy using Lammps calculator
at = read("/{the xyz file that information about the atoms}", "0", format="extxyz")

# It is the directory that has the eim potential parameters (change as needed)
calc = LiFCalculator("/global/cfs/cdirs/m3606/akash/ml_resources_jornada_group/lammps/potentials/ffield.eim")
at.calc = calc

dyn = FIRE(at, trajectory="/{the directory of the trajectory file in which the relaxation is saved}")
dyn.run(fmax=5e-4)
print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
T = 1000  # MD Temperature in K
dt = 1  # Timestep in fs
ntimesteps = 1000
MaxwellBoltzmannDistribution(at, temperature_K=T)
dyn = VelocityVerlet(at, dt * units.fs, trajectory="/{the directory of the trajectory file in which the VelocityVerlet md simulation is saved}")

pe_list = []
ke_list = []
    # Function to print the kinetic energy at each step
def print_energies(atoms):
    # Get results from both calculators
    ekin = atoms.get_kinetic_energy()  # in eV
    epot = atoms.get_potential_energy()
    pe_list.append(epot)
    ke_list.append(ekin)

# Attach the print function to the MD engine
dyn.attach(print_energies, interval=1, atoms=at)

dyn.run(steps=1000)
csv_file_path = '/{the directory of the output file in which the potential energy is saved}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in pe_list:
        writer.writerow([row])
csv_file_path2 = '/{the directory of the outpput file in which the kinetic energy is saved}'
with open(csv_file_path2, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in ke_list:
        writer.writerow([row])
print("Computed the energies!")
        
################################### This code below calculates the mean-sqaured displacement by initiating the Lammps calculator
at = read("/{the file directory that has the xyz information}", "0", format="extxyz")

calc = LiFCalculator("/global/cfs/cdirs/m3606/akash/ml_resources_jornada_group/lammps/potentials/ffield.eim")
at.calc = calc

dyn = FIRE(at, trajectory="/{the directory of the file to save the directory}")
dyn.run(fmax=5e-4)
print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
T = 1000  # MD Temperature in K
dt = 1  # Timestep in fs
ntimesteps = 1000
MaxwellBoltzmannDistribution(at, temperature_K=T)
dyn = VelocityVerlet(at, dt * units.fs, trajectory="/{the directory of the file to save the directory of VelocityVerlet md simulation}")

msd_values = []
# Store initial positions for MSD calculation
initial_positions = at.get_positions()
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
dyn.attach(calculate_msd, interval=1, atoms = at)
dyn.run(steps=ntimesteps )
msd_list = np.array(msd_values).tolist()
print(msd_list)   
csv_file_path = '/{the directory of the output file paths}'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in msd_list:
        print(row)
        writer.writerow([row])
print(f"Data successfully written to the csv")