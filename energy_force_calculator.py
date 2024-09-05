from ase.io import read, write
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.calculator import all_changes
from ase.atoms import Atoms
from ase.data import atomic_numbers, atomic_masses
import numpy as np
import matplotlib.pyplot as plt
from ase.optimize import FIRE
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase import units
from mace.calculators.mace import MACECalculator
from time import perf_counter
from numpy.typing import NDArray
from numpy import int32
from ase.calculators.calculator import Calculator, all_changes
from ase.geometry.analysis import Analysis
from ase.ga.utilities import get_rdf
import subprocess

class LiFCalculator(LAMMPSlib):
    implemented_properties = ["energy", "energies", "forces"]
    default_parameters = dict(
        atom_types=None,
        atom_type_masses=None,
        log_file=None,
        lammps_name="",
        keep_alive=True,
        lammps_header=[
            "units metal",
            "atom_style atomic",
            "atom_modify map array sort 0 0",
        ],
        amendments=None,
        post_changebox_cmds=None,
        boundary=True,
        create_box=True,
        create_atoms=True,
        read_molecular_info=False,
        comm=None,
    )

    def __init__(self, eim_file: str, *args, **kwargs):
        pairstyle_cmd = f"pair_style hybrid/overlay eim lj/cut 5.0"
        paircoeffs_cmd = [f"pair_coeff * * eim Li F {eim_file} Li F", "pair_coeff 1 1 lj/cut 1.0 1.0", "pair_coeff 1 2 lj/cut 1.0 1.0", "pair_coeff 2 2 lj/cut 1.0 1.0"]
        append_cmds = ["neighbor 2.0 bin"]
        cmds = [pairstyle_cmd] + paircoeffs_cmd + append_cmds
        super().__init__(
            lmpcmds=cmds,
            keep_alive=True,
        )

    def calculate(self, atoms=None, properties=["energy", "forces"], system_changes=all_changes):
        tic = perf_counter()
        super().calculate(atoms, properties, system_changes)
        toc = perf_counter()

class DualMACECalculator(Calculator):
    implemented_properties = ["energy", "energies", "forces", "free_energy", "rdf"]
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
        # Attach the first calculator to the atoms object and calculate properties
        
        atoms.set_calculator(self.calc1)
        energy1 = atoms.get_potential_energy()
        forces1 = atoms.get_forces()
        #print(f"force1 include {forces1}")
        # Attach the second calculator to the atoms object and calculate properties
        atoms.set_calculator(self.calc2)
        energy2 = atoms.get_potential_energy()
        forces2 = atoms.get_forces()
        #print(f"forces2 include {forces2}")
       
        properties = ["energy", "energies", "forces", "free_energy", "rdf"]
        # Calculate the sum of energies and forces
        total_energy = energy1 + energy2
        total_forces = forces1 + forces2
        self.results['energy'] = total_energy
        self.results['forces'] = total_forces
        print(f"Total forces {self.results['forces']}")
        
        print(f"The properties contains the following: {properties}.")
        return self.results


if __name__ == "__main__":
    at = read("/dataset file path (.xyz) files", "0", format="extxyz")
    
    ################################ Using two MACE Calculators
    # Define paths to the two models
    model_path1 = "/model 1 path"
    model_path2 = "/model 2 path"  # Example path for the second model
    device = "cuda"  # Ensure we are using the correct device
 
    dual_calc = DualMACECalculator(model_path1, model_path2, device)

    at.dual_calc = dual_calc
    at.dual_calc.calculate(at)
    print(f"The outcome of dual_calc is {dual_calc.calculate(at)}.")
    print(f"Unrelaxed: Total_energy {dual_calc.results['energy']:.3f} eV, \n")
  
   
    at.calc = dual_calc
    dyn = FIRE(at, trajectory="/relaxed trajectory (.traj) file")
    dyn.run(fmax=5e-4)
    print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
    write("/output file (.xyz) file", at, format="extxyz")
    # quit()
    T = 300  # MD Temperature in K
    dt = 1  # Timestep in fs
    ntimesteps = int(1e4)  # No. of timesteps
    MaxwellBoltzmannDistribution(at, temperature_K=T)
    dyn = VelocityVerlet(at, dt * units.fs, trajectory="/relaxed file (.traj) file")  # 5 fs time step.
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn.run(ntimesteps)
 
    ############################### Using LAMMPS Calculator
    calc = LiFCalculator("/global/cfs/cdirs/m3606/akash/ml_resources_jornada_group/lammps/potentials/ffield.eim",)
    at.calc = calc
    at.calc.calculate(at)
    print(f"Unrelaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
    dyn = FIRE(at, trajectory="/relaxed trajectory (.traj) file")
    dyn.run(fmax=5e-4)

    write("/putput file (.xyz) file", at, format="extxyz")
    print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
    T = 300  # MD Temperature in K
    dt = 1  # Timestep in fs
    ntimesteps = int(1e4)  # No. of timesteps
    MaxwellBoltzmannDistribution(at, temperature_K=T)
   
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(at, dt * units.fs, trajectory="/unrelaxed trajectory (.traj) file")  
    dyn.run(ntimesteps)
    print("Completed running VelocityVerlet!")
  
    ############################## Using one MACE Calculator
    model_path = "/model path"
    device = "cuda"
    calc = MACECalculator(model_path, device)
    at.calc = calc
    at.calc.calculate(at)
   
    #print(f"The rdf from calc1 :{rdf_list1}")
    print(f"Unrelaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
  
    dyn = FIRE(at, trajectory="/relaxed trajectory (.traj file)")
    dyn.run(fmax=5e-4)
    print(f"Relaxed: Total_energy {at.calc.results['energy']:.3f} eV, \n")
    write("/output_file directory (.xyz) file", at, format="extxyz")
    # quit()
    T = 300  # MD Temperature in K
    dt = 1  # Timestep in fs
    ntimesteps = int(1e4)  # No. of timesteps
    MaxwellBoltzmannDistribution(at, temperature_K=T)
    #mace calculation of rdf    
    dyn = VelocityVerlet(at, dt * units.fs, trajectory="/Unrelaxed trajectory (.traj) file")  
    dyn.run(ntimesteps)
    print("Completed running velocityverlet!")
 
 