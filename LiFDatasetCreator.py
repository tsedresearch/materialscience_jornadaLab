from ase.atoms import Atoms
from ase.io import read, write
from ase.geometry.cell import Cell, cellpar_to_cell
from ase.build import make_supercell
from ase.data import covalent_radii
import numpy as np
from numpy.random import normal
from numpy.typing import NDArray
from numpy import (
    float_,
    int_,
    vstack,
    argmax,
    argmin,
    array,
    arange,
    allclose,
    dot,
    eye,
    zeros,
    round,
    diag,
    random,
)
from numpy.linalg import norm
from typing import List, Optional, Dict
from fractions import Fraction
from math import cos, pi, acos


from itertools import combinations, combinations_with_replacement, permutations, product

IntArray = NDArray[int_]
FloatArray = NDArray[float_]

# Set Seed
SEED = 0
random.seed = SEED
def generate_random_structures_ratlle(
    base_at: Atoms = read("/{the directory of the structure file}", format = "vasp"),
    N: int = 1,  # Number of Structures to generate
    max_disp: float = 0.1,
    max_strain: float = 0.05,  # In fraction units
    max_angle_strain: float = 0.01,  # In fraction units
    supercell_matrix: Optional[IntArray] = None,
    is_2d = True
):
    """Generate random structures by rattling

    Args:
        base_at (Atoms): Base structure as ASE atoms object
        N (int, optional): Number of structures to generate. Defaults to 100.
        max_disp (float, optional): Maximum displacement of atoms in Angstroms. Defaults to 0.1.
        max_strain (float, optional): Maximum strain of cell lengths. Defaults to 0.05.
        max_angle_strain (float, optional): Maximum strain of cell lengths. Defaults to 0.01.
        is_2d (bool, optional): If 2d don't strain c, alpha, beta. Defaults to True.

    Returns:
        List[Atoms]: List of Atoms objects
    """
    print(base_at.get_cell())
    if supercell_matrix is not None:
        base_at = make_supercell(base_at, supercell_matrix)
        
        
    cell_params: FloatArray = base_at.cell.cellpar()

    list_of_cands: List[Atoms] = []

    for i in range(N):
        at: Atoms = base_at.copy()
        new_cell_params: FloatArray = cell_params.copy()
        new_cell_params[:3] += new_cell_params[0:3] * normal(
            loc=0, scale=max_strain, size=3
        )

        if is_2d:
            # If 2D don't change c
            new_cell_params[2] = cell_params[2]
        new_cell_params[3:] += new_cell_params[3:] * normal(
            loc=0, scale=max_angle_strain, size=3
        )
        new_cell: Cell = cellpar_to_cell(new_cell_params)
        at.set_cell(new_cell, scale_atoms=True)
        # Replacing Rattle, so we can rely on one seed, can also be done by using seed, tis all the same
        # at.rattle(stdev=max_disp, seed=i)
        positions: FloatArray = at.arrays["positions"]
        at.set_positions(positions + normal(scale=max_disp, size=positions.shape))
        at.center()  # Add centering so dipole position can be consistent
        list_of_cands.append(at)
    print(f"supercell_matrix {supercell_matrix}")

    return list_of_cands

candidates = generate_random_structures_ratlle(supercell_matrix=2*np.eye(3))

write("/{the output file directory}", candidates, format='extxyz')
print("Wrote all candidates to all_candidates.xyz")
