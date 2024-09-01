'''
This piece of code converts the traj file to xyz file.
'''
from ase.io import read, write
from ase.io.trajectory import Trajectory
import numpy as np
import matplotlib.pyplot as plt
from ase.neighborlist import neighbor_list
traj_path = "/{the directory of the traj file}"


traj = Trajectory(traj_path)
write(
    "/{name of the file that we want to save the converted xyz file}",
   traj,
    format="extxyz",
)

print("Converted the file!")




