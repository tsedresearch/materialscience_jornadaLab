# Python Scripts Used
Using machine learning interatomic potentials to study the photodeterioration of materials

 <details>
 <summary> Dataset Creation </summary><br />
        Using the [LiFDatasetCreator.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/LiFDatasetCreator.py), we can create rattled structures of materials. <br />
           => to change the number of structures to be created, assign the "N" argument in the **generate_random_strucutres_ratlle()** function to the number of the rattled strucutres that you want to generate.
  </details>
  
   <details>
  <summary>Conversion of Trajectory Files to XYZ-files</summary><br />
         To convert files that have .traj extension to xyz file, use the [convert_traj.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/convert_traj.py)
   </details>
   
  <details> 
 <summary>Calculating Forces and Energies</summary> <br />
    After creating the rattled structures, to calculate the forces and the energies, use [energy_force_calculator.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/energy_force_calculator.py). This file takes one xyz file and creates three different xyz output file: one calculated with single MACE calculator, one calculated with double MACE calculators, one LAMMPS calculator. To understand more about the types of potentials to select and the parameters to put for LAMMPS calculation, check out the following documentation [LAMMPS documentation](https://docs.lammps.org/).
  </details>     

<details>
 <summary>     
Energy Extractor
 </summary>
   <br />
        This file takes an xyz file as an input and extracts the energy in the file and writes it to a csv file. The file is [energy_extractor.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/energy_extractor.py).
</details>

 <details>
   <summary>
 Force Extractor </summary> <br />
        => This file takes an xyz file as an input and extracts the force in the file and writes it to a csv file. The file is [force_extractor.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/force_extractor.py).
 </details>

 <details>
   <summary>
Training MACE Models</summary><br />
         => To train MACE models, use the [train_mace_model.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/train_mace_model.py). To learn more about the parameters, visit the following [github page](https://github.com/ACEsuit/mace).
 </details>

  <details>
  <summary>
  Evaluating MACE Models </summary><br />
          To evalute the trained MACE models, use [evaluate_mace.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/evaluate_mace.py). You can find the trained model inside the "results" directory specified in the training script.
  <br />
  </details>
  
<details>
<summary>Computing Obervables (Potential Energy, Kinetic Energy, Radial Distribution FUnction, Mean-squared Displacement)</summary><br />
       => To compute with LAMMPS, use [lammps_LiFCalculator.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/lammps_LiFCalculator.py). <br /><br />
       => To compute observables with one MACE model using MACE calculator, use [MACE_Computer.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/MACE_COMPUTER.py). <br /><br />
       => To compute observables with two MACE models using MACE calculator, use [DUAL_CALC_Computer.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/DUAL_CALC_Computer.py). <br />
       We can specify the optimizer algorithm so that the we can decide what controls are kept constant and what controls are not. <br />
       <br />
  </details>
  <details>
  <summary>Combining Forces and Energies </summary><br />
        To sum the forces and energies in two different xyz files, we parse the xyz files, extract them and add them. We then write them back to the the rattled structure file. We then end up having forces and energies that are the sum of those in the two different xyz files. We can do this using the following script, [combiner_forces_energy.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/combiner_forces_energy.py).
  </details>  

  <details>
    <summary>Removing Free Energy from MACE-evaluated Dataset </summary> <br/>
    To remove free energy from the MACE-evaluted Dataset not to get confused later, we can use the following script [free_energy_remover.py](https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/free_energy_remover.py)
  </details>
   
     

 
   
   
     
           
                                        
