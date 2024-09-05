# Python Scripts Used
Using machine learning interatomic potentials to study the photodeterioration of materials

**Table of Content**
  - Dataset Creation
  - Conversion of trajectory files to xyz files
  - Calculating forces 
  - Calculating energies
  - Extracting forces
  - Extracting energies
  - Training MACE Models
  - Evaluating MACE Models
  - Computing observables like Potential Energy, Kinetic Energy, Radial Distribution Function and Mean-squared Displacement
  - Combining forces and energies
  - Removing free energy from MACE-evaluted files

  **Dataset Creation** <br />
  <br />
        Using the LiFDatasetCreator.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/LiFDatasetCreator.py), we can create rattled structures of materials. <br />
           => to change the number of structures to be created, assign the "N" argument in the generate_random_strucutres_ratlle() function to the number of the rattled strucutres that you want to generate.

  **Conversion of Trajectory Files to XYZ-files**<br />
  <br />
         To convert files that have .traj extension to xyz file, use the convert_traj.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/convert_traj.py)

  **Calculating Forces and Energies** <br />
  <br />
    After creating the rattled structures, to calculate the forces and the energies, use energy_force_calculator.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/energy_force_calculator.py). This file takes one xyz file and creates three different xyz output file: one calculated with single MACE calculator, one calculated with double MACE calculators, one LAMMPS calculator.
      
     
 **Energy Extractor**<br />
 <br /> 
        This file takes an xyz file as an input and extracts the energy in the file and writes it to a csv file. The file is energy_extractor.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/energy_extractor.py)

 **Force Extractor**<br />
 <br />
        => This file takes an xyz file as an input and extracts the force in the file and writes it to a csv file. The file is force_extractor.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/force_extractor.py)
    
 **Training MACE Models**<br />
 <br />
         => To train MACE models, use the train_mace_model.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/train_mace_model.py). To learn more about the parameters, visit the following github page (https://github.com/ACEsuit/mace)

  **Evaluating MACE Models**<br />
  <br />
          To evalute the trained MACE models, use evaluate_mace.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/evaluate_mace.py). You can find the trained model inside the "results" directory specified in the training script.

   **Computing Obervables (Potential Energy, Kinetic Energy, Radial Distribution FUnction, Mean-squared Displacement)**<br />
       => To compute with LAMMPS, use lammps_LiFCalculator.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/lammps_LiFCalculator.py). <br /><br />
       => To compute observables with one MACE model using MACE calculator, use MACE_Computer.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/MACE_COMPUTER.py). <br /><br />
       => To compute observables with two MACE models using MACE calculator, use DUAL_CALC_Computer.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/DUAL_CALC_Computer.py) <br />
       We can specify the optimizer algorithm so that the we can decide what controls are kept constant and what controls are not. <br />
       <br />
  **Combining Forces and Energies**<br />
  <br />
        To sum the forces and energies in two different xyz files, we parse the xyz files, extract them and add them. We then write them back to the the rattled structure file. We then end up having forces and energies that are the sum of those in the two different xyz files. We can do this using the following script, combiner_forces_energy.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/combiner_forces_energy.py).
     
   
     

 
   
   
     
           
                                        
