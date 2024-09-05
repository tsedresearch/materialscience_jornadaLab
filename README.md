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

  Dataset Creation
    Using the LiFDatasetCreator.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/LiFDatasetCreator.py), 
    we can create rattled structures of materials. 
           => to change the number of structures to be created, assign the "N" argument in the generate_random_strucutres_ratlle() function to the number of the rattled strucutres that you want to generate.

  Conversion of trajectory files to xyz files
     => to convert files that have .traj extension to xyz file, use the convert_traj.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/convert_traj.py)

  ??? Calculating forces and energies
     => 
     
  Energy Extractor
    => This file takes an xyz file as an input and extracts the energy in the file and writes it to a csv file. The file is energy_extractor.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/energy_extractor.py)

  Force Extractor
    => => This file takes an xyz file as an input and extracts the force in the file and writes it to a csv file. The file is force_extractor.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/force_extractor.py)
    
  Training MACE Models
     => To train MACE models, use the train_mace_model.py (https://github.com/tsedresearch/materialscience_jornadaLab/blob/main/train_mace_model.py). To learn more about the parameters, visit the following github page (https://github.com/ACEsuit/mace)

 
   
   
     
           
                                        
