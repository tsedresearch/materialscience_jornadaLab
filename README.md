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
           
                                        
