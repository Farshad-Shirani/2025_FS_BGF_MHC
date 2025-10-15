# 2025-FS-BGF-MHC

# Interspecific Competition can Set Range Limits for Adaptively Dispersing Species

Written by Farshad Shirani (f.shirani@gatech.edu), October 15, 2025.

This repository contains MATLAB codes used in the article:

"Coevolution of Species’ Borders: Interactions Between Interspecific Competition, Gene Flow, and Matching Habitat Choice", F. Shirani and B. G. Freeman, (2025).

All codes are written in MATLAB R2024b.

## Description
These MATLAB codes numerically compute solutions of systems of partial differential equations (PDEs) that simulate adaptive evolution of geographic ranges of two biological species in a one-dimensional geographic space. The species disperse both randomly and adaptively, through matching habitat choice.

## How to run the simulations:

- The code associated with each of the simulations in the abovementioned paper are provided in separate folders. The core of the code in different folders is the same, only some changes are applied to the `initializeSimulation` and `Main_Simulation` files, based on the specifics of each simulation. In particular, all codes included in the subfolder named `AuxiliaryFunctions` are the same for all simulations. These codes are associated with the numerical scheme that is used to solve the PDEs, and should not be changed, unless absolutely needed.

- Each folder contains a function named `initializeSimulation`. Model parameters, simulation parameters, and discretization parameters of the numerical scheme are first set by the user through this function. The model parameters are defined as global parameters. Therefore, the values of model parameters should only be changed through this function.

- Each folder contains a script named `Main_Simulation`, which is the main code that performs the simulations based on the parameters set in the function `initializeSimulation`. When a simulation is complete, the resulting solutions can be saved using the `save` command provided at the end of the script. The computed solutions are stored in a structure array named `populations`. The parameters are also stored in three different structures, named `modelParameters`, `simulationParameters`, and `discretizationParamaters`. The path given to the `save` command should include the subfolder name `Results`, so that all simulation results are organized in this subfolder. The `Results` folder is currently loaded with the simulation results that the authors have performed and used in the above-mentioned paper.

- The results can be plotted using the script `Plotting` available in each folder. The path of the results to be plotted should be given in the `load` command at the beginning of the script. The parameter `incrementSize` inside the script should be set to the desired incremental time steps that the resulting curves should be plotted.

## Which simulation (folder) is associate with which figure?

- The code in the folder `Range Dynamics s` is associated with Figures 1, 2, and 3 in the paper.

- The code in the folder `Environmental Fluctuations` is is associated with Figure 4 in the paper.
