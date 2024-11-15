# msfr-multiregion-dev
OpenFOAM-8 solvers for the Molten Salt Fast Reactor with an additional layer of hastelloy to mimic the outer reflector region.

# OpenFOAM Required Version
OpenFOAM v8 (CFD Direct version)

## Installation
To install the solver, navigate to `msfr-multiregion-dev/OpenFOAM-nrg/src` and run `./Allwmake` (OpenFOAM environment must be active).

## Test
Run `msfrSimpleStructureFoam`. If you encounter the error 'cannot find file "......./src/system/controlDict', the installation was successful.

## General Comments
- The mesh included in the test case is the 2D-EVOL geometry (attached in the email).
- We don't simulate the fertile blanket or the inner neutron absorber.
- We simulate the primary fluid loop and the reflector layer outside.

**Note:** This solver is not a `chtMultiRegionSolver` (this is something we plan to do). Fluid and solid are treated as two materials at the same phase separated by an interface. A brief description of the solver can be found in the attached paper (Section 3).

## Running the Test Cases
Navigate to `msfr-multiregion-dev/OpenFOAM-nrg/tutorials`. There are two types of solvers:
- `SIMPLE` - steady state solver
- `PIMPLE` - transient solver

In both cases, run the `./Allrun` script in the OF case folder. To clean the case, run `./Allclean`.

You may need to give yourself permission to run: `chmod 777 Allrun ; chmod 777 Allclean`.

### Structure of the Allrun Script
1. `blockMesh` - build the 2D mesh
2. `topoSet` - create the pump region
3. `setFields` - set the nuclear parameters as internalFields (see 0 folder)
4. `decomposePar` - perform the case decomposition to run in parallel
5. `renumberMesh` - renumber the mesh cells to optimize the parallel run
6. `msfrSimpleStructureFoam` - run the solver
7. `reconstructPar` - reconstruct the case

**Note:** The final time step in `controlDict` has been set to 5 (SIMPLE) and 1 (PIMPLE) just to check everything works.

## Comments About the Solver (After Run)
All nuclear parameters are treated as fields initialized using `setFields` (see 0 folder after running `./Allrun`). `dynamiccode` allows setting `semiImplicitSources` for pump (momentum) and heat exchanger (heat transfer).

- In the `constant` folder, the `fvOptions` file allows setting momentum and heat sources at every time-step during runtime. Here in `fvOptions`, you can set the pump velocity and the heat transfer coefficient at the HX.
- `constant/nuclearProperties` - determine the number of delayed neutrons and decay heat precursors, the energy groups, the properties for the Doppler effect for transport.

The `postProcessing` folder gives you averaged quantities in time:
- Pump: mass flow rate at the outlet of the pump and pump average velocity
- Core: core average power
- HX: power at the HX
- Reactor parameters: k_effective, reactivity and Q

In the `System` folder:
- Specify the core nominal power in `controlDict`
- Specify the values of the nuclear parameters in `setFieldsDict` in both fluid and structure
- `topoSet` creates the pump region

## About the Transient (PIMPLE) Solver (What Changes from the SIMPLE Case)
Different probes are stored in `postProcessing` (after solver run). For the PIMPLE case, you set up the transient you want to simulate in `fvOptions`. Currently, we can simulate:
- ULOHS (unprotected loss of heat sink) by setting, for example, an exponential decay in time (x is the time variable)
- OVERSPEED of the pump, by setting, for example, an exponential increase in the velocity of the pump
- ULOFF (unprotected loss of flow), by setting, for example, an exponential decrease in the velocity of the pump

## Data Visualization
In the `msfr-multiregion-dev` folders, you can find some Jupyter notebooks (*.ipynb) which we use to visualize the data. You can open them using Jupyter notebook (in Ubuntu: `pip install jupyterlab ; pip install notebook`).

- `msfr_transients.ipynb` - visualize the transients you are going to impose on the transient case
- `RoutinesForPostProcess/plot_postProcess_*.ipynb` - plot the residuals (you can modify this to plot also the rest of the contents of `postProcessing`)

## Python Routines for Parametric Solutions
In the `RoutinesForParametricSolution` folder, you can find some Python routines to run the cases in a parametric way.
The Singular Value Decomposition is used to store the results (only some variables) in a compressed way to save storage space.
The `run.py` script includes the RITUAL for the Polimi Cluster, contact stefano.riva@polimi.it for more information.
