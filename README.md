# openFuelCell2

[openFuelCell2](https://github.com/openFuelCell2/openFuelCell2.github.io) is a computational fluid dynamics (CFD) toolbox for simulating electrochemical devices such as fuel cells and electrolysis. The solver is based on the open-source library, OpenFOAM®.

## About the code

The source code was developed from a previous open-source repository called [openFuelCell](http://openfuelcell.sourceforge.net/). It was also inspired by the standard solver "reactingTwoPhaseEulerFoam" in OpenFOAM®. The solver can consider coupled problems with multi-region and multi-physics issues, including single and two phase flows, multiple species components, charge transfer, and electrochemical reactions in different regions. More applications and solvers will be available in the future.

## How to use

The code is compiled with the OpenFOAM libraries, either [ORG](https://openfoam.org/) or [COM](https://www.openfoam.com/) versions. The default branch is compatable with the COM version, while the other branches are also provided for different OpenFOAM environments. The available environments will include: OpenFOAM-v2012, OpenFOAM-v2106, OpenFOAM-v2306, OpenFOAM-v6, OpenFOAM-v8. Note: the main branch is only compatible with OpenFOAM-v2306, while the others are under preparation.

```bash
# Download the source code
# Setup the corresponding openfoam environment
# Switch to the corresponding branch

# Change dictionary to the repository
cd openFuelCell2/src

# Compile the source code with
./Allwmake

# Or compile in parallel
./Allwmake -j n
```

You can also clear the libraries and executable files with

```bash

cd openFuelCell2/src

./Allwclean

```

## Computational domains

---

Take the cross-section of a fuel cell as an example. The computational domain gives,

<div align="center">
  <img src="images/computationDomain.jpg" height="70%" width="70%">
</div>

In a PEM fuel cell or other types, there are several domains/regions: air, fuel, electrolyte, and interconnect. This can be found from the repository [openFuelCell](http://openfuelcell.sourceforge.net/). However, additional domains/regions, e.g. phiEA, phiEC, and phiI are also necessary to account for electron/ion and dissolved water transfer.

To consider the coupling transfer problems in a PEM fuel cell, a global region, also called as parent mesh, and several local regions, also called as child meshes, are used. In the global region, only the energy equation is solved. In the local regions, corresponding partial differential equations will be discretized and solved. During the simulation, material properties, e.g. density, thermal conductivity, etc., are mapped from local regions to the global region, while the temperature field is mapped from the global region to the local regions.

The local regions can be classified as three different types, namely fluid, solid, and electric regions. See the [code](src/libSrc/fuelCellSystems/regions).

- Fluid region:

  This region represents the space where fluid flows by. In a fuel cell or electrolyzer, it consists with gas channels and/or porous regions. In this region, the following processes are addressed:

  - Fluid flow (single/two phase)
  - Species transfer
  - Electrochemical reaction
  - Heat and mass transfer

  For example, in a fuel cell, the following parts apply to this region,

  - Air flow paths + porous electrodes
  - Fuel flow paths + porous electrodes
  - Cooling channels

- Solid region:

  This represents the solid components in fuel cell or electrolyzer. In the current solver, no equations will be solved here. However, stress analysis during assembly and thermal effects may be implemented in future applications.

  For example. in a fuel cell, the following components apply to this region,

  - Electrolyte/membrane
  - Interconnect/Bipolar-plate
  - Endplate

- Electric region:

  This region accounts for the electric-conductive components. It is designed to consider electron/ion transfer specifically. However, it is found that the proton transfer region is the same as the region where dissolved water transfer takes place. Therefore, a switcher is enabled in the code to turn on/off the dissolved water transfer model. The following equations will be solved,

  - Potential equations (Poisson equations)
  - Dissolved water transfer equations (diffusion and electro-osmotic-drag)

  For example, in a fuel cell, the following components belong to this region:

  - Bipolar-plate, GDL, MPL, CL -> electron transfer regions
  - Catalyst coated membrane, CCM, -> proton transfer and dissolved water transfer region

- Global region:

  The heat source/sinks in local regions will be mapped to this region. And the obtained temperature is mapped back to the local regions. The heat source/sink include:

  - Joule heat from the electron/proton regions.
  - Condensation/evaporation in the fluid regions.
  - Electrochemical reactions in the fluid regions.

## Recent updates

---

- [Oct. 2020] A new branch for openFOAM-ESI
  > The majority part of this update was conducted by Mr. Steffen Hess. The code is able to compile in the OpenFOAM-ESI environment.
- [Nov. 2020] The new branch for openFOAM-ESI
  > Some bugs were found and fixed:
    1. The compiling sequence.
    2. The locations of gravity fields, g. The files "g" move to constant/.
    3. The functionObjects library is missing.
    4. Remove some warnings: apply new functions in OpenFOAM-ESI.
- [Nov. 2021] The new branch for openFOAM-2106
  > Some bugs were found and fixed:
    1. The method **heatTransfer(T, cellListIO)** in class "TwoResistanceHeatTransferPhaseSystem" is fixed.
  > Dictionary structure is rearranged:
    1. src: source files
    2. appSrc: application source files
    3. libSrc: libraries source files
    4. run: test cases
- [Jun. 2022] Clean the source code for releasing
  > Bugs are found and fixed:
    1. The previous solver might predict results with singularities in phiI region. This is fixed.
  > The source code is updated:
    1. Change the phase name to none if single-phase flow is simulated. This makes the variable names change from _A.air_ to _A_.
    2. Moving the correction of diffusivity from MultiComponentPhaseModel to **diffusivityModelList**.
    3. Moving the definition of phase properties from phaseProperties to regionProperties.
    4. Avoiding redundant output of diffusivity coefficients.
    5. Applying a different method to Update the value of phi in phiI region --> use setReference of phiEqn. This seems to make the solution more stable.
    6. Making the "porousZone" flexible to read. If the file doesn't exist, no porous zones are applied.
    7. Update the test cases: rewrite the scripts.
    8. Change the header of each file. openFuelCell is included.
- [Dec. 2022] Update the repository
    1. Fixed a bug in diffusivityList
    2. Update the tutorial
- [Sep. 2023] Update the repository for public release
    1. A new branch is included for OpenFOAM-v2306
    2. The prescribed mean current density and voltage are now defined as a function of time. (Assailable functions can be found in OpenFOAM/primitives/functions/Function1).
    3. Include the radiation model in solid region.
    4. Copy thermoTools to the repo. (need to remove this in next update.)
    5. In test cases, when it comes to two-phase flow, a steadyState scheme is now used, specially for ddt term of species transfer.
    6. Update the preprocessing script for an easier usage.

## Related publications

- Journal

  - Zhang, Shidong, Steffen Hess, Holger Marschall, Uwe Reimer, Steven Beale, and Werner Lehnert. "openFuelCell2: A New Computational Tool for Fuel Cells, Electrolyzers, and other Electrochemical Devices and Processes." Computer Physics Communications, Forthcoming (2023).

  - Zhang, Shidong, Shangzhe Yu, Roland Peters, Steven B. Beale, Holger Marschall, Felix Kunz, and Rüdiger-A. Eichel. "A new procedure for rapid convergence in numerical performance calculations of electrochemical cells." Electrochimica Acta (2023): 143275.

  - Yu, Shangzhe, Shidong Zhang, Dominik Schäfer, Roland Peters, Felix Kunz, and Rüdiger-A. Eichel. "Numerical Modeling and Simulation of the Solid Oxide Cell Stacks and Metal Interconnect Oxidation with OpenFOAM." Energies 16, no. 9 (2023): 3827.

  - Zhang, Shidong, Steven B. Beale, Uwe Reimer, Martine Andersson, and Werner Lehnert. "Polymer electrolyte fuel cell modeling-A comparison of two models with different levels of complexity." International Journal of Hydrogen Energy 45, no. 38 (2020): 19761-19777.

  - Zhang, Shidong. "Low-Temperature Polymer Electrolyte Fuel Cells." In Electrochemical Cell Calculations with OpenFOAM, pp. 59-85. Springer, Cham, 2022.

- Conference

  - Hess, Steffen, Shidong Zhang, Thomas Kadyk, Werner Lehnert, Michael Eikerling, and Steven B. Beale. "Numerical Two-Phase Simulations of Alkaline Water Electrolyzers." ECS Transactions 112, no. 4 (2023): 419.

  - Zhang, Shidong, Kai Wang, Shangzhe Yu, Nicolas Kruse, Roland Peters, Felix Kunz, and Rudiger-A. Eichel. "Multiscale and Multiphysical Numerical Simulations of Solid Oxide Cell (SOC)." ECS Transactions 111, no. 6 (2023): 937.

  - Zhang, Shidong, Steven B. Beale, Yan Shi, Holger Janßen, Uwe Reimer, and Werner Lehnert. "Development of an Open-Source Solver for Polymer Electrolyte Fuel Cells." ECS Transactions 98, no. 9 (2020): 317.

- Thesis

  - Zhang, Shidong. Modeling and Simulation of Polymer Electrolyte Fuel Cells. No. FZJ-2020-02318. Elektrochemische Verfahrenstechnik, 2020.

## Developers

---

The code is firstly developed by [Shidong Zhang](s.zhang@fz-juelich.de) for the PhD thesis, supervised by Prof. [Werner Lehnert](w.lehnert@fz-juelich.de) and Prof. [Steven Beale](s.beale@fz-juelich.de). The detailed model description and simulation results can be found in the thesis, `Modeling and Simulation of Polymer Electrolyte Fuel Cells` by FZJ. The following individuals also contribute to the optimization of the code,

- Steffen Hess (s.hess@fz-juelich.de), Forschungszentrum Juelich, IEK-14

- Prof. Steven B. Beale (s.beale@fz-juelich.de), Forschungszentrum Juelich, IEK-13

To be continued...
