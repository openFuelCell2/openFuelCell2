# **fuelCell0Foam -- clean version**

__This is uploaded to have a repository as backup__

This is going to be cleaned.

The source code was developed from an open-source repository [openFuelCell](http://openfuelcell.sourceforge.net/) and the standard solver "reactingTwoPhaseEulerFoam" in OpenFOAM. It can be used to consider coupling transport phenomena in electrochemical devices, e.g. fuel cells, electrolyzers, etc. More applications will be available in the future.

At this point, the source code needs to be cleaned up and necessary comments should be added. This code will hopefully be made open-source for fuel cell/electrolyzer communities. I hope this can be done as soon as possible.

There will be additional people who are interested in helping cleaning the code, I will keeping adding them in the list.

Names:

- Prof. Steven Beale (s.beale@fz-juelich.de), Forschungszentrum Juelich, IEK-14

## How the code works:

### To compile the code:

The code should be based on the foundation version of [OpenFOAM](https://openfoam.org/), version-6. The code will be updated with the latest OpenFOAM. Anyone who has installed OpenFOAM-v6 on his/her linux machine is able to compile the code and runs the test case.

```
git clone XXX fuelCell0Foam

cd fuelCell0Foam

./Allwmake
```

You can also clear the libraries and executable files with

```
./Allwclean
```
### The cross-section of a fuel cell (as an example):

[Test cases](tutorial/) have been prepared for users in tutorial. The computational domain ![Computation domain](images/computationDomain.jpg)

In a PEM fuel cell, there are several domains/regions: air, fuel, electrolyte, and interconnect. This can be found from the repository [openFuelCell](http://openfuelcell.sourceforge.net/). However, additional domains/regions, e.g. phiEA, phiEC, and phiI are also necessary for electron/proton and dissolved water transfer.

To consider the coupling transfer problems in a PEM fuel cell, a master region, also called as parent mesh, and several sub-regions, also called as child meshes, are used. In the master region, enthalpy equation is solved. In the sub-regions, corresponding partial differential equations, PDEs, will be solved.

The sub-regions can be classified as three different types, fluid, solid, and electron/proton transfer regions. See the [code](fuelCellSystems/fuelCellSystem/regions).

- Fluid region:

  This region represents the space where fluid flows by. In fuel cells, it may consist with gas channels and/or porous regions. In these regions, several transfer phenomena should be obtained numerically:

  - Fluid flow
  - Species transfer
  - Two-phase flow
  - Electrochemical reaction
  - Heat transfer

  For example, in a fuel cell, the following regions belong to this type:

  - Air
  - Fuel
  - Cooling

- Solid region:

  This represents the solid space. In the present applications, no eqautions will be solved here. However, in future applications, stress analysis during assembly and thermal effects may be implemented.

  In a fuel cell, the following regions belong to this type:

  - Electrolyte/membrane
  - Interconnect/Bipolar-plate
  - Endplate

- Electron/proton region:

  This considers the electric-conductive regions. This region was designed to consider electron/proton transfer specifically. However, it was found that the proton transfer region is the same as the region where dissolved water transfer takes place. Therefore, a switcher is set in the code to turn on/off the dissolved water transfer model. The following eqautions will be solved:

  - Potential equations (Possion equations)
  - Dissolved water transfer equations (diffusion and electro-osmotic-drag)

  In a fuel cell, the following regions belong to this type:

  - Bipolar-plate, GDL, MPL, CL -> electron transfer regions
  - Catalyst coated membrane, CCM, -> proton transfer and dissolved water transfer region

- Master region:

  The heat source/sinks in sub-regions will be mapped to this region. And the obtained temperature is mapped back to sub-regions. The heat source/sink include:

  - Joule heat from the electron/proton regions.
  - Phase change in the fluid regions.
  - Electrochemical reaction in the fluid regions.

To be continued...

### Code structures:

- 1st level:

```
  fuelCellFoam.C  [Main code]
  EEqns.H         [Enthalpy equation]
  createFields.H  [Create necessary fields]
  fuelCellSystems [All of the necessary models]
  include         [Courant, delta time, etc.]
  tools           [Tools for topoSet and mesh decomposition]
  tutorial        [Test cases]
```

- 2nd level:

  - fuelCellSystems

```
    derivedFvPatchFields          [Two phase models boundary conditions]
    functionObjects               [Two phase models functions]
    interfacialCompositionModels  [Two phase models, interfacial composition, mass transfer, saturation, and surface tension]
    interfacialModels             [Two phase models, interfacial momentum transfer terms]
    multiPhaseCompressibleTurbulenceModels  [Turbulent model for two phase flow, not necessary at this point.]
    regionCourantNo               [Courant number]
    fuelCellSystem                [Models for fuel cells or other applications]
```

  - tools

```
    decomposeParID              [Generate cell IDs for manual decomposition]
    masterRegionToCell          [Create cell sets from the regions in master region]
    topoSet                     [Recompile topoSet to include masterRegionToCell]
```

- 3rd level (fuelCellSystem):

```
      activationOverpotentialModels     [Activational overpotential models, e.g. Butler-Volmer, Tafel]
      diameterModels                    [Droplets/Bubbles diamter]
      dissolvedModel                    [Dissolved water transfer in membrane]
      nernstModel                       [Nernst potential]
      phasePair                         [Types of paired objects]
      PhaseSystems                      [Templates for different problems in multi-phase flows]
      porosityModel                     [Porosity models in porous regions]
      regions                           [Different types of regions]
      solvers                           [Selectable solvers, e.g. single phase, two phase, or drift-flux]
      BlenderInterfacialModels          [Blender functions for interfacial momentum transfer terms]
      diffusivityModels                 [Species diffusion models, Fick's law]
      phaseModel                        [Multi-models for each phase]
      phaseSystem                       [Phase systems for fluid flow]
      sigmaModels                       [Electric conductivity]
      populationBalanceModel            [Model for bubble columns]
      reactionThermo                    [Thermo models for reactions]
```

To be continued...

### Case structure:

```
  0/                                    [Initial and boundary conditions]
    air/                                  [Air region]
      alpha.air                             [Fields of gas phase saturation]
      alpha.water                           [Fields of liquid phase saturation]
      T.air                                 [Fields of gas phase temperature]
      T.water                               [Fields of liquid phase temperature]
      H2O.air                               [Fields of water vapor mass fraction]
      N2.air                                [Fields of nitrogen mass fraction]
      O2.air                                [Fields of oxygen mass fraction]
      p                                     [Fields of pressure]
      p_rgh                                 [Fields of pressure - rho*g*h]
      U.air                                 [Fields of gas phase velocity]
      U.water                               [Fields of liquid phase velocity]
    fuel/                                 [Fuel region]
      alpha.fuel                            [Fields of gas phase saturation, always 1 in single phase flow]
      T.fuel                                [Fields of temperature]
      H2.fuel                               [Fields of hydrogen mass fraction]
      H2O.fuel                              [Fields of water vapor mass fraction]
      p                                     [Fields of pressure
      p_rgh                                 [Fields of pressure - rho*g*h]
      U.fuel                                [Fields of velocity]
    electrolyte/                          [Electrolyte region]
      p                                     [Fields of pressure, constant in solid region]
      T                                     [Fields of temperature]
    interconnect/                         [Interconnect region]
      p                                     [Fields of pressure, constant in solid region]
      T                                     [Fields of temperature]
    phiEA/                                [PhiEA region, potential field at anode side]
      phi                                   [Fields of potential]
    phiEC/                                [phiEC region, potential field at cathode side]
      phi                                   [Fields of potential]
    phiI/                                 [phiI region, potential field at CCM]
      lambda                                [Fields of water content]
      phi                                   [Fields of potential]
  config/                               [Information for cell sets]
  Makefile                              [Make file]
  parallel.csh                          [Script for parallel pre-processing]
  pre                                   [Script for pre-processing]
  pResidualPlot                         [Residual plot for parallel runs]
  sResidualPlot                         [Residual plot for serial runs]
  constant                              [Model and parameter selections, polyMesh]
    cellProperties                        [General settings of the application]
    regionProperties                      [Types of regions]
    polyMesh/                             [Finite volume mesh information]
    air/                                  [Air region]
      g                                     [Gravity]
      phaseProperties                       [Properties of each phase]
      regionProperties                      [Properties of air region]
      thermophysicalProperties.air          [Thermo information of the gas phase]
      thermophysicalProperties.water        [Thermo information of the liquid phase]
      combustionProperties.air              [Electrochemical reaction settings]
      porousZones                           [Porosity regions]
      turbulenceProperties.air              [Turbulent models for the gas phase. Laminar in use]
      turbulenceProperties.water            [Turbulent models for the liquid phase. Laminar in use]
      polyMesh/                             [Finite volume mesh information]
    fuel/                                 [Fuel region]
      g                                     [Gravity]
      phaseProperties                       [Properties of each phase]
      regionProperties                      [Properties of fuel region]
      thermophysicalProperties.fuel         [Thermo information of the gas phase]
      combustionProperties.fuel             [Electrochemical reaction settings]
      porousZones                           [Porosity regions]
      turbulenceProperties.fuel             [Turbulent models for the gas phase. Laminar in use]
      polyMesh/                             [Finite volume mesh information]
    electrolyte/                          [Electrolyte region]
      regionProperties                      [Properties of electrolyte region]
      thermophysicalProperties              [Thermo information]
      polyMesh/                             [Finite volume mesh information]
    interconnect/                         [Electrolyte region]
      regionProperties                      [Properties of interconnect region]
      thermophysicalProperties              [Thermo information]
      polyMesh/                             [Finite volume mesh information]
    phiEA/                                [PhiEA region]
      regionProperties                      [Properties of phiEA region]
      polyMesh/                             [Finite volume mesh information]
    phiEC/                                [PhiEC region]
      regionProperties                      [Properties of phiEC region]
      polyMesh/                             [Finite volume mesh information]
    phiI/                                 [PhiI region]
      regionProperties                      [Properties of phiI region]
      polyMesh/                             [Finite volume mesh information]
  system/                               [Finite volume method simulation settings]
    controlDict                           [Simulation control]
    blockMeshDict                         [Mesh generation]
    decomposeParDict                      [Mesh decomposition]
    fvSchemes                             [Schemes in temporal and spatial discretization]
    fvSolution                            [Matrix solution and solution control]
    air/                                  [Finite volume method settings to air region]
      decomposeParDict                      [Mesh decomposition]
      fvSchemes                             [Schemes in temporal and spatial discretization]
      fvSolution                            [Matrix solution and solution control]
    fuel/                                 [Finite volume method settings to fuel region]
      decomposeParDict                      [Mesh decomposition]
      fvSchemes                             [Schemes in temporal and spatial discretization]
      fvSolution                            [Matrix solution and solution control]
    electrolyte/                          Finite volume method settings to electrolyte region]
      decomposeParDict                      [Mesh decomposition]
    interconnect/                         Finite volume method settings to interconnect region]
      decomposeParDict                      [Mesh decomposition]
    phiEA/                                [Finite volume method settings to phiEA region]
      decomposeParDict                      [Mesh decomposition]
      fvSchemes                             [Schemes in temporal and spatial discretization]
      fvSolution                            [Matrix solution and solution control]
    phiEC/                                [Finite volume method settings to phiEC region]
      decomposeParDict                      [Mesh decomposition]
      fvSchemes                             [Schemes in temporal and spatial discretization]
      fvSolution                            [Matrix solution and solution control]
    phiI/                                 [Finite volume method settings to phiI region]
      decomposeParDict                      [Mesh decomposition]
      fvSchemes                             [Schemes in temporal and spatial discretization]
      fvSolution                            [Matrix solution and solution control]
```

To be continued...
