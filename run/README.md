# openFuelCell2

## What tutorial cases are provided?

For each application, which is able to be simulated by the code, a tutorial case is provided.

- conjugate heat transfer: CHT
- Low temperature PEM electrolyzer: PEMEC
- Low temperature PEM fuel cell: PEMFC
- High temperature PEM electrolyzer: HTPEMEC
- High temperature PEM fuel cell: HTPEMFC
- Hydrogen pump: hydrogenPump
- Solid oxide fuel cell: SOFC
- Solid oxide electrolyzer: SOEC

Each tutorial:

- Common OpenFoam cases structure:
  - 0.orig : Backup of 0
  - constant
  - system
  - _Makefile_: To generate mesh, run the solver, view the solution...
  - _Allclean_: To clean up the tutorial
  - _parallel.sh_: To generate processors for parallel simulation
  - _preprocessing.sh_: To generate mesh for serial simulation

## How to run the tutorial cases?

- In serial

    ```bash
    # Generate the meshes with 'blockMesh'
    make mesh
    # or Generate the meshes with 'salome'
    # make salomeMesh
    # Run in serial
    make srun
    ```

- In parallel

    ```bash
    # Generate the meshes with 'blockMesh'
    make mesh
    # or Generate the meshes with 'salome'
    # make salomeMesh
    # Edit values of 'nx' and 'ny' in constant/cellProperties
    export NPROCS=nx*ny # (value of 'nx' times 'ny')
    # Generate the cellID
    # For manual decomposition
    make decompose
    # Decomposition for multiple regions
    make parallel
    # Run in parallel
    make run  #( the value of 'nx' times 'ny' needs to be changed in 'Makefile')
    # or
    # mpirun -np nx*ny fuelCell0Foam -parallel -fileHandler collated | tee log.run
    ```

To view the result:

- In serial

  - Residual plot

    ```bash
    make plot
    # or gnuplot ResidualPlot
    ```

  - Simulation results

    ```bash
    # Convert openfoam result to vtk
    make view
    # execute paraview to view
    paraview
    ```

- In parallel

  - Residual plot

    ```bash
    # plot the residuals with gnuplot
    make plot
    # or gnuplot ResidualPlot
    ```

  - Simulation results

    ```bash
    # Reconstruct the mesh
    make reconstruct
    # Convert the results to vtk
    make view
    # Execute paraview to view
    paraview
    ```

## Where can I find more information?

We are happy to provide more information. If you have problems, please send us Emails:

- Contact persons:
  - Shidong Zhang (s.zhang@fz-juelich.de)
    - Main Developer
  - Steffen Hess (s.hess@fz-juelich.de)
    - Developer
  - Steven B. Beale (s.beale@fz-juelich.de)
    - Supervisor