# __PEM fuel cell application__

This is a test case for proton exchange membrane / polymer electrolyte fuel cell.

Chemical reaction:

- Anode side:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{H}_{2}\to\textrm{2H}^{+}+\textrm{2e}^{-}" title="\Large \textrm{H}_{2}\to\textrm{2H}^{+}+\textrm{e}^{-}" />
- Cathode side:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{0.5O}_{2}+\textrm{2H}^{+}+\textrm{2e}^{-}\to\textrm{H}_{2}\textrm{O}" title="\Large \textrm{0.5O}_{2}+\textrm{2H}^{+}+\textrm{2e}^{-}\to\textrm{H}_{2}\textrm{O}" />
- Overall Reaction:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{H}_{2}+\textrm{0.5O}_{2}\to\textrm{H}_{2}\textrm{O}" title="\Large \textrm{H}_{2}+\textrm{0.5O}_{2}\to\textrm{H}_{2}\textrm{O}" />

___The water may exist in vapor, liquid, or both phases in PEM fuel cells, depending on operating temperature___

## Operating conditions

    ```
    Temperature:            353 K, 80 oC
    Pressure:               1 bar/1 bar (anode/cathode)
    gases:                  (H2 & H2O)/(air& H2O) (anode/cathode)
    Relative humidity       90%/50%
    Active area             0.8 cm2

    Mean current density    8000 A/m3
    ```

To run the case:

- In serial

    ```bash

    # Generate the computational meshes
    make mesh

    # Run in serial
    make srun

    ```

- In parallel

    ```bash

    # Generate the computational meshes
    make mesh

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
    #or
    gnuplot ResidualPlot

    ```

  - Simulation results

    ```bash
    make view
    paraview

    ```

- In parallel

  - Residual plot

    ```bash
    make plot
    or
    gnuplot ResidualPlot
    ```

  - Simulation results

    ```bash
    make reconstruct
    make view
    paraview
    ```
