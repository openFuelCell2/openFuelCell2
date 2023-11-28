# __Conjugate heat transfer application__

This is a test case for a heat exchanger. hot <--> solid <--> cold

## Operating conditions:

    ```none
    Temperature(K):         hot(973.15), cold (293.15)
    Pressure:               1 bar/1 bar (hot/cold)
    gases:                  air/H2 (hot/cold)
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
