# __HT-PEM fuel cell application__

This is a test case for high temperature proton exchange membrane / polymer electrolyte fuel cell.

Chemical reaction:

- Anode side:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{H}_{2}\to\textrm{2H}^{+}+\textrm{2e}^{-}" title="\Large \textrm{H}_{2}\to\textrm{2H}^{+}+\textrm{e}^{-}" />
- Cathode side:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{0.5O}_{2}+\textrm{2H}^{+}+\textrm{2e}^{-}\to\textrm{H}_{2}\textrm{O}" title="\Large \textrm{0.5O}_{2}+\textrm{2H}^{+}+\textrm{2e}^{-}\to\textrm{H}_{2}\textrm{O}" />
- Overall Reaction:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\textrm{H}_{2}+\textrm{0.5O}_{2}\to\textrm{H}_{2}\textrm{O}" title="\Large \textrm{H}_{2}+\textrm{0.5O}_{2}\to\textrm{H}_{2}\textrm{O}" />

# Operating conditions:

```
Temperature:            433 K, 160 oC
Pressure:               1 bar/1 bar (anode/cathode)
Stoichiometric factor:  2/2 (anode/cathode)
gases:                  (H2 & H2O)/(air& H2O) (anode/cathode)
Active area             0.8 cm2

Mean current density    8000 A/m3
```

To run the case:

- In serial

```
make mesh

make srun
```

- In parallel


```
make mesh

Edit values of 'nx' and 'ny' in constant/cellProperties

export NPROCS=nx*ny (value of 'nx' times 'ny')

decomposeParID -coordinate '((0 -0.001 0) (0.04 0.001 0))' -allRegions
decomposeParID -coordinate '((0 -0.001 0) (0.04 0.001 0))'

make parallel

make prun ( the value of 'nx' times 'ny' needs to be changed in 'Makefile')
or
mpirun -np nx*ny fuelCell0Foam -parallel -fileHandler collated | tee log.run
```

To view the result:

- In serial

  - Residual plot
```
make plot
or
gnuplot ResidualPlot
```
  - Simulation results
```
make view
paraview
```

- In parallel

  - Residual plot
```
make plot
or
gnuplot ResidualPlot
```
  - Simulation results
```
make reconstruct
make view
paraview
```
