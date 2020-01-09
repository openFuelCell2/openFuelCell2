__fuelCell0Foam -- clean version__

This is a test case for proton exchange membrane / polymer electrolyte fuel cell.

Operating conditions:

'''
Temperature:            353 K, 80 oC
Pressure:               1 bar/1 bar (anode/cathode)
Stoichiometric factor:  2/2 (anode/cathode)
gases:                  H2/air (anode/cathode)
Relative humidity       2/2
Active area             1 cm2

Mean current density    8000 A/m3
'''

To run the case:

'''
fluent3DToFoam  coarse.msh
./pre

fuelCell0Foam > log
'''

To run the case in parallel:

'''
fluent3DToFoam coarse.msh
./pre

Edit values of 'nx' and 'ny' in constant/cellProperties

decomposeParID -coordinate '((0 -0.001 0) (0.05 0.001 0))' -allRegions
decomposeParID -coordinate '((0 -0.001 0) (0.05 0.001 0))'

./parallel-collated.csh
mpirun -np XX fuelCell0Foam -parallel -fileHandler collated > log
'''
