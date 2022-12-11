/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
#----------------------------------------------------------------------#       #
# Solver      |   openFuelCell                                         #       #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #       #
#----------------------------------------------------------------------#       #
# Source code |   https://jugit.fz-juelich.de/s.zhang/fuelcellfoam     #       #
# Update from |   15.06.2022                                           #       #
#----------------------------------------------------------------------#       #
                                                                               #
################################################################################

# Plot probe values
# plot the times needed per iteration
gnuplot<<\EOF
    set terminal pngcairo size 608,813  enhanced font 'Verdana,10'
    set output "airProbePressure.eps"
    set xlabel "iteration"
    set ylabel "pressure in Pa"
    set grid
    set samples 1000
    set key right top
    plot \
    "postProcessing/probePressureAir/air/0/p_rgh" u 1:2 with lines t "pressureProbe1", \
    "postProcessing/probePressureAir/air/0/p_rgh" u 1:3 with lines t "pressureProbe2", \
EOF

# Plot probe values - alpha.water -
# plot the times needed per iteration
gnuplot<<\EOF
    set terminal pngcairo size 608,813  enhanced font 'Verdana,10'
    set output "airWaterProbe.eps"
    set xlabel "iteration"
    set ylabel "s in [-]"
    set grid
    set samples 1000
    set key right top
    plot \
    "postProcessing/probesAir/air/0/alpha.water" u 1:2 with lines t "H2OSaturationProbe1", \
EOF

# Plot probe values - j -
# plot the times needed per iteration
gnuplot<<\EOF
    set terminal pngcairo size 608,813  enhanced font 'Verdana,10'
    set output "jProbe.eps"
    set xlabel "iteration"
    set ylabel "current density in A/m^3"
    set grid
    set samples 1000
    set key right top
    plot \
    "postProcessing/probesAir/air/0/j" u 1:2 with lines t "jProbe", \
EOF


# Plot probe values - O2 -
# plot the times needed per iteration
gnuplot<<\EOF
    set terminal pngcairo size 608,813  enhanced font 'Verdana,10'
    set output "O2H2Probe.eps"
    set xlabel "iteration"
    set ylabel "mass fraction in [-]"
    set grid
    set samples 1000
    set key right top
    plot \
    "postProcessing/probesAir/air/0/O2.air" u 1:2 with lines t "O2CathodeProbe", \
    "postProcessing/probesFuel/fuel/0/H2" u 1:2 with lines t "H2AnodeProbe", \
EOF


# plot residuals

gnuplot<<\EOF
    set terminal pngcairo size 608,813  enhanced font 'Verdana,10'
    set output "residuals.eps"
    set multiplot layout 2, 1
    set xlabel "timeStep"
    set ylabel "Residuals"
    set format y "10^{%L}"
    set grid
    set samples 1000
    set key right top
    set logscale y
    set ytics nomirror
    set key outside top
    # set xrange [0:1]
    # set yrange [105:140]
    plot \
    "postProcessing/air/residualsAir/0/solverInfo.dat" u 1:3 with lines t "O2.airInit", \
    "postProcessing/air/residualsAir/0/solverInfo.dat" u 1:4 with lines t "O2.airFinal", \
    "postProcessing/air/residualsAir/0/solverInfo.dat" u 1:8 with lines t "prgh.airInit", \
    "postProcessing/air/residualsAir/0/solverInfo.dat" u 1:9 with lines t "prgh.airFinal", \
    "postProcessing/fuel/residualsFuel/0/solverInfo.dat" u 1:3 with lines t "prgh.fuelInit", \
    "postProcessing/fuel/residualsFuel/0/solverInfo.dat" u 1:4 with lines t "prgh.fuelFinal", \
EOF


