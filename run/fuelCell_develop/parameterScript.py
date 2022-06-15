#!/usr/bin/python3

import numpy as np

#current
current = 6000  # A/m^3

nProc = 4

# ================== Case geometry parameter ====================#

l = 0.04            # length of the channel in m
wRip = 0.0005       # width of the rip , due to symmetry half width
wCh = 0.001         # width of the flow channel
hCh = 0.001         # height of the flow channel


# ==================   Thermodynamical data ================= #

# Temperature
T = 353.15   # K
T_celcius = T - 273.15     # °C
Tfuel = 353.15

# Densities
rhoO2= 3.65383184 - 0.01147382*T + 0.0000119237*T**2
rhoAir = 3.3284522 - 0.01051861*T + 0.0000109931*T**2
rhoH2 = 0.22836866 - 0.00071277*T + 0.0000007368*T**2

# Surface tension
sigma = 0.07275 * (1 - 0.002*(T - 291))
theta = 30

pRefAir = 101325   # 1 bar
pRefFuel = 101325   # 1 bar
RAir = 287.1    # J/kg/K
RD = 461.5      # J/kg/K   Wasserdampf
phi = 1         # relative humidity Avg. in germany
wH2OA = 0       # mass fraction water in air


# =================== Calculate inlet values ====================#

relHumidity = 1

# Magnus-Formula over water - Temperatur -50 up to +100°C
p_sat = 611.2 * np.exp((17.62 * T_celcius) / (T_celcius + 243.5))

# Partial pressure water vapor (phi=1 bzw. 100% dann p = psat)
pi = phi * p_sat

# Mass of water vapor per m^3 dry air
mD = (phi * p_sat) / (RD * T)                     

rhoWD = pi / (RD * T)
wH2OA = rhoWD / (rhoWD + rhoAir)

wLeft = 1-wH2OA

O2Masspercentage = 0.23         # Mass fraction of oxygen in air
N2Masspercentage = 0.77         # Mass fraction of nitrogen in air

# Mass fraction oxygen and nitrogen
wO = wLeft * O2Masspercentage
wN = wLeft * N2Masspercentage

# Mass fraction fuel side (Hydrogen + water)
wH = 0.267
wH2OH = 0.733


Aact = (wCh + 2*wRip) * l     # active area
lambaA = 2                    # stoichiometric factor
lambaH = 2                    # stoichiometric factor
Ain = wCh * hCh               # inlet area
zA = 4
zH = 2                        # number of electrons transfer
F = 96485.332                 # Faraday constant  C/mol

MH2O = 18.01528 * 10**(-3)


# air inlet velocity
MO = 2 * 15.9994 * 10**(-3)         # molar mass of each reactant (oxygen or hydrogen)
MN = 2 * 14.00 * 10**(-3)         # molar mass of each reactant (oxygen or hydrogen)

XO = (wO / MO) / (((wO / MO) + (wN / MN) + (wH2OA / MH2O)))       # mole fraction of each reactant in the inlet gas

# Calculate the air inlet velocity for this stoichiometry 
uAInlet = (MO * Aact * current * lambaA) / (rhoO2 * XO * Ain * zA * F)


MH = 2 * 1.00777 * 10**(-3)                        # molar mass of each reactant (oxygen(air) or hydrogen)  kg/mol
XH =  (wH / MH) / ((wH /MH) + (wH2OH / MH2O))      # mole fraction of each reactant in the inlet gas


# Calculate the air inlet velocity for this stoichiometry 
uHInlet = (MH * Aact * current * lambaH) / (rhoH2 * XH * Ain * zH * F)

f = open("simulationParameter.in", 'w')
f.write("//- Geometrical data \n")
f.write("hCh     " + str(hCh  ) + ";\n")
f.write("wCh     " + str(wCh  ) + ";\n")
f.write("length     " + str(l  ) + ";\n")
f.write("wRip     " + str(wRip  ) + ";\n")
f.write("\n")
f.write("//- Operating data \n")
f.write("current     " + str(current  ) + ";\n")
f.write("theta     " + str(theta  ) + ";\n")
f.write("uAInlet     " + str(uAInlet  ) + ";\n")
f.write("uHInlet     " + str(uHInlet  ) + ";\n")
f.write("T     " + str(T  ) + ";\n")
f.write("wH2OA     " + str(wH2OA  ) + ";\n")
f.write("wH2OH     " + str(wH2OH  ) + ";\n")
f.write("wH     " + str(wH  ) + ";\n")
f.write("wO     " + str(wO  ) + ";\n")
f.write("wN     " + str(wN  ) + ";\n")
f.write("Tfuel     " + str(Tfuel  ) + ";\n")
f.write("pRefAir     " + str(pRefAir  ) + ";\n")
f.write("pRefFuel     " + str(pRefFuel  ) + ";\n")
f.write("sigma     " + str(sigma  ) + ";\n")
f.write("\n")
f.write("//- Simulation data \n")
f.write("nProc     " + str(nProc  ) + ";\n")
f.close()

    
