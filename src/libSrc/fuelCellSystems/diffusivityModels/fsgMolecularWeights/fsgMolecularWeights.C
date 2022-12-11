/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fsgMolecularWeights.H"

namespace Foam
{
namespace diffusivityModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// [kg/kmol]

const fsgMolecularWeightTable::fsgMolecularWeight
fsgMolecularWeightTable::fsgMolecularWeights[fsgMolecularWeightTable::nMolecules] = 
{
    {"E",     0},
    {"C",    12.01115},
    {"H",     1.00797},
    {"O",    15.99940},
    {"N",    14.00670},
    {"S",    32.06400},
    {"F",    18.99840},
    {"Cl",   35.45300},
    {"Br",   79.90090},
    {"I",   126.90440},
    {"He",    4.00260},
    {"Ne",   20.18300},
    {"Ar",   39.94800},
    {"Kr",   83.80000},
    {"Xe",  131.30000},
    {"H2",    2.01594},
    {"D2",    4.02820},
    {"N2",   28.01340},
    {"O2",   31.99880},
    {"Air",  28.96440},
    {"CO",   28.01055},
    {"CO2",  44.00995},
    {"N2O",  44.01280},
    {"NH3",  17.03061},
    {"H2O",  18.01534},
    {"SF4", 108.05760},
    {"Cl2",  70.90600},
    {"Br2", 159.80180},
    {"SO2",  64.06280}
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fsgMolecularWeightTable::fsgMolecularWeightTable()
{
    for (int i=0; i<nMolecules; i++)
    {
        insert
        (
            word(fsgMolecularWeights[i].name),
            fsgMolecularWeights[i].molWeight
        );
    }
}


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

fsgMolecularWeightTable fsgMolecularWeights;

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //
