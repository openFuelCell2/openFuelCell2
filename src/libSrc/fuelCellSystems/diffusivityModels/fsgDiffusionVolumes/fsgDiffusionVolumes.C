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

#include "fsgDiffusionVolumes.H"

namespace Foam
{
namespace diffusivityModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// [cm^3]

const fsgDiffusionVolumeTable::fsgDiffusionVolume
fsgDiffusionVolumeTable::fsgDiffusionVolumes[fsgDiffusionVolumeTable::nMolecules] = 
{
    {"E",     0},
    {"C",    15.9},
    {"H",     2.31},
    {"O",     6.11},
    {"N",     4.54},
    {"S",    22.9},
    {"F",    14.7},
    {"Cl",   21.0},
    {"Br",   21.9},
    {"I",    29.8},
    {"He",    2.67},
    {"Ne",    5.98},
    {"Ar",   16.2},
    {"Kr",   24.4},
    {"Xe",   32.7},
    {"H2",    6.12},
    {"D2",    6.84},
    {"N2",   18.5},
    {"O2",   16.3},
    {"Air",  19.7},
    {"CO",   18.0},
    {"CO2",  26.7},
    {"N2O",  35.9},
    {"NH3",  20.7},
    {"H2O",  13.1},
    {"SF4",  71.3},
    {"Cl2",  38.4},
    {"Br2",  69.0},
    {"SO2",  41.8}
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fsgDiffusionVolumeTable::fsgDiffusionVolumeTable()
{
    for (int i=0; i<nMolecules; i++)
    {
        insert
        (
            word(fsgDiffusionVolumes[i].name),
            fsgDiffusionVolumes[i].volume
        );
    }
}


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

fsgDiffusionVolumeTable fsgDiffusionVolumes;

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //
