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

Class
    Foam::fsgDiffusionVolumeTable

Description
    A table of Fuller-Schettler-Giddings diffusion volumes [cm^3]
    for available molecules

    data source:
        Fuller, Ensley, and Giddings, Diffusion of Halogenated Hydrocarbons
        in Helium, Journal of Physical Chemistry v73 n11 November 1969

SourceFiles
    fsgDiffusionVolumes.C

\*---------------------------------------------------------------------------*/

#ifndef fsgDiffusionVolumes_H
#define fsgDiffusionVolumes_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

/*---------------------------------------------------------------------------*\
                           Class fsgDiffusionVolumes Declaration
\*---------------------------------------------------------------------------*/

class fsgDiffusionVolumeTable
:
    public HashTable<scalar>
{

public:

    // Public types and data

        //- Structure to hold the molecule name and diffusion volume pair
        struct fsgDiffusionVolume
        {
            char name[4];
            scalar volume;
        };

        static const int nMolecules = 29;

        //- Static table of the volumes of all known molecules
        static const fsgDiffusionVolume fsgDiffusionVolumes[nMolecules];


    // Constructors

        //- Construct from fsgDiffusionVolumes_
        fsgDiffusionVolumeTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// Fuller-Schettler-Giddings diffusion volumes table for available molecules

extern fsgDiffusionVolumeTable fsgDiffusionVolumes;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
