/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "regionType.H"
#include "fuelCellSystem.H"

namespace Foam
{
    defineTypeNameAndDebug(regionType, 0);
    defineRunTimeSelectionTable(regionType, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::regionType::regionType
(
    const fvMesh& mesh,
    const word& regionName
)
:
    fvMesh
    (
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ
        )
    ),

    mesh_(mesh),

    dict_
    (
        IOobject
        (
            "regionProperties",
            this->time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    faceRegionAddressingIO_
    (
        IOobject
        (
            "faceRegionAddressing",
            this->time().findInstance(this->meshDir(), "faces"),
            polyMesh::meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    ),

    cellMapIO_
    (
        IOobject
        (
            "cellRegionAddressing",
            this->time().findInstance(this->meshDir(), "faces"),
            polyMesh::meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    ),

    patchesMapIO_
    (
        IOobject
        (
            "boundaryRegionAddressing",
            this->time().findInstance(this->meshDir(), "faces"),
            polyMesh::meshSubDir,
            *this,
            IOobject::MUST_READ
        )
    ),

    faceMap_(faceRegionAddressingIO_.size(), 1),

    faceMask_(faceRegionAddressingIO_.size(), 1)
{

    forAll(faceRegionAddressingIO_, i)
    {
        faceRegionAddressing_.insert
        (
            faceRegionAddressingIO_[i],
            i
        );
    }

    forAll(cellMapIO_, i)
    {
        cellMap_.insert
        (
            cellMapIO_[i],
            i
        );
    }

    forAll(patchesMapIO_, i)
    {
        patchesMap_.insert
        (
            patchesMapIO_[i],
            i
        );
    }

    forAll(faceMap_, i)
    {
        faceMap_[i] = mag(faceRegionAddressingIO_[i]) - 1;
        faceMask_[i] = sign(faceRegionAddressingIO_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionType::~regionType()
{}


// ************************************************************************* //
