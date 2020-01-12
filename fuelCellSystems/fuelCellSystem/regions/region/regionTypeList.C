/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "regionTypeList.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "fuelCellSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypeList::regionTypeList
(
    const fvMesh& mesh
)
:
    PtrList<regionType>(),
    mesh_(mesh),
    region_(mesh.time())
{
    reset(region_);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypeList::~regionTypeList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionTypeList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No porosity models active" << endl;
    }

    return a;
}


void Foam::regionTypeList::reset(const regionProperties& rp)
{
    wordList regionNames;

    forAllConstIter(HashTable<wordList>, rp, iter)
    {
        const wordList& regions = iter();

        forAll(regions, regionI)
        {
            if (findIndex(regionNames, regions[regionI]))
            {
                regionNames.append(regions[regionI]);
            }
        }
    }

    this->setSize(regionNames.size());

    label i = 0;

    forAllConstIter(HashTable<wordList>, rp, iter)
    {
        const word& modelType = iter.key();
        const wordList& regions = iter();

        if (regions.size())
        {
            forAll(regions, regionI)
            {
                Info << "Creating " << regions[regionI] << endl;

                this->set
                (
                    i++,
                    regionType::New
                    (
                        mesh_,
                        regions[regionI],
                        modelType
                    )
                );
            }
        }
    }
}


void Foam::regionTypeList::correct()
{
    forAll(*this, i)
    {
        this->operator[](i).correct();
    }
}


void Foam::regionTypeList::setRDeltaT()
{
    forAll(*this, i)
    {
        this->operator[](i).setRDeltaT();
    }
}


void Foam::regionTypeList::solve()
{
    forAll(*this, i)
    {
        this->operator[](i).solve();
    }
}


void Foam::regionTypeList::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    forAll(*this, i)
    {
        this->operator[](i).mapToCell(fuelCell);
    }
}


void Foam::regionTypeList::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    forAll(*this, i)
    {
        this->operator[](i).mapFromCell(fuelCell);
    }
}

// ************************************************************************* //
