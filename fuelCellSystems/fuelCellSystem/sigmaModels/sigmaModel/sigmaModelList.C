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

#include "sigmaModelList.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaModelList::sigmaModelList
(
    const fvMesh& mesh
)
:
    PtrList<sigmaModel>(),
    mesh_(mesh),
    dict_
    (
        IOdictionary 
        (
            IOobject
            (
                "sigmaModels",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    )
{
    reset(dict_);

    active(true);
}


Foam::sigmaModelList::sigmaModelList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<sigmaModel>(),
    mesh_(mesh),
    dict_(dict)
{
    reset(dict_);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmaModelList::~sigmaModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sigmaModelList::active(const bool warn) const
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


void Foam::sigmaModelList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
//             const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            this->set
            (
                i++,
                sigmaModel::New(mesh_, modelDict)
            );
        }
    }
}


void Foam::sigmaModelList::correct
(
    volScalarField& sigmaField
)
{
    forAll(*this, i)
    {
        this->operator[](i).correct(sigmaField);
    }
}

// ************************************************************************* //
