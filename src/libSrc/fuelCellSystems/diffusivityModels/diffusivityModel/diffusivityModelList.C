/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
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

#include "diffusivityModelList.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Static data  * * * * * * * * * * * * * * //

const Foam::word Foam::diffusivityModelList::dictName("diffusivityModel");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModelList::diffusivityModelList
(
    const fvMesh& mesh,
    word phaseName
)
:
    PtrList<diffusivityModel>(),
    mesh_(mesh),
    dict_
    (
        IOdictionary 
        (
            IOobject
            (
                IOobject::groupName(dictName, phaseName),
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    ),
    diff_(mesh.nCells(), 0.0)
{
    reset(dict_);

    active();
}


Foam::diffusivityModelList::diffusivityModelList
(
    const fvMesh& mesh,
    const dictionary& dict,
    word phaseName
)
:
    PtrList<diffusivityModel>(),
    mesh_(mesh),
    dict_(dict),
    diff_(mesh.nCells(), 0.0)
{
    reset(dict);

    active();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusivityModelList::~diffusivityModelList()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::diffusivityModelList::active() const
{
    if (!this->size())
    {
        Info<< "No diffusivity models active" << endl;
    }

    return this->size();
}


void Foam::diffusivityModelList::reset(const dictionary& dict)
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
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            this->set
            (
                i++,
                diffusivityModel::New(name, mesh_, diff_, modelDict)
            );
        }
    }
}


void Foam::diffusivityModelList::correct
(
    const PtrList<volScalarField>& X,
    PtrList<volScalarField>& diffSp,
    label inertIndex
)
{
    forAll(X, i)
    {
        if (i == inertIndex)
        {
            continue;
        }

        Info<< "\nCorrect diffusivity for specie: " << X[i].group() << nl << endl;

        forAll(*this, j)
        {
            if(this->operator[](j).isFixed())
            {
                this->operator[](j).evaluate();

                this->operator[](j).writeData();

                // copy calculated zone diffusivity from diffModel to diffSp[i]
                for (auto k : this->operator[](j).cells())
                {
                    diffSp[i][k] = this->operator[](j).diff()[k];
                }
            }
            else if(!this->operator[](j).isBinary())
            {
                FatalErrorIn("ERROR: multicomponent diffusivity")
                    << "requires fixed or binary model"
                    << exit(FatalError);
            }
            else
            {
                // species diffusivities in mixture from binary diffusivities
                // ----------------------------------------------------------
                // D_{i} = (1-x[i])/sum_{l!=i}(x[l]/D[i,l])

                // pairwise binary diff calculation and accumulation
                // -------------------------------------------
                //initialize sum(x[l]/D[il]
                scalarField sum(diffSp[i].size(), 0);

                forAll(X, l)
                {
                    if (l != i)
                    {
                        this->operator[](j).setSpecies
                        (
                            X[i].group(),
                            X[l].group()
                        );
                        this->operator[](j).evaluate();
                        this->operator[](j).writeData();

                        for (auto k : this->operator[](j).cells())
                        {
                            if(this->operator[](j).diff()[k] != 0)
                            {
                                sum[k] += X[l][k]/this->operator[](j).diff()[k];
                            }
                        }
                    }
                }

                // diffSpIn[zone] <-- (1-x[a])/sum
                for (auto k : this->operator[](j).cells())
                {
                    if (sum[k] != 0)
                    {
                        diffSp[i][k] = (1 - X[i][k])/sum[k];
                    }
                }
            } //isBinary
        } //m
        diffSp[i].correctBoundaryConditions();
    } //!inert
}

// ************************************************************************* //
