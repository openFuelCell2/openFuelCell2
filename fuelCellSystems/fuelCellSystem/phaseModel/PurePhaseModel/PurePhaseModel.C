/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "PurePhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PurePhaseModel<BasePhaseModel>::PurePhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index)
{
    //- Size 1
    iDmdt_.resize(1);

    iDmdt_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(word("iDmdt" + phaseName), this->water),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("iDmdt", dimDensity/dimTime, 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PurePhaseModel<BasePhaseModel>::~PurePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::PurePhaseModel<BasePhaseModel>::pure() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::PurePhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    FatalErrorInFunction
        << "Cannot construct a species fraction equation for a pure phase"
        << exit(FatalError);

    return tmp<fvScalarMatrix>();
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::Y() const
{
    // Y_ has never been set, so we are returning an empty list

    return Y_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::X() const
{
    // X_ has never been set, so we are returning an empty list

    return X_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::iDmdt() const
{
    return iDmdt_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::PurePhaseModel<BasePhaseModel>::Y(const word& name) const
{
    FatalErrorInFunction
        << "Cannot get a species fraction by name from a pure phase"
        << exit(FatalError);

    return NullObjectRef<volScalarField>();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::PurePhaseModel<BasePhaseModel>::X(const word& name) const
{
    FatalErrorInFunction
        << "Cannot get a species fraction by name from a pure phase"
        << exit(FatalError);

    return NullObjectRef<volScalarField>();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::PurePhaseModel<BasePhaseModel>::iDmdt(const word& name) const
{
    return iDmdt_[0];
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::YRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return Y_;
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::XRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return X_;
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::iDmdtRef()
{
    return iDmdt_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::YActive() const
{
    // Y_ has never been set, so we are returning an empty list

    return Y_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::XActive() const
{
    // X_ has never been set, so we are returning an empty list

    return X_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::YActiveRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return Y_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel>::XActiveRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return X_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PurePhaseModel<BasePhaseModel>::dmdt() const
{
    return iDmdt_[0];
}

// ************************************************************************* //
