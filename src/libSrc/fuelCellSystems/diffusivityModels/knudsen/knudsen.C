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

#include "knudsen.H"
#include "addToRunTimeSelectionTable.H"

#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diffusivityModels
    {
        defineTypeNameAndDebug(knudsen, 0);
        addToRunTimeSelectionTable(diffusivityModel, knudsen, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::knudsen::knudsen
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    const dictionary& dict
)
:
    diffusivityModel(name, mesh, diff, dict),
    Tname_(dict_.lookup("Tname")),
    dPore_(dict_.lookup("dPore")),
    MW_(dict_.lookup("MW"))
{}


Foam::diffusivityModels::knudsen::knudsen
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    word Tname,
    const dimensionedScalar& dPore,
    const dimensionedScalar& MW
)
:
    diffusivityModel(name, mesh, diff),
    Tname_(Tname),
    dPore_(dPore),
    MW_(MW)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::diffusivityModels::knudsen::writeData()
{
    if (firstIndex_)
    {
        Info<< "diffusivityModels::knudsen:" << nl
            << "dPore = " << dPore_ << nl
            << "MW = " << MW_ << endl;

        firstIndex_ = false;
    }
}


void Foam::diffusivityModels::knudsen::evaluate()
{
    //  D_{knudsen} = (poreDiameter/2)*97*sqrt(T/MW)
    //  where
    //      poreDiameter = [m]
    //      T ............ [K]
    //      MW ........... [kg/kmol]
    //  Geankoplis, Christie J, Transport Processes and Unit Operations,
    //  second edition (1983), Allyn and Bacon Series in Engineering,
    //  ISBN 0-205-07788-9, page 452.

    // When using OpenFOAM-1.6.x, ...
    const volScalarField& T =
        mesh_.lookupObject<volScalarField>(Tname_);

    forAll(cells_, i)
    {
        diff_[cells_[i]] =
            48.5*dPore_.value()
            *Foam::sqrt(T.internalField()[cells_[i]]/MW_.value());
    }
}

// ************************************************************************* //

