/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
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

#include "nernstModel.H"
#include "phaseModel.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::nernstModel>
Foam::nernstModel::New
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
{
    word nernstModelType
    (
        word(dict.lookup("type"))
      + "<"
      + phase1.thermo().type()
      + ","
      + phase2.thermo().type()
      + ">"
    );

    Info<< "Selecting nernstModel for "
        << phase1.mesh().name() << ": "
        << nernstModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(nernstModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown nernstModelType type "
            << nernstModelType << endl << endl
            << "Valid nernstModelType types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(phase1, phase2, dict);
}


// ************************************************************************* //
