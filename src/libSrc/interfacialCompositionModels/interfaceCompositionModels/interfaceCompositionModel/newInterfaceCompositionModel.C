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

#include "interfaceCompositionModel.H"
#include "phasePair.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceCompositionModel>
Foam::interfaceCompositionModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    word interfaceCompositionModelType
    (
        word(dict.lookup("type"))
      + "<"
      + pair.phase1().thermo().type()
      + ","
      + pair.phase2().thermo().type()
      + ">"
    );

    Info<< "Selecting interfaceCompositionModel for "
        << pair << ": " << interfaceCompositionModelType << endl;

    auto* ctorPtr =
        dictionaryConstructorTable(interfaceCompositionModelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown interfaceCompositionModelType type "
            << interfaceCompositionModelType << endl << endl
            << "Valid interfaceCompositionModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr(dict, pair);
}


// ************************************************************************* //
