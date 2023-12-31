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

#include "phaseTransferModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseTransferModel> Foam::phaseTransferModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    word phaseTransferModelType(dict.lookup("type"));

    Info<< "Selecting phaseTransferModel for "
        << pair << ": " << phaseTransferModelType << endl;

    auto* ctorPtr =
        dictionaryConstructorTable(phaseTransferModelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown phaseTransferModelType type "
            << phaseTransferModelType << endl << endl
            << "Valid phaseTransferModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr(dict, pair);
}


// ************************************************************************* //
