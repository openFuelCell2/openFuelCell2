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

#include "masterRegionToCell.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterRegionToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, masterRegionToCell, word);
    addToRunTimeSelectionTable(topoSetSource, masterRegionToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::masterRegionToCell::usage_
(
    masterRegionToCell::typeName,
    "\n    Usage: masterRegionToCell zone\n\n"
    "    Select all cells in the cellZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterRegionToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            mesh_.time().timeName(),
            mesh_.time(),
            Foam::IOobject::MUST_READ
        )
    );

    forAll(mesh.cellZones(), i)
    {
        const cellZone& zone = mesh.cellZones()[i];

        if (zoneName_.match(zone.name()))
        {
            const labelList& cellLabels = mesh.cellZones()[i];

            Info<< "    Found matching zone " << zone.name()
                << " with " << cellLabels.size() << " cells." << endl;

            hasMatched = true;

            forAll(cellLabels, i)
            {
                // Only do active cells
                if (cellLabels[i] < mesh.nCells())
                {
                    if (cellMap_.found(cellLabels[i]))
                    {
                        addOrDelete(set, cellMap_[cellLabels[i]], add);
                    }
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningInFunction
            << "Cannot find any cellZone named " << zoneName_ << endl
            << "Valid names are " << mesh.cellZones().names() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterRegionToCell::masterRegionToCell
(
    const polyMesh& mesh,
    const word& zoneName
)
:
    topoSetSource(mesh),
    cellMapIO_
    (
        IOobject
        (
            "cellRegionAddressing",
            mesh.time().findInstance(mesh.meshDir(), "faces"),
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    ),
    zoneName_(zoneName)
{
    forAll(cellMapIO_, i)
    {
        cellMap_.insert
        (
            cellMapIO_[i],
            i
        );
    }
}


Foam::masterRegionToCell::masterRegionToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    masterRegionToCell(mesh, dict.lookup("name"))
{}


Foam::masterRegionToCell::masterRegionToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    masterRegionToCell(mesh, word(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterRegionToCell::~masterRegionToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::masterRegionToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells of cellZone " << zoneName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of cellZone " << zoneName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
