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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porousZones::powerLaw::apply
(
    scalarField& Udiag,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];

            Udiag[cellI] +=
                V[cellI]*rho[cellI]*C0*pow(magSqr(U[cellI]), C1m1b2);
        }
    }
}


template<class RhoFieldType>
void Foam::porousZones::powerLaw::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];

            AU[cellI] =
                AU[cellI] + I*(rho[cellI]*C0*pow(magSqr(U[cellI]), C1m1b2));
        }
    }
}


// ************************************************************************* //
