/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

inline const Foam::word& Foam::porousZone::name() const
{
    return name_;
}


inline bool Foam::porousZone::active() const
{
    return active_;
}


inline const Foam::labelList& Foam::porousZone::zoneIds() const
{
    return cellZoneIDs_;
}


inline const Foam::wordRe& Foam::porousZone::zoneName() const
{
    return zoneName_;
}


inline const Foam::dictionary& Foam::porousZone::dict() const
{
    return dict_;
}


inline Foam::scalar Foam::porousZone::porosity() const
{
    return porosity_;
}
// ************************************************************************* //
