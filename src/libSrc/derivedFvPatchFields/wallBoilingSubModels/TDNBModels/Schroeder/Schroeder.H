/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd
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

Class
    Foam::wallBoilingModels:TDNBModels:::Schroeder

Description
    Departure from nulceate boiling correlation.

    References:
    \verbatim
        Schroeder-Richter D. and Bartsch G. Analytical calculation of
        DNB-superheating by a pos-tulated thermo-mechanical effect of
        nucleate boiling.
        International Journal of Multiphase Flow, 20(6):1143–1167, 1994.

    \endverbatim


    \verbatim
        THEORETICAL CRITICAL HEAT FLUX PREDICTION BASEDNON-EQUILIBRIUM
        THERMODYNAMICS CONSIDERATIONSTHE SUBCOOLED BOILING PHENOMENON

        Germán Thelera and Daniel Freisba TECNA
        Estudios y Proyectos de Ingenierı́a S.A.
        Encarnación Ezcurra 365, C1107CLA Buenos Aires, Argentina
        Westinghouse
        Electric Germany GmbH
        Dudenstraße 44, 68167 Mannheim, Germany
    \endverbatim


SourceFiles
    Schroeder.C

\*---------------------------------------------------------------------------*/

#ifndef Schroeder_H
#define Schroeder_H

#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{

/*---------------------------------------------------------------------------*\
                         Class Schroeder Declaration
\*---------------------------------------------------------------------------*/

class Schroeder
:
    public TDNBModel
{

    // Private data:

        //- Isoentropic expansion factor for ideal gases
        // 5/3 monoatomic
        // 7/5 diatomic
        scalar kg_;

public:

    //- Runtime type information
    TypeName("Schroeder");

    // Constructors

        //- Construct from a dictionary
        Schroeder(const dictionary& dict);


    //- Destructor
    virtual ~Schroeder();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> TDNB
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;


        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TDNBModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
