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

#include "porousZoneList.H"
#include "volFields.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Static data  * * * * * * * * * * * * * * //

const Foam::word Foam::porousZoneList::dictName("porousZones");

// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

Foam::IOobject Foam::porousZoneList::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        dictName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Creating porous model list from " << io.name() << nl << endl;

        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        Info<< "No porous models present... " << nl << endl;

        io.readOpt(IOobject::NO_READ);
    }

    return io;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZoneList::porousZoneList
(
    const fvMesh& mesh
)
:
    PtrList<porousZone>(),
    mesh_(mesh),
    dict_(createIOobject(mesh))
{
    reset(dict_);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousZoneList::~porousZoneList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::porousZoneList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No porosity models active" << endl;
    }

    return a;
}


void Foam::porousZoneList::reset(const dictionary& dict)
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
                porousZone::New(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::porousZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        porousZone& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::porousZoneList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::porousZoneList::addResistance
(
    fvVectorMatrix& UEqn
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}


void Foam::porousZoneList::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, rho, mu);
    }
}


void Foam::porousZoneList::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, AU, correctAUprocBC);
    }
}

void Foam::porousZoneList::dpcds
(
    volTensorField& pc,
    volTensorField& dpcds
)
{
    forAll(*this, i)
    {
        this->operator[](i).dpcds(pc, dpcds);
    }
}


void Foam::porousZoneList::correctU
(
    volVectorField& U1,
    volVectorField& U2,
    const volVectorField& U,
    const volTensorField& dpcds
)
{
    forAll(*this, i)
    {
        this->operator[](i).correctU(U1, U2, U, dpcds);
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const porousZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
