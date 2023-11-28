#include "createRDeltaT.H"

const uniformDimensionedVectorField& g = g_;

#include "readhRef.H"
#include "gh.H"

volVectorField& U = phase_.URef();
//surfaceScalarField& phi = phase_.phiRef();
// const surfaceScalarField& alphaPhi = phase_.alphaPhi();
surfaceScalarField& phi = phase_.alphaRhoPhiRef();

rhoThermo& thermo = phase_.thermoRef();

volScalarField& p = thermo.p();
volScalarField rho = thermo.rho();
const volScalarField& psi = thermo.psi();
volScalarField dmdt = phase_.dmdt();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    *this
);

dimensionedScalar rhoMin = dimensionedScalar::lookupOrDefault
(
    "rhoMin",
    pimple_.dict(),
    dimDensity,
    Zero
);

dimensionedScalar rhoMax = dimensionedScalar::lookupOrDefault
(
    "rhoMax",
    pimple_.dict(),
    dimDensity,
    GREAT
);

p_rgh_ = p - rho*gh;

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Foam::fvModels& fvModels = this->fvModels(mesh);
Foam::fvConstraints& fvConstraints = this->fvConstraints(mesh);