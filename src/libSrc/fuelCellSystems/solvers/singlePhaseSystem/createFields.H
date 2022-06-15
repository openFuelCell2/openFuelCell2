#include "createRDeltaT.H"

const uniformDimensionedVectorField& g =
  runTime.lookupObject<uniformDimensionedVectorField>("g");

#include "readhRef.H"
#include "gh.H"

volVectorField& U = phase_.URef();
//surfaceScalarField& phi = phase_.phiRef();
// const surfaceScalarField& alphaPhi = phase_.alphaPhi();
surfaceScalarField& phi = phase_.alphaRhoPhiRef();

rhoThermo& thermo = phase_.thermoRef();

volScalarField& p = thermo.p();
volScalarField& rho = thermo.rho();
const volScalarField& psi = thermo.psi();
volScalarField dmdt = phase_.dmdt();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    *this
);

dimensionedScalar rhoMin
(
    "rhoMin",
    dimDensity,
    Zero,
    pimple_.dict()
);

dimensionedScalar rhoMax
(
    "rhoMax",
    dimDensity,
    GREAT,
    pimple_.dict()
);

p_rgh_ = p - rho*gh;

fv::options& fvOptions = this->fvOptions();

#include "createDpdt.H"
