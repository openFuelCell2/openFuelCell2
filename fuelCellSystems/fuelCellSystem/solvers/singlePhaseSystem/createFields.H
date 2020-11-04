const uniformDimensionedVectorField& g =
  runTime.lookupObject<uniformDimensionedVectorField>("g");

#include "readhRef.H"
#include "gh.H"

volVectorField& U = phase_.URef();
surfaceScalarField& phi = phase_.phiRef();
// const surfaceScalarField& alphaPhi = phase_.alphaPhi();
surfaceScalarField& rhoPhi = phase_.alphaRhoPhiRef();

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

fv::options& fvOptions = this->fvOptions();
