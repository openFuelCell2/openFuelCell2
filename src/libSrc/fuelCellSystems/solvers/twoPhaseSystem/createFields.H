const uniformDimensionedVectorField& g =
  runTime.lookupObject<uniformDimensionedVectorField>("g");

#include "readhRef.H"

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    *this
);

#include "gh.H"

volScalarField& p = phase1_.thermoRef().p();

fv::options& fvOptions = this->fvOptions();
