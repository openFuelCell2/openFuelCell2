const uniformDimensionedVectorField& g = g_;

#include "readhRef.H"

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    *this
);

#include "gh.H"

volScalarField& p = phase1_.thermoRef().p();

Foam::fvModels& fvModels = this->fvModels(mesh);
Foam::fvConstraints& fvConstraints = this->fvConstraints(mesh);
