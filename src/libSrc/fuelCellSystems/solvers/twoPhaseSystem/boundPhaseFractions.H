/*---------------------------------------------------------------------------*\
  Complementary physical phase fractions; residualAlpha regularizes closures.
\*---------------------------------------------------------------------------*/

#ifndef boundPhaseFractions_H
#define boundPhaseFractions_H

#include "fvCFD.H"

namespace Foam
{

// Emergency bounding is not a conservative flux correction. Report its volume
// change explicitly, including changes to the complementary phase. These are
// mesh-volume integrals, not pore-volume inventories or time-integrated sources.
inline void boundPhaseFractions
(
    volScalarField& alpha1,
    volScalarField& alpha2,
    const scalar alphaMax1,
    const scalar alphaMax2,
    const bool report = false
)
{
    const scalar lower = max(scalar(0), scalar(1) - alphaMax2);
    const scalar upper = min(scalar(1), alphaMax1);

    if (lower > upper || alphaMax1 < 0 || alphaMax2 < 0)
    {
        FatalErrorInFunction
            << "Incompatible alphaMax values: " << alphaMax1
            << " and " << alphaMax2 << exit(FatalError);
    }

    const scalarField before1(alpha1.primitiveField());
    const scalarField before2(alpha2.primitiveField());

    alpha1.clamp_range(lower, upper);
    alpha2 = scalar(1) - alpha1;

    const scalarField& V = alpha1.mesh().V();
    const scalarField delta1(alpha1.primitiveField() - before1);
    const scalarField delta2(alpha2.primitiveField() - before2);
    const scalar change1 = gSum(V*delta1);
    const scalar change2 = gSum(V*delta2);
    const scalar absoluteChange1 = gSum(V*mag(delta1));
    const scalar absoluteChange2 = gSum(V*mag(delta2));

    // alpha2 is normally stale here: its complement update is expected. Only
    // an appreciable correction to the transported alpha warrants a warning.
    if (report || absoluteChange1 > 1e-12*gSum(V))
    {
        Info<< "Phase fraction bounding: " << alpha1.name()
            << " deltaV=" << change1 << " L1deltaV=" << absoluteChange1
            << "; " << alpha2.name()
            << " updateDeltaV=" << change2
            << " updateL1deltaV=" << absoluteChange2
            << " [m3]; maxSumError="
            << gMax(mag(alpha1.primitiveField()
                      + alpha2.primitiveField() - scalar(1))) << endl;
    }
}

} // End namespace Foam

#endif
