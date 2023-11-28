if (!pimple_.simpleRho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh_));

if (pimple_.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF_.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

MRF_.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh_, rho, U, phiHbyA, rhorAUf, MRF_);

{
    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh_))
      + fvc::div(phiHbyA)
     ==
        fvModels.source(psi, p_rgh_, rho.name())
      + dmdt
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p_rgh_));

        pEqn.solve(mesh.solution().solverDict(p_rgh_.select(pimple_.finalInnerIter())));

        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh_ + rho*gh, pMin);

    p_rgh_ = p - rho*gh;
}

p_rgh_.relax();

U = HbyA - rAU*fvc::grad(p_rgh_);
U.correctBoundaryConditions();
fvConstraints.constrain(U);

//#include "rhoEqn.H"
//#include "compressibleContinuityErrs.H"

thermo.correctRho(psi*p - psip0) ;
//rho = thermo.rho();

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p_rgh_);
}
