if (!pimple_.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple_.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + MRF_.zeroFilter
        (
            fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi, rhoUf)
        )
    )
  + phig
);

fvc::makeRelative(phiHbyA, rho, U);
MRF_.makeRelative(fvc::interpolate(rho), phiHbyA);

volScalarField rhorAtU("rhorAtU", rho*rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh_, rho, U, phiHbyA, rhorAtU, MRF_);

{
    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p_rgh_)*mesh.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(p_rgh_);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh_))
      + fvc::div(phiHbyA)
     ==
        fvOptions(psi, p_rgh_, rho.name())
      + dmdt
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p_rgh_));

        pEqn.solve(mesh.solver(p_rgh_.select(pimple_.finalInnerIter())));

        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh_ + rho*gh, pMin);

    p_rgh_ = p - rho*gh;
}

//#include "rhoEqn.H"
//#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_rgh_.relax();

U = HbyA - rAtU*fvc::grad(p_rgh_);
U.correctBoundaryConditions();
fvOptions.correct(U);

thermo.correctRho(psi*p - psip0, rhoMin, rhoMax) ;
rho = thermo.rho();

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
