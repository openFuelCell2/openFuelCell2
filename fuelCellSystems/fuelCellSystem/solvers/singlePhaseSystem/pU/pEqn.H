volScalarField rAU
(
    "rAU",
    1.0
   /(
        UEqn.A()
      + byDt(rho)
    )
);

const surfaceScalarField rAUf
(
    fvc::interpolate(rAU)
);

// --- Pressure corrector loop
while (pimple_.correct())
{
    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF_.correctBoundaryFlux(U, phi);

    // Combined buoyancy and force fluxes
    const surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    const surfaceScalarField phigF
    (
        rAUf*ghSnGradRho
    );

    // Predicted velocities
    volVectorField HbyA
    (
        "HbyA",
        U
    );
    HbyA =
        rAU
       *(
            UEqn.H()
          + byDt(rho)*U.oldTime()
        );

    // Predicted fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA) - phigF
    );

    MRF_.makeRelative(phiHbyA);

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField() - phi.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp;
    if (phase_.compressible())
    {
        pEqnComp =
            (
                fvc::ddt(rho) + fvc::div(rhoPhi)
              - fvc::Sp(fvc::div(phi), rho)
            )/rho
          + (psi/rho)*correction(fvm::ddt(p_rgh_));
    }

    // Add option sources
    {
        if (fvOptions.appliesToField(rho.name()))
        {
            tmp<fvScalarMatrix> optEqn = fvOptions(rho);
            if (pEqnComp.valid())
            {
                pEqnComp.ref() -= (optEqn & rho)/rho;
            }
            else
            {
                pEqnComp = fvm::Su(- (optEqn & rho)/rho, p_rgh_);
            }
        }
    }

    // Add additional mass 
    // Mass transfer between regions
    {
        if (pEqnComp.valid())
        {
            pEqnComp.ref() -= dmdt/rho;
        }
        else
        {
            pEqnComp = fvm::Su(- dmdt/rho, p_rgh_);
        }
    }

    // Cache p prior to solve for density update
    const volScalarField p_rgh_0(p_rgh_);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh_)
        );

        {
            fvScalarMatrix pEqn(pEqnIncomp);

            if (pEqnComp.valid())
            {
                pEqn += pEqnComp();
            }

            ::solve
            (
                pEqn,
                mesh.solver(p_rgh_.select(pimple_.finalInnerIter()))
            );
        }

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            rhoPhi = phi*fvc::interpolate(rho);

            // Set the phase dilatation rates
            if (pEqnComp.valid())
            {
                phase_.divU(-pEqnComp & p_rgh_);
            }

            // Optionally relax pressure for velocity correction
            p_rgh_.relax();

            surfaceScalarField mSfGradp
            (
                "mSfGradp",
                pEqnIncomp.flux()/rAUf
            );

            U = HbyA + fvc::reconstruct(rAUf*mSfGradp - phigF);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh_ + rho*gh, pMin);

    // Limit p_rgh
    p_rgh_ = p - rho*gh;

    // Update densities from change in p_rgh
    rho += psi*(p_rgh_ - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho*gh;
    p_rgh_.correctBoundaryConditions();
}
