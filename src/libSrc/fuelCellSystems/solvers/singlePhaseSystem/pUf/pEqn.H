surfaceScalarField rhof
(
    "rhof",
    fvc::interpolate
    (
       rho.oldTime()
    )
);

surfaceScalarField rAUf
(
    "rAUf",
    1.0
  / (
        byDt(rhof)
      + fvc::interpolate(UEqn.A())
    )
);

while (pimple_.correct())
{

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho*gh;

    surfaceScalarField rhof0(fvc::interpolate(rho));

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
        rAUf * ghSnGradRho
    );

    // Predicted fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phi
    );

    phiHbyA =
        rAUf
       *(
            rhof*byDt(MRF_.absolute(phi.oldTime()))
          + fvc::flux(UEqn.H())
        )
      - phigF;

    MRF_.makeRelative(phiHbyA);

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - phi.boundaryField()
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
                pEqnComp = fvm::Su(-(optEqn & rho)/rho, p_rgh_);
            }
        }
    }

    // Add additional mass transfer
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
    volScalarField p_rgh_0("p_rgh_0", p_rgh_);

    while (pimple_.correctNonOrthogonal())
    {
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

        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            rhoPhi = phi*fvc::interpolate(rho);

            U = fvc::reconstruct(MRF_.absolute(phi));
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

