surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

surfaceScalarField alphaRhof10
(
    "alphaRhof10",
    fvc::interpolate
    (
        max(alpha1.oldTime(), phase1_.residualAlpha())
       *rho1.oldTime()
    )
);

surfaceScalarField alphaRhof20
(
    "alphaRhof20",
    fvc::interpolate
    (
        max(alpha2.oldTime(), phase2_.residualAlpha())
       *rho2.oldTime()
    )
);

volScalarField rAU
(
    "rAU",
    1.0
   /(
        UEqn.A()
      + byDt(rho_)
    )
);

const surfaceScalarField rAUf
(
    fvc::interpolate(rAU)
);

// --- Pressure corrector loop
while (pimple.correct())
{
    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho_*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF_.correctBoundaryVelocity(U_);

    // Combined buoyancy and force fluxes
    const surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho_)*mesh.magSf()
    );

    const surfaceScalarField phigF
    (
        rAUf*ghSnGradRho
    );

    // Predicted velocities
    volVectorField HbyA
    (
        "HbyA",
        U_
    );
    HbyA =
        rAU
       *(
            UEqn.H()
          + byDt(rho_)*U_.oldTime()
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
            phiHbyA.boundaryField() - phi_.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp1, pEqnComp2;
    if (phase1_.compressible())
    {
        pEqnComp1 =
            (
                fvc::ddt(alpha1, rho1) + fvc::div(phase1_.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh_));
    }
    if (phase2_.compressible())
    {
        pEqnComp2 =
            (
                fvc::ddt(alpha2, rho2) + fvc::div(phase2_.alphaRhoPhi())
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh_));
    }

    tmp<fvScalarMatrix> pEqnComp;
    pEqnComp =
        (
            fvc::ddt(rho_) + fvc::div(rhoPhi_)
          - fvc::Sp(fvc::div(phi_), rho_)
        )/rho_
        + (scalar(1)/p_rgh_)*correction(fvm::ddt(p_rgh_));

    // Add option sources
    {
        if (fvModels.addsSupToField(rho1.name()))
        {
            tmp<fvScalarMatrix> optEqn1 = fvModels.source(alpha1, rho1);
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= (optEqn1 & rho1)/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- (optEqn1 & rho1)/rho1, p_rgh_);
            }
        }
        if (fvModels.addsSupToField(rho2.name()))
        {
            tmp<fvScalarMatrix> optEqn2 = fvModels.source(alpha2, rho2);
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= (optEqn2 & rho2)/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- (optEqn2 & rho2)/rho2, p_rgh_);
            }
        }

        if (fvModels.addsSupToField(rho_.name()))
        {
            tmp<fvScalarMatrix> optEqn = fvModels.source(rho_);
            pEqnComp.ref() -= (optEqn & rho_)/rho_;
        }
    }

    // Add additional mass 
    // Mass transfer between regions
    {
        volScalarField dmdt1 = phase1_.dmdt();
        volScalarField dmdt2 = phase2_.dmdt();

        if (pEqnComp1.valid())
        {
            pEqnComp1.ref() -= dmdt1/rho1;
        }
        else
        {
            pEqnComp1 = fvm::Su(- dmdt1/rho1, p_rgh_);
        }

        if (pEqnComp2.valid())
        {
            pEqnComp2.ref() -= dmdt2/rho2;
        }
        else
        {
            pEqnComp2 = fvm::Su(- dmdt2/rho2, p_rgh_);
        }

        pEqnComp.ref() -= (dmdt1 + dmdt2)/rho_;
    }

    // Cache p prior to solve for density update
    const volScalarField p_rgh_0(p_rgh_);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
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

            pEqn.solve(mesh.solution().solverDict(p_rgh_.select(pimple_.finalInnerIter())));
        }

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            rhoPhi_ = phi_*fvc::interpolate(rho_);

            // Optionally relax pressure for velocity correction
            p_rgh_.relax();

            surfaceScalarField mSfGradp
            (
                "mSfGradp",
                pEqnIncomp.flux()/rAUf
            );

            U_ = HbyA + fvc::reconstruct(rAUf*mSfGradp - phigF);
            U_.correctBoundaryConditions();
            fvConstraints.constrain(U_);

            //- Calculate phase velocities
            {
                surfaceScalarField phir(fvc::flux(UdmModel_->Udm()));

                if (phase1_.name() == this->continuous())
                {
                    phi2 = phi_ + phir;
                    phi1 = phi_ - alphaRhof20/alphaRhof10*phir;
                }

                if (phase2_.name() == this->continuous())
                {
                    phi1 = phi_ + phir;
                    phi2 = phi_ - alphaRhof10/alphaRhof20*phir;
                }

                U1 = fvc::reconstruct(MRF_.absolute(phi1));
                U2 = fvc::reconstruct(MRF_.absolute(phi2));

                //- Correct the velocities in the porous media
                porousZone_->correctU(U1, U2, U_, dpcds_);

                U1.correctBoundaryConditions();
                fvConstraints.constrain(U1);
                U2.correctBoundaryConditions();
                fvConstraints.constrain(U2);
           }

           // Set the phase dilatation rates
           if (pEqnComp1.valid())
           {
               phase1_.divU(-pEqnComp1 & p_rgh_);
           }
           if (pEqnComp2.valid())
           {
               phase2_.divU(-pEqnComp2 & p_rgh_);
           }

       }
    }

    // Update and limit the static pressure
    p = max(p_rgh_ + rho_*gh, pMin);

    // Limit p_rgh
    p_rgh_ = p - rho_*gh;

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh_ - p_rgh_0);
    rho2 += psi2*(p_rgh_ - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho_ = mixture_->rho();
    p_rgh_ = p - rho_*gh;
    p_rgh_.correctBoundaryConditions();
}
