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

surfaceScalarField rhof
(
    "rhof",
    fvc::interpolate
    (
       rho_.oldTime()
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

while (pimple.correct())
{

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho_*gh;

    surfaceScalarField rhof0(fvc::interpolate(rho_));

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
        rAUf * ghSnGradRho
    );

    // Predicted fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phi_
    );

    phiHbyA =
        rAUf
       *(
            rhof*byDt(MRF_.absolute(phi_.oldTime()))
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
          - phi_.boundaryField()
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

    // Add mass transfer
    {
        PtrList<volScalarField> dmdts(this->dmdts());
        if (dmdts.set(0))
        {
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= dmdts[0]/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- dmdts[0]/rho1, p_rgh_);
            }
        }
        if (dmdts.set(1))
        {
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= dmdts[1]/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- dmdts[1]/rho2, p_rgh_);
            }
        }
    }

    // Add additional mass transfer
    // Mass transfer between regions
    {
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
    volScalarField p_rgh_0("p_rgh_0", p_rgh_);

    while (pimple.correctNonOrthogonal())
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

            pEqn.solve(mesh.solution().solverDict(p_rgh_.select(pimple_.finalInnerIter())));
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            rhoPhi_ = phi_*fvc::interpolate(rho_);

            U_ = fvc::reconstruct(MRF_.absolute(phi_));
            U_.correctBoundaryConditions();
            fvConstraints.constrain(U_);

            surfaceScalarField phir(fvc::flux(UdmModel_->Udm()));

            //- Calculate phase velocities
            {
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

                //- Update velocities in porous media
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
