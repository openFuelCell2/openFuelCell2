const surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
const surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

PtrList<volScalarField> rAUs;
rAUs.append
(
    new volScalarField
    (
        IOobject::groupName("rAU", phase1_.name()),
        1.0
       /(
            U1Eqn.A()
          + byDt(max(phase1_.residualAlpha() - alpha1, scalar(0))*rho1)
        )
    )
);
rAUs.append
(
    new volScalarField
    (
        IOobject::groupName("rAU", phase2_.name()),
        1.0
       /(
            U2Eqn.A()
          + byDt(max(phase2_.residualAlpha() - alpha2, scalar(0))*rho2)
        )
    )
);
const volScalarField& rAU1 = rAUs[0];
const volScalarField& rAU2 = rAUs[1];

const surfaceScalarField alpharAUf1
(
    fvc::interpolate(max(alpha1, phase1_.residualAlpha())*rAU1)
);
const surfaceScalarField alpharAUf2
(
    fvc::interpolate(max(alpha2, phase2_.residualAlpha())*rAU2)
);

// Drag coefficients
const volScalarField Kd(this->Kd());
const volScalarField rAUKd1(rAU1*Kd);
const volScalarField rAUKd2(rAU2*Kd);
const surfaceScalarField phiKd1(fvc::interpolate(rAUKd1));
const surfaceScalarField phiKd2(fvc::interpolate(rAUKd2));

// Explicit force fluxes
PtrList<surfaceScalarField> phiFs(this->phiFs(rAUs));
const surfaceScalarField& phiF1 = phiFs[0];
const surfaceScalarField& phiF2 = phiFs[1];

// --- Pressure corrector loop
while (pimple_.correct())
{
    volScalarField rho("rho", this->rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    correctBoundaryFlux();

    // Combined buoyancy and force fluxes
    const surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    const surfaceScalarField phigF1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
        )
      + phiF1
    );

    const surfaceScalarField phigF2
    (
        alpharAUf2
       *(
           ghSnGradRho
         - alphaf1*fvc::interpolate(rho2 - rho1)*(g & mesh.Sf())
        )
      + phiF2
    );

    // Predicted velocities
    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1_.name()),
        U1
    );
    HbyA1 =
        rAU1
       *(
            U1Eqn.H()
          + byDt(max(phase1_.residualAlpha() - alpha1, scalar(0))*rho1)
           *U1.oldTime()
        );

    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2_.name()),
        U2
    );
    HbyA2 =
        rAU2
       *(
            U2Eqn.H()
         +  byDt(max(phase2_.residualAlpha() - alpha2, scalar(0))*rho2)
           *U2.oldTime()
        );

    // Correction force fluxes
    PtrList<surfaceScalarField> ddtCorrByAs(this->ddtCorrByAs(rAUs));

    // Predicted fluxes
    const surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1_.name()),
        fvc::flux(HbyA1) - phigF1 - ddtCorrByAs[0]
    );

    const surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2_.name()),
        fvc::flux(HbyA2) - phigF2 - ddtCorrByAs[1]
    );

    ddtCorrByAs.clear();

    // Drag fluxes
    PtrList<surfaceScalarField> phiKdPhis(this->phiKdPhis(rAUs));
    const surfaceScalarField& phiKdPhi1 = phiKdPhis[0];
    const surfaceScalarField& phiKdPhi2 = phiKdPhis[1];

    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 - phiKdPhi1) + alphaf2*(phiHbyA2 - phiKdPhi2)
    );

    MRF_.makeRelative(phiHbyA);

    phiKdPhis.clear();

    // Construct pressure "diffusivity"
    const surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
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

    // Add additional mass 
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

            if (pEqnComp1.valid())
            {
                pEqn += pEqnComp1();
            }

            if (pEqnComp2.valid())
            {
                pEqn += pEqnComp2();
            }

            pEqn.solve(mesh.solution().solverDict(p_rgh_.select(pimple_.finalInnerIter())));
        }

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp
            (
                "mSfGradp",
                pEqnIncomp.flux()/rAUf
            );

            // Partial-elimination phase-flux corrector
            {
                const surfaceScalarField phi1s
                (
                    phiHbyA1 + alpharAUf1*mSfGradp
                );

                const surfaceScalarField phi2s
                (
                    phiHbyA2 + alpharAUf2*mSfGradp
                );

                surfaceScalarField phir
                (
                    ((phi1s + phiKd1*phi2s) - (phi2s + phiKd2*phi1s))
                   /(1 - phiKd1*phiKd2)
                );

                phi1 = phi_ + alphaf2*phir;
                phi2 = phi_ - alphaf1*phir;
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

            // Optionally relax pressure for velocity correction
            p_rgh_.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            // Partial-elimination phase-velocity corrector
            {
                const volVectorField Us1
                (
                    HbyA1
                  + fvc::reconstruct(alpharAUf1*mSfGradp - phigF1)
                );

                const volVectorField Us2
                (
                    HbyA2
                  + fvc::reconstruct(alpharAUf2*mSfGradp - phigF2)
                );

                const volVectorField U
                (
                    alpha1*(Us1 + rAUKd1*U2) + alpha2*(Us2 + rAUKd2*U1)
                );

                const volVectorField Ur
                (
                    ((1 - rAUKd2)*Us1 - (1 - rAUKd1)*Us2)/(1 - rAUKd1*rAUKd2)
                );

                U1 = U + alpha2*Ur;
                U1.correctBoundaryConditions();
                fvConstraints.constrain(U1);

                U2 = U - alpha1*Ur;
                U2.correctBoundaryConditions();
                fvConstraints.constrain(U2);
            }
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh_ + rho*gh, pMin);

    // Limit p_rgh
    p_rgh_ = p - rho*gh;

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh_ - p_rgh_0);
    rho2 += psi2*(p_rgh_ - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho = this->rho();
    p_rgh_ = p - rho*gh;
    p_rgh_.correctBoundaryConditions();
}
