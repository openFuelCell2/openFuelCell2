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

// Drag coefficient
const surfaceScalarField Kdf("Kdf", this->Kdf());

// Diagonal coefficients
PtrList<surfaceScalarField> AFfs(this->AFfs());

PtrList<surfaceScalarField> rAUfs;
rAUfs.append
(
    new surfaceScalarField
    (
        IOobject::groupName("rAUf", phase1_.name()),
        1.0
       /(
            byDt(alphaRhof10)
          + fvc::interpolate(U1Eqn.A())
          + AFfs[0]
        )
    )
);
rAUfs.append
(
    new surfaceScalarField
    (
        IOobject::groupName("rAUf", phase2_.name()),
        1.0
       /(
            byDt(alphaRhof20)
          + fvc::interpolate(U2Eqn.A())
          + AFfs[1]
        )
    )
);
const surfaceScalarField& rAUf1 = rAUfs[0];
const surfaceScalarField& rAUf2 = rAUfs[1];

AFfs.clear();

// Explicit force fluxes
PtrList<surfaceScalarField> phiFfs(this->phiFfs(rAUfs));
const surfaceScalarField& phiFf1 = phiFfs[0];
const surfaceScalarField& phiFf2 = phiFfs[1];

while (pimple_.correct())
{
    volScalarField rho("rho", this->rho());     // rho = this->rho();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p - rho*gh;

    surfaceScalarField rhof1(fvc::interpolate(rho1));
    surfaceScalarField rhof2(fvc::interpolate(rho2));

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    correctBoundaryFlux();

    const surfaceScalarField alpharAUf1
    (
        IOobject::groupName("alpharAUf", phase1_.name()),
        max(alphaf1, phase1_.residualAlpha())*rAUf1
    );

    const surfaceScalarField alpharAUf2
    (
        IOobject::groupName("alpharAUf", phase2_.name()),
        max(alphaf2, phase2_.residualAlpha())*rAUf2
    );

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
          - alphaf2*(rhof1 - rhof2)*(g & mesh.Sf())
        )
      + phiFf1
    );

    const surfaceScalarField phigF2
    (
        alpharAUf2
       *(
            ghSnGradRho
          - alphaf1*(rhof2 - rhof1)*(g & mesh.Sf())
        )
      + phiFf2
    );

    // Predicted fluxes
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1_.name()),
        phi1
    );

    phiHbyA1 =
        rAUf1
       *(
            alphaRhof10*byDt(MRF_.absolute(phi1.oldTime()))
          + fvc::flux(U1Eqn.H())
        )
      - phigF1;

    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2_.name()),
        phi2
    );

    phiHbyA2 =
        rAUf2
       *(
            alphaRhof20*byDt(MRF_.absolute(phi2.oldTime()))
          + fvc::flux(U2Eqn.H())
        )
      - phigF2;

    // Drag fluxes
    PtrList<surfaceScalarField> phiKdPhifs(this->phiKdPhifs(rAUfs));
    const surfaceScalarField& phiKdPhif1 = phiKdPhifs[0];
    const surfaceScalarField& phiKdPhif2 = phiKdPhifs[1];

    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 - phiKdPhif1) + alphaf2*(phiHbyA2 - phiKdPhif2)
    );

    MRF_.makeRelative(phiHbyA);

    phiKdPhifs.clear();

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

            if (pEqnComp1.valid())
            {
                pEqn += pEqnComp1();
            }

            if (pEqnComp2.valid())
            {
                pEqn += pEqnComp2();
            }

            //pEqn.relax();

            pEqn.solve(mesh.solution().solverDict(p_rgh_.select(pimple_.finalInnerIter())));
        }

        if (pimple_.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            phi_ = phiHbyA + pEqnIncomp.flux();

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
                ((phi2s + rAUf2*Kdf*phi1s) - (phi1s + rAUf1*Kdf*phi2s))
               /(1 - rAUf1*rAUf2*sqr(Kdf))
            );

            phi1 = phi_ - alphaf2*phir;
            phi2 = phi_ + alphaf1*phir;

            U1 = fvc::reconstruct(MRF_.absolute(phi1));
            U1.correctBoundaryConditions();
            fvConstraints.constrain(U1);

            U2 = fvc::reconstruct(MRF_.absolute(phi2));
            U2.correctBoundaryConditions();
            fvConstraints.constrain(U2);

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
