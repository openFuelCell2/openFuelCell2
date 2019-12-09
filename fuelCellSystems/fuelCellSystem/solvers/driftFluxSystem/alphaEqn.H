{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const dictionary& alphaControls = mesh_.solverDict(alpha1.name());

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    bool LTS = fv::localEulerDdt::enabled(mesh_);

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    // Construct the dilatation rate source term
    tmp<volScalarField::Internal> tdgdt;

    if (phase1_.divU().valid() && phase2_.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1_.divU()()()
          - alpha1()
           *phase2_.divU()()()
        );
    }
    else if (phase1_.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1_.divU()()()
        );
    }
    else if (phase2_.divU().valid())
    {
        tdgdt =
        (
          - alpha1()
           *phase2_.divU()()()
        );
    }

    alpha1.correctBoundaryConditions();

    surfaceScalarField phir("phir", phi1 - phi2);

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dimless/dimTime, 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        if (tdgdt.valid())
        {
            scalarField& dgdt = tdgdt.ref();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] > 0.0)
                {
                    Sp[celli] -= dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                    Su[celli] += dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                }
                else if (dgdt[celli] < 0.0)
                {
                    Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
                }
            }
        }

        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phi_,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        phase1_.correctInflowOutflow(alphaPhi1);

        if (nAlphaSubCycles > 1)
        {
            tmp<volScalarField> trSubDeltaT;

            if (LTS)
            {
                trSubDeltaT =
                    fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
            }

            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhi10(alphaPhi1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhi10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    UniformField<scalar>(phase1_.alphaMax()),
                    zeroField()
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_.alphaPhiRef() = alphaPhi10;
                }
                else
                {
                    phase1_.alphaPhiRef() += alphaPhi10;
                }
            }

            phase1_.alphaPhiRef() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhi1,
                Sp,
                Su,
                UniformField<scalar>(phase1_.alphaMax()),
                zeroField()
            );

            phase1_.alphaPhiRef() = alphaPhi1;
        }

        phase1_.alphaRhoPhiRef() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhiRef() = phi_ - phase1_.alphaPhi();
        phase2_.correctInflowOutflow(phase2_.alphaPhiRef());
        phase2_.alphaRhoPhiRef() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.maxMin(SMALL, 1);

        // Update the phase-fraction of the other phase
        alpha2 = scalar(1) - alpha1;
        alpha2.maxMin(SMALL, 1);
    }
}
