Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix UEqn(U, rho.dimensions()*U.dimensions()*dimVol/dimTime);

{
    {
        UEqn =
        (
            phase_.UfEqn()
         ==
            fvOptions(rho, U)
        );
        porousZone_->addResistance(UEqn);
        UEqn.relax();
        fvOptions.constrain(UEqn);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

