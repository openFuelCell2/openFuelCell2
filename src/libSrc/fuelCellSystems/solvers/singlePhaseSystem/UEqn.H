Info<< "Constructing momentum equations" << endl;

tmp<fvVectorMatrix> tUEqn
(
    phase_.UEqn()
 ==
    fvOptions(rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

porousZone_->addResistance(UEqn);
UEqn.relax();
fvOptions.constrain(UEqn);
U.correctBoundaryConditions();
fvOptions.correct(U);

if (pimple_.momentumPredictor())
{
    ::solve(UEqn == -fvc::grad(p_rgh_));

    fvOptions.correct(U);
}
