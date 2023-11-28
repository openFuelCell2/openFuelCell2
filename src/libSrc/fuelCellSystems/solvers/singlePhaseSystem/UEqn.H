Info<< "Constructing momentum equations" << endl;

tmp<fvVectorMatrix> tUEqn
(
    phase_.UEqn()
 ==
    fvModels.source(rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

porousZone_->addResistance(UEqn);
UEqn.relax();
fvConstraints.constrain(UEqn);
U.correctBoundaryConditions();
fvConstraints.constrain(U);

if (pimple_.momentumPredictor())
{
    ::solve(UEqn == -fvc::grad(p_rgh_));

    fvConstraints.constrain(U);
}
