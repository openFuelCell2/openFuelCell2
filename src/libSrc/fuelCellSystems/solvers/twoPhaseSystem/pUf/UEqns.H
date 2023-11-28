Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

{
    momentumTransfer(); // !!! Update coefficients shouldn't be necessary
                              //     This should be done on demand

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(momentumTransferf());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    {
        U1Eqn =
        (
            phase1_.UfEqn()
         ==
           *momentumTransfer[phase1_.name()]
          + fvModels.source(alpha1, rho1, U1)
        );
        porousZone_->addResistance(U1Eqn);
        U1Eqn.relax();
        fvConstraints.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvConstraints.constrain(U1);
    }

    {
        porousZone_->dpcds(pc_, dpcds_);

        U2Eqn =
        (
            phase2_.UfEqn()
         ==
           *momentumTransfer[phase2_.name()]
          + fvModels.source(alpha2, rho2, U2)
          + alpha2 * (dpcds_ & fvc::grad(alpha2))
        );
        porousZone_->addResistance(U2Eqn);
        U2Eqn.relax();
        fvConstraints.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvConstraints.constrain(U2);
    }
}
