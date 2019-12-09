Info<< "Constructing momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    {
        U1Eqn =
        (
            phase1_.UEqn()
         ==
           *momentumTransfer[phase1_.name()]
          + fvOptions(alpha1, rho1, U1)
        );
        porousZone_->addResistance(U1Eqn);
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        fvOptions.correct(U1);
    }

    {
        porousZone_->dpcds(pc_, dpcds_);

        U2Eqn =
        (
            phase2_.UEqn()
         ==
           *momentumTransfer[phase2_.name()]
          + fvOptions(alpha2, rho2, U2)
          + alpha2 * (dpcds_ & fvc::grad(alpha2))
        );
        porousZone_->addResistance(U2Eqn);
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        fvOptions.correct(U2);
    }
}
