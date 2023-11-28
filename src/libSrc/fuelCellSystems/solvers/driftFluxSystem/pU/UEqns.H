Info<< "Constructing momentum equations" << endl;

fvVectorMatrix UEqn(U_, rho_.dimensions()*U_.dimensions()*dimVol/dimTime);

{
    {
        porousZone_->dpcds(pc_, dpcds_);

        UEqn =
        (
           fvm::ddt(rho_, U_) + fvm::div(rhoPhi_, U_)
         + MRF_.DDt(rho_, U_)
         + fvc::div(UdmModel_->tauDm())
         + turbulence_->divDevTau(U_)
         ==
           fvModels.source(rho_, U_)
         + alpha2 * (dpcds_ & fvc::grad(alpha2))
        );
        porousZone_->addResistance(UEqn);
        UEqn.relax();
        fvConstraints.constrain(UEqn);
        fvConstraints.constrain(U_);
    }
}
