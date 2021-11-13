Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix UEqn(U_, rho_.dimensions()*U_.dimensions()*dimVol/dimTime);

{
    {
        porousZone_->dpcds(pc_, dpcds_);

        UEqn =
        (
           fvm::ddt(rho_, U_) + fvm::div(rhoPhi_, U_)
         + MRF_.DDt(rho_, U_)
         + fvc::div(UdmModel_->tauDm())
         + turbulence_->divDevRhoReff(U_)
         ==
           fvOptions(rho_, U_)
         + alpha2 * (dpcds_ & fvc::grad(alpha2))
        );
        porousZone_->addResistance(UEqn);
        UEqn.relax();
        fvOptions.constrain(UEqn);
        U_.correctBoundaryConditions();
        fvOptions.correct(U_);
    }
}
