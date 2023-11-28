{
    phasePairKey pair(phase1_.name(), phase2_.name());
    phaseModel& phase = phaseModels_[pair.other(continuous_)];

    const volScalarField& alpha = phase;

    volVectorField& U = phase.URef();

    rhoThermo& thermo = phase.thermoRef();

    volScalarField& rho = thermo.rho();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(heatTransfer());

    phaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();

    if (!phase.isothermal())
    {
        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvModels.source(alpha, rho, thermo.he())
        );

        EEqn.relax();
        fvConstraints.constrain(EEqn);
        EEqn.solve();
        fvConstraints.constrain(thermo.he());
    }
}
