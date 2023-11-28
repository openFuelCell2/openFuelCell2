{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(this->massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    if (!phase1_.pure())
    {
        UPtrList<volScalarField>& Y1 = phase1_.YActiveRef();

        forAll(Y1, i)
        {
            fvScalarMatrix Y1iEqn
            (
                phase1_.YiEqn(Y1[i])
             ==
               *massTransfer[Y1[i].name()]
              + fvModels.source(alpha1, rho1, Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve("Yi");

            Y1[i].max(SMALL);
            Y1[i].min(1.0);

            Info << "Y: " << Y1[i].name() << ":"
                 << gMin(Y1[i]) << " -> " << gMax(Y1[i])
                 << endl;
        }
    }

    if (!phase2_.pure())
    {
        UPtrList<volScalarField>& Y2 = phase2_.YActiveRef();

        forAll(Y2, i)
        {
            fvScalarMatrix Y2iEqn
            (
                phase2_.YiEqn(Y2[i])
             ==
               *massTransfer[Y2[i].name()]
              + fvModels.source(alpha2, rho2, Y2[i])
            );

            Y2iEqn.relax();
            Y2iEqn.solve("Yi");

            Y2[i].max(SMALL);
            Y2[i].min(1.0);

            Info << "Y: " << Y2[i].name() << ":"
                 << gMin(Y2[i]) << " -> " << gMax(Y2[i])
                 << endl;
        }
    }
}
