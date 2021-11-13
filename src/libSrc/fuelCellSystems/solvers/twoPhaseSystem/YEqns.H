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
              + fvOptions(alpha1, rho1, Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve(mesh.solver("Yi"));

            Y1[i].max(SMALL);
            Y1[i].min(1.0);
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
              + fvOptions(alpha2, rho2, Y2[i])
            );

            Y2iEqn.relax();
            Y2iEqn.solve(mesh.solver("Yi"));

            Y2[i].max(SMALL);
            Y2[i].min(1.0);
        }
    }
}
