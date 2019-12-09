{
    if (!phase_.pure() && Y)
    {
        UPtrList<volScalarField>& Y1 = phase_.YActiveRef();

        forAll(Y1, i)
        {
            fvScalarMatrix Y1iEqn
            (
                phase_.YiEqn(Y1[i])
             ==
                fvOptions(rho, Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve(mesh.solver("Yi"));

            Y1[i].max(SMALL);
            Y1[i].min(1.0);

            Info << "Y: " << Y1[i].name() << ":"
                 << gMin(Y1[i]) << " -> " << gMax(Y1[i])
                 << endl;
        }
    }
}
