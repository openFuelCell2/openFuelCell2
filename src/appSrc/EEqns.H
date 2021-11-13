{
    //- Solving temperature
    Info << "Solving for temperature:\n" << endl;

    fvScalarMatrix E1Eqn
    (
        fvm::ddt(rhoCpCell, TCell)
      + fvm::div(rhoCpPhiCell, TCell)
      + fvm::SuSp(-fvc::div(rhoCpPhiCell), TCell)
      - fvm::laplacian(kCell, TCell)
      - QdotCell
    );

    // Apply under-relaxation: not allowed
    // Under-relaxation in solid conduction (which dominates) creates a large
    // thermal inertia.  HJ, 25/May/2011
    E1Eqn.relax();

    E1Eqn.solve();

//     TCell.max(273);
//     TCell.min(373);

Info<< "T"
        << ": min = " << Foam::min(TCell.primitiveField())
        << ", mean = " << Foam::average(TCell.primitiveField())
        << ", max = " << Foam::max(TCell.primitiveField())
        << endl;
}
