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

    // We have to relax the equation sometimes...
    // But this greatly decreases convergence speed.
    E1Eqn.relax();

    E1Eqn.solve();

Info<< "T"
        << ": min = " << Foam::min(TCell.primitiveField())
        << ", mean = " << Foam::average(TCell.primitiveField())
        << ", max = " << Foam::max(TCell.primitiveField())
        << endl;
}
