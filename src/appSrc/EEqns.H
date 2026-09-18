{
    //- Solving temperature
    Info << "Solving for temperature:\n" << endl;

    // Semi-implicit linearization of the electrochemical heat source.
    //
    // Qdot is strongly temperature dependent: the catalyst heat release is
    // j*(sign*eta - T*dS/(|n|F)), whose explicit temperature derivative is
    // exactly dQdotdT. That derivative is already assembled across all regions
    // by fuelCell.mapToCell() (see fluid::mapToCell and
    // activationOverpotentialModel::dQdotdT) and zeroed every step, but until
    // now it never entered the energy equation: Qdot was applied as a purely
    // explicit source. During an unconverged electrochemical coupling iteration
    // the reaction current transiently overshoots by orders of magnitude, so the
    // explicit source made the temperature matrix ill-conditioned and the T
    // solve diverged (residual growing instead of falling).
    //
    // Linearising Q(T) ~ Q_old + dQdT*(T - T_old) gives
    //     ... - laplacian(k,T) - dQdT*T = Q_old - dQdT*T_old
    // so the implicit part moves dQdT onto the diagonal (increasing diagonal
    // dominance, since dQdT = -j*dS/(|n|F) is negative for water splitting) and
    // the explicit correction -dQdT*T_old is added to the source. This does not
    // alter the converged solution; it only stabilises the march towards it.
    fvScalarMatrix E1Eqn
    (
        fvm::ddt(rhoCpCell, TCell)
      + fvm::div(rhoCpPhiCell, TCell)
      + fvm::SuSp(-fvc::div(rhoCpPhiCell), TCell)
      - fvm::laplacian(kCell, TCell)
      - QdotCell
      - fvm::Sp(dQdotdTCell, TCell)
      + dQdotdTCell*TCell.oldTime()
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
