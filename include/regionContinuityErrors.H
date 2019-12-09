{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
    (
        fvc::domainIntegrate(mag(rho - thermo.rho()))/totalMass
    ).value();

    scalar globalContErr =
    (
        fvc::domainIntegrate(rho - thermo.rho())/totalMass
    ).value();

    cumulativeContErr[i] += globalContErr;

    Info<< "time step continuity errors (" << mesh.name() << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr[i]
        << endl;
}

{
    volScalarField contErr(fvc::div(phi));

    scalar sumLocalContErr = runTime.deltaTValue()*
        mag(contErr)().weightedAverage(mesh.V()).value();

    scalar globalContErr = runTime.deltaTValue()*
        contErr.weightedAverage(mesh.V()).value();
    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
}
