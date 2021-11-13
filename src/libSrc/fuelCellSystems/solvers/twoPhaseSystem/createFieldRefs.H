const volScalarField& alpha1 = phase1_;
const volScalarField& alpha2 = phase2_;

volVectorField& U1 = phase1_.URef();
surfaceScalarField& phi1 = phase1_.phiRef();
const surfaceScalarField& alphaPhi1 = phase1_.alphaPhi();

volVectorField& U2 = phase2_.URef();
surfaceScalarField& phi2 = phase2_.phiRef();
const surfaceScalarField& alphaPhi2 = phase2_.alphaPhi();

rhoThermo& thermo1 = phase1_.thermoRef();
rhoThermo& thermo2 = phase2_.thermoRef();

volScalarField& rho1 = thermo1.rho();
const volScalarField& psi1 = thermo1.psi();

volScalarField& rho2 = thermo2.rho();
const volScalarField& psi2 = thermo2.psi();

volScalarField dmdt1 = phase1_.dmdt();
volScalarField dmdt2 = phase2_.dmdt();
