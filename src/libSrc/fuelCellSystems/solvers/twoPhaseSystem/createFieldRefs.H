const volScalarField& alpha1 = phase1_;
const volScalarField& alpha2 = phase2_;

volVectorField& U1 = phase1_.URef();
surfaceScalarField& phi1 = phase1_.phiRef();
const tmp<surfaceScalarField> tAlphaPhi1(phase1_.alphaPhi());
const surfaceScalarField& alphaPhi1 = tAlphaPhi1();

volVectorField& U2 = phase2_.URef();
surfaceScalarField& phi2 = phase2_.phiRef();
const tmp<surfaceScalarField> tAlphaPhi2(phase2_.alphaPhi());
const surfaceScalarField& alphaPhi2 = tAlphaPhi2();

rhoThermo& thermo1 = phase1_.thermoRef();
rhoThermo& thermo2 = phase2_.thermoRef();

volScalarField& rho1 = thermo1.rho();
const volScalarField& psi1 = thermo1.psi();

volScalarField& rho2 = thermo2.rho();
const volScalarField& psi2 = thermo2.psi();

volScalarField dmdt1 = phase1_.dmdt();
volScalarField dmdt2 = phase2_.dmdt();
