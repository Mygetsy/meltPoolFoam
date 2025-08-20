const dimensionedScalar rhoSol = mixture.thermo().rhoSolid();
const dimensionedScalar rhoGas = mixture.thermo().rhoGas();
const dimensionedScalar rhoLiq =  mixture.thermo().rhoLiquid();
const dimensionedScalar rhoLiqPrime = mixture.thermo().rhoLiquid()*mixture.thermo().betaLiquid();
const surfaceScalarField rhof = fvc::interpolate(rho, "interpolate(rho)");

const surfaceScalarField h_f = fvc::interpolate(h, "interpolate(h)");
const surfaceScalarField hLiq_f = fvc::interpolate(mixture.thermo().hLiq(T, alpha1), "interpolate(h)");
const surfaceScalarField hSol_f = fvc::interpolate(mixture.thermo().hSol(T, alpha1), "interpolate(h)");
const surfaceScalarField hGas_f = fvc::interpolate(mixture.thermo().hGas(T, alpha1), "interpolate(h)");
const surfaceScalarField liquidFraction_f = fvc::interpolate(liquidFraction, "interpolate(liquidFraction)");

const surfaceScalarField phiRhoLH = rhoLiqPrime*(hLiq_f - h_f);
const surfaceScalarField kappaByRhoCp = mixture.kappaf()/(mixture.Cpf()*rhof + liquidFraction_f*phiRhoLH);

const surfaceScalarField hDeltaRhoSolidLiquid = rhoSol*(h_f - hSol_f) - rhoLiq*(h_f - hLiq_f);
const surfaceScalarField hDeltaRhoGasSolid = rhoGas*(h_f - hGas_f) - rhoSol*(h_f - hSol_f);
