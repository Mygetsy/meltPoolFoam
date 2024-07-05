using constant::mathematical::twoPi;
using constant::physicoChemical::R;     // universal gas constant
using constant::physicoChemical::sigma; // Stefan--Boltzmann constant

radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
evaporativeCooling =
    evaporationCoeff*mixture.vapourPressure(ambientPressure)*Hvap*Foam::sqrt(molarMass/twoPi/R/T);

const surfaceScalarField rhof = fvc::interpolate(rho, "interpolate(rho)");
const surfaceScalarField phiRhoLBeta = fvc::interpolate(liquidFraction, "interpolate(liquidFraction)")*mixture.thermo().rhoLiquid()*mixture.thermo().betaLiquid();
const surfaceScalarField kappaByRhoCp = mixture.kappaf()/(
                                     mixture.Cpf()*rhof + phiRhoLBeta*(fvc::interpolate(mixture.thermo().hLiq(T, alpha1), "interpolate(h)") - fvc::interpolate(h, "interpolate(h)")));


// Test code start
// const volScalarField phiRhoLBeta = liquidFraction*mixture.thermo().rhoLiquid()*mixture.thermo().betaLiquid();
// const volScalarField kappaByRhoCp = mixture.kappa()/(mixture.Cp()*rho + phiRhoLBeta*(mixture.thermo().hLiq(T, alpha1) - h));
// const volScalarField hDeltaRhoSolidLiquid = (mixture.thermo().rhoSolid() - mixture.thermo().rhoLiquid())*h
//                                             + (mixture.thermo().rhoLiquid() * mixture.thermo().hLiq(T, alpha1))
//                                             -(mixture.thermo().rhoSolid() * mixture.thermo().hSol(T, alpha1));
// const volScalarField hDeltaRhoGasSolid = (mixture.thermo().rhoGas() - mixture.thermo().rhoSolid())*h
//                                          - (mixture.thermo().rhoGas() * mixture.thermo().hGas(T, alpha1))
//                                          - (mixture.thermo().rhoSolid() * mixture.thermo().hSol(T, alpha1));
//
// const surfaceScalarField kappaByRhoCphDeltaRhoSolidLiquid = fvc::interpolate(kappaByRhoCp*hDeltaRhoSolidLiquid, "interpolate(kappa)");
// const surfaceScalarField kappaByRhoCphDeltaRhoGasSolid = fvc::interpolate(kappaByRhoCp*hDeltaRhoGasSolid, "interpolate(kappa)");
// const surfaceScalarField rhokappaByRhoCp = fvc::interpolate(rho*kappaByRhoCp, "interpolate(kappa)");
// Test code end

const volScalarField& redistribution = mixture.surfaceHeatSourceRedistribution();
const volScalarField magGradAlphaM = mag(mixture.gradAlphaM());
const volScalarField evaporativeCoolingPrime = evaporativeCooling*molarMass*Hvap/R/Cp/sqr(T);
const surfaceScalarField hDeltaRhoSolidLiquid = (mixture.thermo().rhoSolid() - mixture.thermo().rhoLiquid())*fvc::interpolate(h, "interpolate(h)")
                                                + ((mixture.thermo().rhoLiquid() * fvc::interpolate(mixture.thermo().hLiq(T, alpha1), "interpolate(h)"))
                                                -(mixture.thermo().rhoSolid() * fvc::interpolate(mixture.thermo().hSol(T, alpha1), "interpolate(h)")));

const surfaceScalarField hDeltaRhoGasSolid = (mixture.thermo().rhoGas() - mixture.thermo().rhoSolid())*fvc::interpolate(h, "interpolate(h)")
                                              - ((mixture.thermo().rhoGas() * fvc::interpolate(mixture.thermo().hGas(T, alpha1), "interpolate(h)"))
                                              -(mixture.thermo().rhoSolid() * fvc::interpolate(mixture.thermo().hSol(T, alpha1), "interpolate(h)")));



energySource =
(
   - fvc::laplacian(kappaByRhoCp*hDeltaRhoSolidLiquid, liquidFraction) // Initial
//  - fvc::laplacian(kappaByRhoCphDeltaRhoSolidLiquid, liquidFraction)
// + fvc::laplacian(kappaByCp*mixture.HsPrimeAlphaGf(), alpha1)
  - fvc::laplacian(kappaByRhoCp*hDeltaRhoGasSolid, alpha1) // Initial
//  - fvc::laplacian(kappaByRhoCphDeltaRhoGasSolid, alpha1)
  + (*laserHeatSource)()
  - redistribution*magGradAlphaM*(radiativeCooling + evaporativeCooling)
)/rho;

pp(T);

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + fvm::div(rhoPhi, h)
  - fvm::Sp(contErr, h)
 ==
    fvm::laplacian(rhof*kappaByRhoCp, h) // Initial
//    fvm::laplacian(rhokappaByRhoCp, h)
//  - fvm::laplacian(Hfus*kappaByCp*liquidFractionPrimeEnthalpyf, h)
//  + fvc::laplacian(Hfus*kappaByCp*liquidFractionPrimeEnthalpyf, h)
  + rho*energySource
  + redistribution*magGradAlphaM
   *(evaporativeCoolingPrime*h - fvm::Sp(evaporativeCoolingPrime, h))
  //+ fvc::div(rho*U & turbulence->devReff())
);
hEqn.solve();
pp(h);
pp(T);

if (gMin(h) < 0)
{
    runTime.writeNow();
    returnReduceOr(true);   // MPI barrier
    FatalError
        << "Enthalpy becomes negative!"
        << exit(FatalError);
}

mixture.correctThermo();
