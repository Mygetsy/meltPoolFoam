using constant::mathematical::twoPi;
using constant::physicoChemical::R;     // universal gas constant
using constant::physicoChemical::sigma; // Stefan--Boltzmann constant

radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
evaporativeCooling =
    evaporationCoeff*mixture.vapourPressure(ambientPressure)*Hvap*Foam::sqrt(molarMass/twoPi/R/T);

const surfaceScalarField kappaByCp = mixture.kappaf()/mixture.Cpf();
const volScalarField& redistribution = mixture.surfaceHeatSourceRedistribution();
const volScalarField magGradAlphaM = mag(mixture.gradAlphaM());
const volScalarField evaporativeCoolingPrime = evaporativeCooling*molarMass*Hvap/R/Cp/sqr(T);
const surfaceScalarField hDeltaRhoSolidLiquidByRho = ((mixture.thermo().rhoSolid() - mixture.thermo().rhoLiquid())*fvc::interpolate(h, "interpolate(h)")
                                      + ((mixture.thermo().rhoLiquid() * fvc::interpolate(mixture.thermo().hLiq(T, alpha1), "interpolate(h)"))
                                        -(mixture.thermo().rhoSolid() * fvc::interpolate(mixture.thermo().hSol(T, alpha1), "interpolate(h)"))))
                                     /fvc::interpolate(rho, "interpolate(rho)")
                                      ;
const surfaceScalarField hDeltaRhoGasSolidByRho = ((mixture.thermo().rhoGas() - mixture.thermo().rhoSolid())*fvc::interpolate(h, "interpolate(h)")
                                      - ((mixture.thermo().rhoGas() * fvc::interpolate(mixture.thermo().hGas(T, alpha1), "interpolate(h)"))
                                        -(mixture.thermo().rhoSolid() * fvc::interpolate(mixture.thermo().hSol(T, alpha1), "interpolate(h)"))))
                                     /fvc::interpolate(rho, "interpolate(rho)")
                                      ;

// const surfaceScalarField liquidFractionPrimeEnthalpyf =
//     fvc::interpolate(mixture.liquidFractionPrimeEnthalpy());

energySource =
(
  - fvc::laplacian(kappaByCp*hDeltaRhoSolidLiquidByRho, liquidFraction)
  // + fvc::laplacian(kappaByCp*mixture.HsPrimeAlphaGf(), alpha1)
  - fvc::laplacian(kappaByCp*hDeltaRhoGasSolidByRho, alpha1)
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
    fvm::laplacian(kappaByCp, h)
//  - fvm::laplacian(Hfus*kappaByCp*liquidFractionPrimeEnthalpyf, h)
//  + fvc::laplacian(Hfus*kappaByCp*liquidFractionPrimeEnthalpyf, h)
  + rho*energySource
  + redistribution*magGradAlphaM
   *(evaporativeCoolingPrime*h - fvm::Sp(evaporativeCoolingPrime, h))
  //+ fvc::div(rho*U & turbulence->devReff())
);
hEqn.solve();
pp(h);

if (gMin(h) < 0)
{
    runTime.writeNow();
    returnReduceOr(true);   // MPI barrier
    FatalError
        << "Enthalpy becomes negative!"
        << exit(FatalError);
}

mixture.correctThermo();
