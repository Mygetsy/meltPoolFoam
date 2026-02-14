wasMelted = Foam::max(wasMelted, liquidFraction);

mixture.correctPassiveFields();

if (writeAllFields)
{
    heatConvection = fvc::div(rhoPhi, h);

    // For debug: check that heatConduction is almost equal to heatConduction2
    heatConduction = fvc::laplacian(mixture.kappaf(), T);

    #include "gradT.H"

    heatConduction2 = fvc::laplacian(rhof*kappaByRhoCp, h)
                     - fvc::laplacian(kappaByRhoCp*hDeltaRhoSolidLiquid, liquidFraction)
                     - fvc::laplacian(kappaByRhoCp*hDeltaRhoGasSolid, alpha1);
}


