scalar maxDeltaTemp
(
    runTime.controlDict().get<scalar>("maxDeltaTemp")
);

const scalar nCellsGlobal = mesh.globalData().nTotalCells();

const scalarField deltaT = (T - T.oldTime());
const scalar meanDeltaTempNum = gSum(deltaT)/nCellsGlobal;
const scalar maxDeltaTempNum = gMax(deltaT);

Info<< "Temperature change mean: " << meanDeltaTempNum
    << " max: " << maxDeltaTempNum << endl;

const scalar meanTemp = gSum(T)/nCellsGlobal;
const scalar maxTemp = gMax(T);

Info<< "Temperature mean: " << meanTemp
    << " max: " << maxTemp << endl;

const scalar meanVelocity = gSum(mag(U.primitiveField()))/nCellsGlobal;
const scalar maxVelocity = gMax(mag(U.primitiveField()));

Info<< "Velocity mean: " << meanVelocity
    << " max: " << maxVelocity << endl;

const scalar metalMass  = fvc::domainIntegrate(rho1 * alpha1).value();
Info<< "Mass in metal: " << metalMass << endl;
