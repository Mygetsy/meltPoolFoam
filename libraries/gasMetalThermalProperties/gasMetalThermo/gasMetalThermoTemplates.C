/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of gasMetalThermalProperties.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "generateGeometricField.H"
#include "geometricUniformField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::scalar gasMetalAverage
(
    const Foam::Polynomial<2>& gasProperty,
    const Foam::Polynomial<2>& liquidProperty,
    const Foam::Polynomial<2>& solidProperty,
    Foam::scalar T,
    Foam::scalar liquidFraction,
    Foam::scalar gasFraction
)
{
    Foam::scalar solidFraction = 1 - liquidFraction - gasFraction;

    return
        solidProperty.value(T)*solidFraction
      + liquidProperty.value(T)*liquidFraction
      + gasProperty.value(T)*gasFraction;
}


} // End namespace

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::Cp
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    const T2 sharpLiquidFraction
    (
        (1 - gasFraction)*pos(liquidFraction - (1 - gasFraction)/2) //!TODO: Mb delete
    );

    return generateGeometricField<T1>
    (
        "Cp",
        mesh_,
        dimGasConstant,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar rho = alphaG*rhoGas_ + phi*rhoLiquid_*(1 + betaLiquid_*(T - Tmelting_)) + (alphaM - phi)*rhoSolid_;
            
            return (alphaG*rhoGas_*gas_.Cp.value(T)
                    + (alphaM - phi)*rhoSolid_*solid_.Cp.value(T)
                    + phi*rhoLiquid_*(1 + betaLiquid_*(T - Tmelting_))*liquid_.Cp.value(T))/rho;
        },
        T, sharpLiquidFraction, gasFraction
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::kappa
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    return generateGeometricField<T1>
    (
        "kappa",
        mesh_,
        dimPower/dimLength/dimTemperature,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;

            return alphaG*gas_.kappa.value(T)
                   + alphaM*( T <= Tmelting_ ? solid_.kappa.value(T) : liquid_.kappa.value(T));
        },
        T, liquidFraction, gasFraction
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T3> Foam::gasMetalThermo::h
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction,
    const word& name
) const
{
    return generateGeometricField<T3>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar rho = alphaG*rhoGas_ + phi*rhoLiquid_*(1 + betaLiquid_*(T - Tmelting_)) + (alphaM - phi)*rhoSolid_;
            
            return (alphaG*rhoGas_*gas_.Cp.integral(0, T)
                    + (alphaM - phi)*rhoSolid_*solid_.Cp.integral(0, T) 
                    + phi*rhoLiquid_*(1 + betaLiquid_*(T - Tmelting_))*(
                    solid_.Cp.integral(0, Tmelting_) + liquid_.Cp.integral(Tmelting_, T) + Hfusion_)
                    )/rho;
        },
        T, liquidFraction, gasFraction
    );
}

template<class T1, class T2>
Foam::tmp<T2> Foam::gasMetalThermo::hGas
(
    const T1& T,
    const T2& gasFraction, //!TODO: How to do not use this field.
    const word& name
) const
{
    return generateGeometricField<T2>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T)
        {
            return gas_.Cp.integral(0, T);
        },
        T
    );
}

template<class T1, class T2>
Foam::tmp<T2> Foam::gasMetalThermo::hSol
(
    const T1& T,
    const T2& gasFraction,//!
    const word& name
) const
{
    return generateGeometricField<T2>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T)
        {
            return solid_.Cp.integral(0, T);
        },
        T
    );
}

template<class T1, class T2>
Foam::tmp<T2> Foam::gasMetalThermo::hLiq
(
    const T1& T,
    const T2& gasFraction,//!
    const word& name
) const
{
    return generateGeometricField<T2>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T)
        {
            return solid_.Cp.integral(0, Tmelting_) + liquid_.Cp.integral(Tmelting_, T) + Hfusion_;
        },
        T
    );
}

template<class T1>
Foam::tmp<T1> Foam::gasMetalThermo::hAtMelting
(
    const T1& gasFraction
) const
{
    return h
    (
        geometricUniformField<scalar>(Tmelting_),
        ((1 - gasFraction)/2)(),
        gasFraction,
        "hAtMelting"
    );
}


template<class T1>
Foam::tmp<T1> Foam::gasMetalThermo::HsPrimeAlphaG
(
    const T1& T
) const
{
    return generateGeometricField<T1>
    (
        "HsPrimeAlphaG",
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T)
        {
            scalar piecewise =
                T <= Tmelting_
              ? solid_.Cp.integral(Tmelting_, T)
              : liquid_.Cp.integral(Tmelting_, T);

            return gas_.Cp.integral(0, T) - solid_.Cp.integral(0, Tmelting_) - piecewise;
        },
        T
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::T
(
    const T1& h,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{   
    return generateGeometricField<T1>
    (
        "T",
        mesh_,
        dimTemperature,
        [this](scalar h, scalar phi, scalar alphaG)
        {

            scalar alphaM = 1 - alphaG;
            scalar hLiquidTmelt = solid_.Cp.integral(0, Tmelting_) + Hfusion_;

            scalar A = (alphaM - phi)*rhoSolid_*solid_.Cp.derivative(0) 
                        + 2*betaLiquid_*phi*rhoLiquid_*liquid_.Cp.value(0);

            scalar B = betaLiquid_*phi*rhoLiquid_*(hLiquidTmelt - h)
                        + (1 - 2*Tmelting_*betaLiquid_)*phi*rhoLiquid_*liquid_.Cp.value(0)
                        + (alphaM - phi)*rhoSolid_*solid_.Cp.value(0)
                        + alphaG*rhoGas_*gas_.Cp.value(0);
            
            scalar C = hLiquidTmelt*phi*rhoLiquid_*(Tmelting_*betaLiquid_ - 1)
                        + Tmelting_*phi*rhoLiquid_*liquid_.Cp.value(0)*(1 - Tmelting_*betaLiquid_)
                        + (
                            (1 - Tmelting_*betaLiquid_)*phi*rhoLiquid_ + (alphaM - phi)*rhoSolid_ + alphaG*rhoGas_
                        )*h;

            return 
            mag(A) > ROOTSMALL
              ? (Foam::sqrt(sqr(B) + 2*A*C) - B)/A
              : C/B;
        },
        h, liquidFraction, gasFraction
    );
}

template<class T1, class T2>
Foam::tmp<T2> Foam::gasMetalThermo::phiCalc
(
    const T1& h,
    const T2& gasFraction
) const
{   
    return generateGeometricField<T2>
    (
        "phiCalc",
        mesh_,
        dimless,
        [this](scalar h, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar hGasTmelt = gas_.Cp.integral(0, Tmelting_);
            scalar hSolidTmelt = solid_.Cp.integral(0, Tmelting_);
            scalar hLiquidTmelt = solid_.Cp.integral(0, Tmelting_) + Hfusion_;

            scalar gasSolidEnthalpy = (alphaG*rhoGas_*hGasTmelt + alphaM*rhoSolid_*hSolidTmelt)
                                    /(alphaG*rhoGas_ + alphaM*rhoSolid_);

            scalar gasLiquidEnthalpy = (alphaG*rhoGas_*hGasTmelt + alphaM*rhoLiquid_*hLiquidTmelt)
                                    /(alphaG*rhoGas_ + alphaM*rhoLiquid_);
            
            scalar phiTmelt = (h*(alphaG*rhoGas_ + alphaM*rhoSolid_) - alphaG*rhoGas_*hGasTmelt - alphaM*rhoSolid_*hSolidTmelt)
                            /(h*(rhoSolid_ - rhoLiquid_) + rhoLiquid_*hLiquidTmelt - rhoSolid_*hSolidTmelt);
            

            if (h < gasSolidEnthalpy){
                return 0.;
                }
            else if ((h > gasLiquidEnthalpy)) {
                return alphaM;
                }
            else {
                return phiTmelt;
                }
        },
        h, gasFraction
    );
}


// ************************************************************************* //
