/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "Colebrook.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::dimensionedScalar Foam::Colebrook::calcRe()
{
	return Uin_*rho_*Dh_/eta_;
}

Foam::dimensionedScalar Foam::Colebrook::calcCf()
{
	label nIters = 0;

	dimensionedScalar Cf("Cf", dimless, 0.01);
	Cf_.value() = Cf.value();
	do 
	{
		Cf.value() = Cf_.value();
		const scalar Re = calcRe().value();
		Cf_.value() = 1.0/pow(max(-A_*log10(ebyDh_/B_ + C_/(Re*sqrt(Cf.value()))) + D_, SMALL),2);
		nIters++;
	}
	while ((mag(Cf.value() - Cf_.value()) > conv_) && (nIters < maxIters_));

	// formula for Cf implemented is for Fanning friction
	// factor which is 4 times lower than Darcy friction factor
	Info<< endl;
	Info<< "Re = " << calcRe() << endl;
	Cf.value() = 4*Cf_.value(); 
	Info<< "Cf =  " << Cf_.value() << endl;
	Info<< endl;
		
	return Cf;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Colebrook::Colebrook()
:
	// Default are taken from https://doi.org/10.1016/j.cryogenics.2008.02.004
	A_{4.0},
	B_{3.7},
	C_{1.25},
	D_{0},
	ebyDh_{1.4e-4},
	// Default are for superfluid helium at p = 1 bar and T = 1.7 K
	rho_{dimensionedScalar("rho", dimDensity, 147)},
	Uin_{dimensionedScalar("Uin", dimVelocity, 16)},
	Dh_{dimensionedScalar("Dh", dimLength, 0.01)},
	eta_{dimensionedScalar("eta", dimDynamicViscosity, 0.0000013591)},
	Cf_{calcCf()}
{}


Foam::Colebrook::Colebrook
(
	const dimensionedScalar& rho,
	const dimensionedScalar& Uin,
	const dimensionedScalar& Dh,
	const dimensionedScalar& eta,
	const scalar ebyDh,
	const scalar A,
	const scalar B,
	const scalar C,
	const scalar D
)
:
	A_{A}, B_{B}, C_{C}, D_{D}, ebyDh_{ebyDh},
	rho_{rho}, Uin_{Uin}, Dh_{Dh}, eta_{eta},
	Cf_{calcCf()}
{ }


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void Foam::Colebrook::operator=(const Colebrook& rhs)
//{
//    if (this == &rhs)
//    {
//        return;  // Self-assignment is a no-op
//    }
//}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
