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

Class
    Foam::Colebrook

Description
	The class for calculation of Darcy friction factor Cf for a rough wall.
	Currently it computes Cf from Colebrook correlation for a rough wall.
	In future, it can be extended to some an abstract class
	and implement specific derived clases, e.g. Colebrook.

	Solved equation is of the form:

	1.0/sqrt(Cf) = -A*log( (e/Dh)/B + C/(Re*sqrt(Cf)) ) + D

	where 
	A, B, C, D are equation constant, 
	e is a wall roughness,
	Dh is a hydraulic diameter,
	Re is the Reynolds number.

SourceFiles
    Colebrook.C

\*---------------------------------------------------------------------------*/

#ifndef Colebrook_H
#define Colebrook_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
//class Istream;
//class Ostream;
//
////Istream& operator>>(Istream&, Colebrook&);
//Ostream& operator<<(Ostream&, const Colebrook&);


/*---------------------------------------------------------------------------*\
                         Class Colebrook Declaration
\*---------------------------------------------------------------------------*/

class Colebrook
//:
//    public baseClassName
{
private:
    // Private Data

		//- A constant
		scalar A_;

		//- B constant
		scalar B_;

		//- C constant
		scalar C_;

		//- D constant
		scalar D_;

		//- Relative roughness
		scalar ebyDh_;

		//- Maximum number of iterations
		const label maxIters_ = 1000;

		//- Convergence criterion
		const scalar conv_ = 1e-12;

        //- Density
        dimensionedScalar rho_;

        //- Specific velocity
        dimensionedScalar Uin_;

        //- Hydraulic diameter
        dimensionedScalar Dh_;

        //- Dynamic viscosity
        dimensionedScalar eta_;

        //- Fanning friction factor
        dimensionedScalar Cf_;

    // Private Member Functions

        //- No copy construct
        //Colebrook(const Colebrook&) = delete;

        //- No copy assignment
        void operator=(const Colebrook&) = delete;

		//- Calculates Reynolds number
		dimensionedScalar calcRe();

		//- Calculates Fanning friction factor Cf
		dimensionedScalar calcCf();


public:

    // Static Data Members


    // Generated Methods

//        //- Default construct
//        Colebrook() = default;
//
//        //- Copy construct
//        Colebrook(const Colebrook&) = default;
//
//        //- Copy assignment
//        Colebrook& operator=(const Colebrook&) = default;


    // Constructors

        //- Default construct
        Colebrook();

        //- Construct from components
        Colebrook
		(
			const dimensionedScalar& rho,
			const dimensionedScalar& Uin,
			const dimensionedScalar& Dh,
			const dimensionedScalar& eta,
			const scalar ebyDh = 1.4e-4,
			const scalar A = 4.0,
			const scalar B = 3.7,
			const scalar C = 1.25,
			const scalar D = 0 
		);

        //- Copy construct
        //Colebrook(const Colebrook&);

        //- Construct from Istream
        //explicit Colebrook(Istream& is);


    // Selectors


    //- Destructor
    ~Colebrook() = default;


    // Member Functions
		

    // Access
		
		//- Returns Darcy friction factor
		const dimensionedScalar& Cf() const
		{
			return Cf_;
		}

    // Check

    // Edit

    // Write


    // Member Operators

        //void operator=(const Colebrook&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

//        friend Istream& operator>>(Istream&, Colebrook&);
        friend Ostream& operator<<(Ostream& os, const Colebrook& cb)
		{
			return os << "Darcy friction factor is Cf = " << cb.Cf() << endl;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "ColebrookI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
