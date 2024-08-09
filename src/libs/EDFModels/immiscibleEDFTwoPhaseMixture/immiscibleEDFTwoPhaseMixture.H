/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::immiscibleConstitutiveTwoPhaseMixture

Description
    An immiscible GNF or VE two-phase mixture transport model

SourceFiles
    immiscibleConstitutiveTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleEDFTwoPhaseMixture_H
#define immiscibleEDFTwoPhaseMixture_H

#include "EDFTwoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class immiscibleConstitutiveTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleEDFTwoPhaseMixture
:
    public EDFTwoPhaseMixture,
    public interfaceProperties
{

public:

    // Constructors

        //- Construct from components
        immiscibleEDFTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~immiscibleEDFTwoPhaseMixture()
    {}


    // Member Functions
 
        //- Correct the transport and interface properties
        virtual void correctAll()
        {
            EDFTwoPhaseMixture::correct();
            interfaceProperties::correct();
        }
        
        //- Correct the transport properties
        virtual void correctFe()
        {
            EDFTwoPhaseMixture::correct();
        }
        
        //- Correct the interface properties
        virtual void correctInterface()
        {
            interfaceProperties::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read()
        {
            return
                EDFTwoPhaseMixture::read()
             && interfaceProperties::read();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
