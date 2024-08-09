/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::EDFTwoPhaseMixture

Description
    A two-phase mixture model handling the electric force term calculation

SourceFiles
    EDFTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef EDFTwoPhaseMixture_H
#define EDFTwoPhaseMixture_H

#include "EDFModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class EDFTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class EDFTwoPhaseMixture
:
    public IOdictionary,
    public twoPhaseMixture
{
protected:

    // Protected data

    autoPtr<EDFModel> phase1_;
    autoPtr<EDFModel> phase2_;

    const volVectorField& U_;
    const surfaceScalarField& phi_;

public:

    TypeName("EDFTwoPhaseMixture");

    // Constructors

    //- Construct from components
    EDFTwoPhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    //- Destructor
    virtual ~EDFTwoPhaseMixture()
    {}

    // Member Functions

    //- Return const-access to phase1 EDFModel
    const EDFModel& phase1() const
    {
        return phase1_();
    }

    //- Return const-access to phase2 EDFModel
    const EDFModel& phase2() const
    {
        return phase2_();
    }

    //- Return the electric force term for the momentum equation
    tmp<volVectorField> Fe() const;

    //- Correct the models
    virtual void correct();

    //- Read base transportProperties dictionary
    virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

