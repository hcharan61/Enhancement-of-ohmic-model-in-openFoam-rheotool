/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "EDFTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

namespace Foam
{
    defineTypeNameAndDebug(EDFTwoPhaseMixture, 0);
}

Foam::EDFTwoPhaseMixture::EDFTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "electricProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),
    phase1_
    (
        new EDFModel(phi)
    ),
    phase2_
    (
        new EDFModel(phi)
    ),
    U_(U),
    phi_(phi)
{
    const dictionary& phase1Dict = subDict("water").subDict("parameters");
    phase1_->read(phase1Dict);

    const dictionary& phase2Dict = subDict("air").subDict("parameters");
    phase2_->read(phase2Dict);
}

Foam::tmp<Foam::volVectorField>
Foam::EDFTwoPhaseMixture::Fe() const
{
    const volScalarField bAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const volScalarField bAlpha2(1.0 - bAlpha1);

    tmp<volVectorField> FePhase1 = phase1_->Fe();
    tmp<volVectorField> FePhase2 = phase2_->Fe();

    tmp<volVectorField> FeTotal = bAlpha1 * FePhase1 + bAlpha2 * FePhase2;

    return FeTotal;
}

void Foam::EDFTwoPhaseMixture::correct()
{
    phase1_->correct();
    phase2_->correct();
}

bool Foam::EDFTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        const dictionary& phase1Dict = subDict("water").subDict("parameters");
        phase1_->read(phase1Dict);

        const dictionary& phase2Dict = subDict("air").subDict("parameters");
        phase2_->read(phase2Dict);

        return true;
    }
    else
    {
        return false;
    }
}

