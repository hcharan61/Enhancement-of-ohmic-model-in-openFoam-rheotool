/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twoPhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseChangeModel, 0);
    defineRunTimeSelectionTable(twoPhaseChangeModel, dictionary);
}

const Foam::word Foam::twoPhaseChangeModel::phaseChangePropertiesName
(
    "phaseChangeProperties"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::twoPhaseChangeModel::createIOobject
(
    const immiscibleConstitutiveTwoPhaseMixture& constitutiveMixture,
    const immiscibleEDFTwoPhaseMixture& edfMixture
) const
{
    IOobject io
    (
        phaseChangePropertiesName,
        constitutiveMixture.U().mesh().time().constant(),
        constitutiveMixture.U().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModel::twoPhaseChangeModel
(
    const word& type,
    const immiscibleConstitutiveTwoPhaseMixture& constitutiveMixture,
    const immiscibleEDFTwoPhaseMixture& edfMixture
)
:
    IOdictionary(createIOobject(constitutiveMixture, edfMixture)),
    constitutiveMixture_(constitutiveMixture),
    edfMixture_(edfMixture),
    twoPhaseChangeModelCoeffs_(optionalSubDict(type + "Coeffs"))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseChangeModel::correct()
{
    // Correct both mixture models
    constitutiveMixture_.correct();
    edfMixture_.correct();
}

bool Foam::twoPhaseChangeModel::read()
{
    if (regIOobject::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::twoPhaseChangeModel> Foam::twoPhaseChangeModel::New
(
    const immiscibleConstitutiveTwoPhaseMixture& constitutiveMixture,
    const immiscibleEDFTwoPhaseMixture& edfMixture
)
{
    IOobject twoPhaseChangeModelIO
    (
        IOobject
        (
            phaseChangePropertiesName,
            constitutiveMixture.U().time().constant(),
            constitutiveMixture.U().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(twoPhaseChangeModels::noPhaseChange::typeName);

    if (twoPhaseChangeModelIO.typeHeaderOk<IOdictionary>(false))
    {
        IOdictionary(twoPhaseChangeModelIO).lookup
        (
            twoPhaseChangeModel::typeName
        ) >> modelType;
    }
    else
    {
        Info<< "No phase change: "
            << twoPhaseChangeModelIO.name()
            << " not found" << endl;
    }

    Info<< "Selecting phaseChange model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << twoPhaseChangeModel::typeName << " type "
            << modelType << nl << nl
            << "Valid twoPhaseChangeModels are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseChangeModel>(cstrIter()(constitutiveMixture, edfMixture));
}

// ************************************************************************* //

