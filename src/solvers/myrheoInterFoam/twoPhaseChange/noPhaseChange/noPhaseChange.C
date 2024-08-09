#include "noPhaseChange.H"
#include "fvScalarMatrix.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace twoPhaseChangeModels
{
    defineTypeNameAndDebug(noPhaseChange, 0);
    addToRunTimeSelectionTable(twoPhaseChangeModel, noPhaseChange, dictionary);
}
}

Foam::twoPhaseChangeModels::noPhaseChange::noPhaseChange
(
    const immiscibleConstitutiveTwoPhaseMixture& mixture,
    const immiscibleEDFTwoPhaseMixture& EDFmixture
)
:
    twoPhaseChangeModel(typeName, mixture, EDFmixture)
{}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::noPhaseChange::mDotAlphal() const
{
    return Pair<tmp<volScalarField>>(volScalarField::null(), volScalarField::null());
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::noPhaseChange::mDotP() const
{
    return Pair<tmp<volScalarField>>(volScalarField::null(), volScalarField::null());
}

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::noPhaseChange::Salpha
(
    volScalarField& alpha
) const
{
    return Pair<tmp<volScalarField::Internal>>(tmp<volScalarField::Internal>(nullptr), tmp<volScalarField::Internal>(nullptr));
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::noPhaseChange::Sp_rgh
(
    const volScalarField& rho,
    const volScalarField& gh,
    volScalarField& p_rgh
) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(p_rgh, dimVolume/dimTime));
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::noPhaseChange::SU
(
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimMass*dimVelocity/dimTime));
}

void Foam::twoPhaseChangeModels::noPhaseChange::correct()
{
    twoPhaseChangeModel::correct();
}

bool Foam::twoPhaseChangeModels::noPhaseChange::read()
{
    return twoPhaseChangeModel::read();
}

