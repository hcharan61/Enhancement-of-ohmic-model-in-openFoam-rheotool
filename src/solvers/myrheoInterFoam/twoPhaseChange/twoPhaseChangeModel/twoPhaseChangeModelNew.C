#include "twoPhaseChangeModel.H"
#include "noPhaseChange.H"
#include "immiscibleConstitutiveTwoPhaseMixture.H"
#include "immiscibleEDFTwoPhaseMixture.H"

Foam::autoPtr<Foam::twoPhaseChangeModel> Foam::twoPhaseChangeModel::New
(
    const immiscibleConstitutiveTwoPhaseMixture& mixture,
    const immiscibleEDFTwoPhaseMixture& EDFmixture
)
{
    IOobject twoPhaseChangeModelIO
    (
        IOobject
        (
            phaseChangePropertiesName,
            mixture.U().time().constant(),
            mixture.U().db(),
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
        Info << "No phase change: "
             << twoPhaseChangeModelIO.name()
             << " not found" << endl;
    }

    Info << "Selecting phaseChange model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << twoPhaseChangeModel::typeName << " type "
            << modelType << nl << nl
            << "Valid twoPhaseChangeModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseChangeModel>(cstrIter()(mixture, EDFmixture));
}

