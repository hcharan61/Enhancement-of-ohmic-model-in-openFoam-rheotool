#ifndef twoPhaseChangeModel_H
#define twoPhaseChangeModel_H

#include "immiscibleConstitutiveTwoPhaseMixture.H"
#include "immiscibleEDFTwoPhaseMixture.H"
#include "fvMatricesFwd.H"
#include "Pair.H"

namespace Foam
{

class twoPhaseChangeModel : public IOdictionary
{
    // Private Member Functions
    IOobject createIOobject(const immiscibleConstitutiveTwoPhaseMixture& mixture) const;

protected:
    const immiscibleConstitutiveTwoPhaseMixture& mixture_;
    const immiscibleEDFTwoPhaseMixture& EDFmixture_;
    dictionary twoPhaseChangeModelCoeffs_;

public:
    TypeName("phaseChangeModel");
    static const word phaseChangePropertiesName;

    declareRunTimeSelectionTable
    (
        autoPtr,
        twoPhaseChangeModel,
        dictionary,
        (
            const immiscibleConstitutiveTwoPhaseMixture& mixture,
            const immiscibleEDFTwoPhaseMixture& EDFmixture
        ),
        (mixture, EDFmixture)
    );

    twoPhaseChangeModel
    (
        const word& type,
        const immiscibleConstitutiveTwoPhaseMixture& mixture,
        const immiscibleEDFTwoPhaseMixture& EDFmixture
    );

    twoPhaseChangeModel(const twoPhaseChangeModel&) = delete;

    static autoPtr<twoPhaseChangeModel> New
    (
        const immiscibleConstitutiveTwoPhaseMixture& mixture,
        const immiscibleEDFTwoPhaseMixture& EDFmixture
    );

    virtual ~twoPhaseChangeModel() {}

    virtual Pair<tmp<volScalarField::Internal>> Salpha(volScalarField& alpha) const = 0;
    virtual tmp<fvScalarMatrix> Sp_rgh(const volScalarField& rho, const volScalarField& gh, volScalarField& p_rgh) const = 0;
    virtual tmp<fvVectorMatrix> SU(const volScalarField& rho, const surfaceScalarField& rhoPhi, volVectorField& U) const = 0;
    virtual void correct() = 0;
    virtual bool read() = 0;
    virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;
    virtual Pair<tmp<volScalarField>> mDotP() const = 0;
};

} // End namespace Foam

#endif

