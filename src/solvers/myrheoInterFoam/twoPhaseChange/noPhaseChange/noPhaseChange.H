#ifndef noPhaseChange_H
#define noPhaseChange_H

#include "twoPhaseChangeModel.H"

namespace Foam
{
namespace twoPhaseChangeModels
{

class noPhaseChange : public twoPhaseChangeModel
{
public:
    // Runtime type information
    TypeName("none");

    // Constructors
    noPhaseChange(const immiscibleConstitutiveTwoPhaseMixture& mixture, const immiscibleEDFTwoPhaseMixture& EDFmixture);

    // Destructor
    virtual ~noPhaseChange() {}

    // Member Functions
    Pair<tmp<volScalarField>> mDotAlphal() const override;
    Pair<tmp<volScalarField>> mDotP() const override;
    Pair<tmp<volScalarField::Internal>> Salpha(volScalarField& alpha) const override;
    tmp<fvScalarMatrix> Sp_rgh(const volScalarField& rho, const volScalarField& gh, volScalarField& p_rgh) const override;
    tmp<fvVectorMatrix> SU(const volScalarField& rho, const surfaceScalarField& rhoPhi, volVectorField& U) const override;
    void correct() override;
    bool read() override;
};

} // End namespace twoPhaseChangeModels
} // End namespace Foam

#endif

