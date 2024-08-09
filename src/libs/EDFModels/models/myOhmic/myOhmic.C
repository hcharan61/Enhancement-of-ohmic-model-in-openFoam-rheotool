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

#include "myOhmic.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EDFEquations
{
    defineTypeNameAndDebug(myOhmic, 0);
    addToRunTimeSelectionTable(EDFEquation, myOhmic, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::myOhmic::OhSpecie::OhSpecie
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    zi_(dict.lookup("z")),
    Di_(dict.lookup("D")),
    c0i_(dict.lookup("c0")),
    /*ci_
 (
   IOobject
   (
     name,
     phi.time().timeName(),
     phi.mesh(),
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
   ),
   phi.mesh()
 ),*/
    namei_(name) 
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::myOhmic::myOhmic
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    solvePhiE_(checkForPhiE(name, phi)),
    sigma_
    (
        IOobject
        (
            "sigma" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    psi_
    (
        IOobject
        (
            "psi" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    phiE_
    (
        IOobject
        (
            "phiE" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::READ_IF_PRESENT,
            solvePhiE_ == false ? (IOobject::NO_WRITE) : (IOobject::AUTO_WRITE)
        ),
        phi.mesh(),
        dimensionedScalar
        (
                "zero",
                psi_.dimensions(),
                pTraits<scalar>::zero
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    ),
    relPerm_(dict.lookup("relPerm")),
    Deff_("0", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.),
    T_(dict.lookup("T")),
    
    extraE_(dict.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero))),
    psiContrib_(dict.lookupOrDefault<bool>("psiContrib", true)),
    sigmaEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("sigmaEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    phiEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    
    maxIterSigma_(phi.mesh().solutionDict().subDict("electricControls").subDict("sigmaEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterPsi_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<int>("maxIter", 50)),
    
    species_(),
    nSpecies_(0)
{
    PtrList<entry> specEntries(dict.lookup("species"));
    nSpecies_ = specEntries.size();
    species_.setSize(nSpecies_);
    Info << nl << endl;
    dimensionedScalar cum("0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0);
      
    forAll (species_, specI)
    {
    
        species_.set
        (
            specI,
            new OhSpecie
            (
                specEntries[specI].keyword(),
                phi,
                specEntries[specI].dict()
            )
        );
        cum += species_[specI].DebyeLengthP(relPerm_, T_);   
    }   
    dimensionedScalar DebL
     (
       sqrt
           (
              Foam::EDFEquation::epsilonK_* relPerm_ * Foam::EDFEquation::kbK_*T_     
            / (cum*Foam::EDFEquation::FK_*Foam::EDFEquation::eK_)
           )
     );
     
    Info << "Debye length: " << DebL.value() << " m." << nl << endl;
    // Compute Deff_ (assuming two species only)
    scalar z0 = mag(species_[0].zi().value());
    scalar z1 = mag(species_[1].zi().value());
       
    Deff_ =  2*( species_[0].Di()*species_[1].Di() )
           / ( species_[0].Di() + species_[1].Di() ); 
    
    // This is a fix to avoid FP error when the electric module
    // is computed after hydrodynamics. This call ensures that on
    // first time-step the electric component is the 1st computed.
    // Issue appeared after reordering solving sequence (v3.0).     
    correct();  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::EDFEquations::myOhmic::Fe() const
{
    // Assume rhoE is either calculated within this context or passed from the PB model
    // This requires access to PB model variables like psi_, zi_, c0i_, and T_
    // If not directly accessible, calculate rhoE here or ensure it's updated before this method is called

    
        volScalarField rhoE( psi_ * dimensionedScalar("norm", epsilonK_.dimensions()/dimArea, 0.) );
    
        forAll (species_, i)
        {
      rhoE 
         += (
               species_[i].zi()*species_[i].c0i()*FK_
             * Foam::exp(-eK_*species_[i].zi()*psi_/(kbK_*T_) )
            );
         }// Placeholder for rhoE calculation logic from PB model if needed here
        // Otherwise, use the rhoE calculated in the PB model directly
     if (solvePhiE_)
     {
       if (psiContrib_)
        {
           return
           (
             (-rhoE + fvc::laplacian(epsilonK_*relPerm_, phiE_, "rhoE")) * ( fvc::grad(phiE_+psi_) - extraE_) 
           );
        }
       else
        {
           return
           (
             (-rhoE + fvc::laplacian(epsilonK_*relPerm_, phiE_, "rhoE")) * ( fvc::grad(phiE_) - extraE_) 
           ); 
        }   
     }
    else
     {
       if (psiContrib_)
        {
           return
           (
             -rhoE  * ( fvc::grad(psi_) - extraE_) 
           );
        }
       else
        {
           return
           (
             -rhoE  * (-extraE_) 
           ); 
        }  
     } 
    

    // Use the calculated rhoE to compute the body force as specified
     
}

void Foam::EDFEquations::myOhmic::correct()
{

       scalar res=GREAT; 
       scalar iter=0; 
         
  
   //- Equation for the conductivity 
   
       while (res > sigmaEqRes_ && iter < maxIterSigma_)
          { 
          
		fvScalarMatrix sigmaEqn
		(
		    fvm::ddt(sigma_)
		  + fvm::div(phi(), sigma_)
		  ==
		    fvm::laplacian(Deff_, sigma_, "laplacian(Deff,sigma)") 
		     
		);
		
		sigmaEqn.relax();
		res=sigmaEqn.solve().initialResidual();
		
		iter++;
          }
          
          //- Equation for the current (rhoE)       
 
  
  //- Equation for the current (rhoE)       
       res=GREAT;
       iter=0;  
       
       if (solvePhiE_)
       {
       while (res > phiEEqRes_ && iter < maxIterPhiE_)
         { 

		
		fvScalarMatrix phiEEqn
		(	  
		     fvm::laplacian(sigma_, phiE_ ) 		    
		);
	        
	        phiEEqn.relax();
		res=phiEEqn.solve().initialResidual();

		iter++;
        }
        }
        //- Equation for the intrinsic potential
        res=GREAT;
       iter=0;
   
       volScalarField souE(psi_ * dimensionedScalar("norm1",dimless/dimArea,0.));
       volScalarField souI(psi_ * dimensionedScalar("norm2",dimless/(dimArea*psi_.dimensions()),0.));
  
       forAll (species_, i)
       {
           
           dimensionedScalar bi(-eK_*species_[i].zi()/(kbK_*T_)); 
           volScalarField ai((species_[i].zi()*species_[i].c0i()*FK_) * Foam::exp(bi*psi_));
       
           souE += (
                     (-1./(relPerm_*epsilonK_)) * (ai - bi*ai*psi_)
                   );
                
           souI += (
                     (-1./(relPerm_*epsilonK_)) * (bi*ai)
                   );
       }
  
        while (res > psiEqRes_ && iter < maxIterPsi_)
        { 

		fvScalarMatrix psiEqn
		(	  
		     fvm::laplacian(psi_)
		  == 
	             fvm::SuSp(souI, psi_)
		  +  souE
		           		    
		);
		
	        psiEqn.relax();
		res=psiEqn.solve().initialResidual();

		iter++;
        } 
   
       
}

// ************************************************************************* //
