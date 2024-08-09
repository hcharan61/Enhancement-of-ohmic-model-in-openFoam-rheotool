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

#ifndef myOhmic_H
#define myOhmic_H

#include "EDFEquation.H"

namespace Foam
{
namespace EDFEquations
{

class myOhmic
:
    public EDFEquation
{
    // Private data
    //- True if phiE_ exists in the directory, false otherwise  
        bool solvePhiE_;

        //- Externally applied electric potential
        volScalarField sigma_;
        //- Intrinsic electric potential
        volScalarField psi_;

        //- Internal electric potential
        volScalarField phiE_;

        //- Model constants

            // Electrical permittivity of the fluid
            dimensionedScalar relPerm_;
            
            // Effective diffusivity (m2/s)     
            dimensionedScalar Deff_;
            // Absolute temperature (K)  
            dimensionedScalar T_;
          
       //- Extra electric field (optional)
       dimensionedVector extraE_;
       bool psiContrib_;

       //- Solution controls for phiE and sigma
          
            // Residuals treshold
            scalar sigmaEqRes_;
            scalar phiEEqRes_;
            scalar psiEqRes_;
            
            // Max number of iterations
            int maxIterSigma_;
            int maxIterPhiE_;
            int maxIterPsi_;
            
       class OhSpecie
          {   
              
            //- Model parameters
            
              // Ionic valence (-)     
                dimensionedScalar zi_;
                
              // Diffusion coefficient (m2/s)     
                dimensionedScalar Di_;
                // Bulk concentration (mol/m3)   
                dimensionedScalar c0i_;
                //volScalarField ci_;
                
              // Species name     
                word namei_;  
          
           // Private Member Functions

              //- Disallow default bitwise copy construct
                OhSpecie(const OhSpecie&);

              //- Disallow default bitwise assignment
                void operator=(const OhSpecie&);
                
              public:
              
               // Constructor
               
                OhSpecie
                (
                    const word& name,
                    const surfaceScalarField& phi,
                    const dictionary& dict
                );
              
               // Destructor
               ~OhSpecie() {}
             
             // Public Member Functions 
             
             // Const-access to zi   
              const dimensionedScalar zi() const
               {
                   return zi_;
               } 
               
             // Const-access to Di   
              const dimensionedScalar Di() const
               {
                   return Di_;
               }
               // Const-access to c0i   
              const dimensionedScalar c0i() const
               {
                   return c0i_;
               }
              /* // Const access to ci   
              const volScalarField& ci() const
               {
                   return ci_;
               }
               
              // Non-const access to ci   
              volScalarField& ci()
               {
                   return ci_;
               }*/
               
             // Const access to namei   
              const word namei() const
               {
                   return namei_;
               }  
               // Const-access to member func used in 
             // computing the Debye length scale  
              const dimensionedScalar DebyeLengthP
               (
                 const dimensionedScalar relPerm,
                 const dimensionedScalar T
               ) const
               {
                   
                   return 
                   (    
                      zi_*zi_*c0i_        
                   );
                             
               }
       };
       
       //- List of species
       PtrList<OhSpecie> species_;
       
       //- Number of species
       int nSpecies_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        myOhmic(const myOhmic&);

        //- Disallow default bitwise assignment
        void operator=(const myOhmic&);
     
public:

    //- Runtime type information
    TypeName("myOhmic");

    // Constructors

        //- Construct from components
        myOhmic
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~myOhmic()
        {}


    // Member Functions

        //- Return the electric force for the momentum equation
         virtual tmp<volVectorField> Fe() const;

        //- Update the electric force
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace EDFEquations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
