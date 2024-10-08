/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "interpolationCellPointBary.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationCellPointBary<Type>::interpolationCellPointBary
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    fieldInterpolation<Type, interpolationCellPointBary<Type>>(psi),
    volPointInterpolation_(psi.mesh()),
    psip_
    (
        volPointInterpolation_.interpolate
        (
            psi,
            "volPointInterpolate(" + psi.name() + ')',
            true        // use cache
        )
    )
{
    // Uses cellPointWeight to do interpolation which needs tet decomposition
    (void)psi.mesh().tetBasePtIs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolationCellPointBary<Type>::updatePsip
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    psip_ = volPointInterpolation_.interpolate
     (
       psi,
       "volPointInterpolate(" + psi.name() + ')',
       true        // use cache
     );
}

// ************************************************************************* //
