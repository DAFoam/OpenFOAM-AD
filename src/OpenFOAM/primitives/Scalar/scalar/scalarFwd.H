/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

Description
    Typedefs for float/double/scalar without requiring scalar.H

    The definition of scalar as a floating-point number depends
    on pre-processor macros WM_SP, WM_SPDP or WM_DP being defined.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_scalarFwd_H
#define Foam_scalarFwd_H

// codi: add CoDiPack header
#include "codi.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- A typedef for float
typedef float floatScalar;

//- A typedef for double
// codi: comment out this
// typedef double doubleScalar;

// codi: replace doubleScalar with the codi's data type
#ifdef CODI_ADF
typedef codi::RealForward doubleScalar; // forward mode AD
#elif defined(CODI_ADR)
typedef codi::RealReverse doubleScalar; // reverse mode AD
#else
#error "WM_AD_MODE must be defined as ADF or ADR in etc/bashrc!"
#endif

#if defined(WM_SP)

typedef floatScalar scalar;
typedef floatScalar solveScalar;

#elif defined(WM_SPDP)

typedef floatScalar scalar;
typedef doubleScalar solveScalar;

#elif defined(WM_DP)

typedef doubleScalar scalar;
typedef doubleScalar solveScalar;

#else
//    #error "PRECISION must be set to WM_SP, WM_SPDP or WM_DP"
#endif

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
