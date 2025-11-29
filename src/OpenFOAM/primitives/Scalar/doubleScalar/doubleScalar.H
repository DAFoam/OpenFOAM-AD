/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Typedef
    Foam::doubleScalar

Description
    Floating-point double precision scalar type.

SourceFiles
    doubleScalar.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_double_H
#define Foam_primitives_double_H

#include "scalarFwd.H"
#include "doubleFloat.H"
#include "direction.H"
#include "pTraits.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Typedef (doubleScalar) in scalarFwd.H

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Largest and smallest scalar values allowed in certain parts of the code.
// See std::numeric_limits max(), min(), epsilon()
// codi: replaced doubleScalar with double
constexpr double doubleScalarGREAT = 1.0e+15;
constexpr double doubleScalarROOTGREAT = 3.0e+8;
constexpr double doubleScalarVGREAT = 1.0e+300;
constexpr double doubleScalarROOTVGREAT = 1.0e+150;
constexpr double doubleScalarSMALL = 1.0e-15;
constexpr double doubleScalarROOTSMALL = 3.0e-8;
constexpr double doubleScalarVSMALL = 1.0e-300;
constexpr double doubleScalarROOTVSMALL = 1.0e-150;


#define Scalar doubleScalar
#define ScalarVGREAT doubleScalarVGREAT
#define ScalarVSMALL doubleScalarVSMALL
#define ScalarROOTVGREAT doubleScalarROOTVGREAT
#define ScalarROOTVSMALL doubleScalarROOTVSMALL
#define ScalarRead readDouble


inline Scalar mag(const Scalar s)
{
    // codi: using ADL to auto-select func namespace
    using std::fabs;
    return fabs(s);
}

inline Scalar hypot(const Scalar x, const Scalar y)
{
    // codi: manual implementation of hypot
    return sqrt(x * x + y * y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    // codi: using ADL to auto-select func namespace
    using std::atan2;
    return atan2(y, x);
}

// codi: AD-aware min for scalar / doubleScalar
inline Scalar min(const Scalar& a, const Scalar& b)
{
    return (a < b) ? a : b;
}

// codi: AD-aware max for scalar / doubleScalar
inline Scalar max(const Scalar& a, const Scalar& b)
{
    return (a > b) ? a : b;
}

// Normal (double-precision) transcendental functions
// codi: using ADL to auto-select func namespace
#define transFunc(func)                                         \
inline Scalar func(const Scalar s)                              \
{                                                               \
    using std::func;                                            \
    return func(s);                                             \
}

// Normal (double-precision) bessel functions.
// May not be available on all systems
#ifdef Foam_no_besselFunc
    // Not available
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        std::cerr<< "No '" << #func << "' function\n";          \
        return 0;                                               \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        std::cerr<< "No '" << #func << "' function\n";          \
        return 0;                                               \
    }
#elif defined(__MINGW32__)
    // Mingw: with '_' prefix
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        return _##func(s);                                      \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        return _##func(n, s);                                   \
    }
#else
    #define besselFunc(func)                                    \
    inline Scalar func(const Scalar s)                          \
    {                                                           \
        return ::func(s.getValue());                            \
    }
    #define besselFunc2(func)                                   \
    inline Scalar func(const int n, const Scalar s)             \
    {                                                           \
        return ::func(n, s.getValue());                         \
    }
#endif
// codi: using passive model for besselFunc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// scalarImpl.H is used for template-like substitution (but using macros)

#define Foam_use_scalarImpl_header
#include "scalarImpl.H"
#undef Foam_use_scalarImpl_header

#undef Scalar
#undef ScalarVGREAT
#undef ScalarVSMALL
#undef ScalarROOTVGREAT
#undef ScalarROOTVSMALL
#undef ScalarRead
#undef transFunc
#undef besselFunc
#undef besselFunc2


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
