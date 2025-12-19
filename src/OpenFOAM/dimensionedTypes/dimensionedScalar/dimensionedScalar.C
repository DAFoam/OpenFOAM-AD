/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// codi: Fix problematic (ambiguous) scalar operators. They conflict with codipack's operators
// here we restrict the input type to dimensionedScalar and scalar.
template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, scalar>::value,
    dimensionedScalar>::type
operator+(const T1& ds1, const T2 s2)
{
    return ds1 + dimensionedScalar(s2);
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, scalar>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator+(const T1 s1, const T2& ds2)
{
    return dimensionedScalar(s1) + ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, scalar>::value,
    dimensionedScalar>::type
operator-(const T1& ds1, const T2 s2)
{
    return ds1 - dimensionedScalar(s2);
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, scalar>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator-(const T1 s1, const T2& ds2)
{
    return dimensionedScalar(s1) - ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, scalar>::value,
    dimensionedScalar>::type
operator*(const T1& ds1, const T2 s2)
{
    return ds1 * dimensionedScalar(s2);
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, scalar>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator*(const T1 s1, const T2& ds2)
{
    return dimensionedScalar(s1) * ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, scalar>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator/(const T1 s1, const T2& ds2)
{
    return dimensionedScalar(s1)/ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, scalar>::value,
    dimensionedScalar>::type
operator/(const T1& ds1, const T2 s2)
{
    return ds1 / dimensionedScalar(s2);
}

// codi: here we need to define operators that between double and dimensionedScalar
template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, double>::value,
    dimensionedScalar>::type
operator+(const T1& ds1, const T2 d2)
{
    return ds1 + dimensionedScalar(scalar(d2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, double>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator+(const T1 d1, const T2& ds2)
{
    return dimensionedScalar(scalar(d1)) + ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, double>::value,
    dimensionedScalar>::type
operator-(const T1& ds1, const T2 d2)
{
    return ds1 - dimensionedScalar(scalar(d2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, double>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator-(const T1 d1, const T2& ds2)
{
    return dimensionedScalar(scalar(d1)) - ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, double>::value,
    dimensionedScalar>::type
operator*(const T1& ds1, const T2 d2)
{
    return ds1 * dimensionedScalar(scalar(d2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, double>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator*(const T1 d1, const T2& ds2)
{
    return dimensionedScalar(scalar(d1)) * ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, double>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator/(const T1 d1, const T2& ds2)
{
    return dimensionedScalar(scalar(d1)) / ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, double>::value,
    dimensionedScalar>::type
operator/(const T1& ds1, const T2 d2)
{
    return ds1 / dimensionedScalar(scalar(d2));
}

// codi: here we need to define operators that between int and dimensionedScalar
template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, int>::value,
    dimensionedScalar>::type
operator+(const T1& ds1, const T2 i2)
{
    return ds1 + dimensionedScalar(scalar(i2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, int>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator+(const T1 i1, const T2& ds2)
{
    return dimensionedScalar(scalar(i1)) + ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, int>::value,
    dimensionedScalar>::type
operator-(const T1& ds1, const T2 i2)
{
    return ds1 - dimensionedScalar(scalar(i2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, int>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator-(const T1 i1, const T2& ds2)
{
    return dimensionedScalar(scalar(i1)) - ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, int>::value,
    dimensionedScalar>::type
operator*(const T1& ds1, const T2 i2)
{
    return ds1 * dimensionedScalar(scalar(i2));
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, int>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator*(const T1 i1, const T2& ds2)
{
    return dimensionedScalar(scalar(i1)) * ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, int>::value &&
    std::is_same<typename std::decay<T2>::type, dimensionedScalar>::value,
    dimensionedScalar>::type
operator/(const T1 i1, const T2& ds2)
{
    return dimensionedScalar(scalar(i1)) / ds2;
}

template<typename T1, typename T2>
typename std::enable_if<
    std::is_same<typename std::decay<T1>::type, dimensionedScalar>::value &&
    std::is_same<typename std::decay<T2>::type, int>::value,
    dimensionedScalar>::type
operator/(const T1& ds1, const T2 i2)
{
    return ds1 / dimensionedScalar(scalar(i2));
}

// codi: Explicit instantiations for dimensionedScalar operators with AD scalar types
// These instantiations fix linking errors for operator* and operator/ between
// scalar (which aliases to AD type) and dimensioned<scalar> (which aliases to dimensioned<AD>)
template dimensionedScalar operator+(const scalar, const dimensionedScalar&);
template dimensionedScalar operator-(const scalar, const dimensionedScalar&);
template dimensionedScalar operator*(const scalar, const dimensionedScalar&);
template dimensionedScalar operator/(const scalar, const dimensionedScalar&);

template dimensionedScalar operator+(const dimensionedScalar&, const scalar);
template dimensionedScalar operator-(const dimensionedScalar&, const scalar);
template dimensionedScalar operator*(const dimensionedScalar&, const scalar);
template dimensionedScalar operator/(const dimensionedScalar&, const scalar);

template dimensionedScalar operator+(const double, const dimensionedScalar&);
template dimensionedScalar operator-(const double, const dimensionedScalar&);
template dimensionedScalar operator*(const double, const dimensionedScalar&);
template dimensionedScalar operator/(const double, const dimensionedScalar&);

template dimensionedScalar operator+(const dimensionedScalar&, const double);
template dimensionedScalar operator-(const dimensionedScalar&, const double);
template dimensionedScalar operator*(const dimensionedScalar&, const double);
template dimensionedScalar operator/(const dimensionedScalar&, const double);

template dimensionedScalar operator+(const int, const dimensionedScalar&);
template dimensionedScalar operator-(const int, const dimensionedScalar&);
template dimensionedScalar operator*(const int, const dimensionedScalar&);
template dimensionedScalar operator/(const int, const dimensionedScalar&);

template dimensionedScalar operator+(const dimensionedScalar&, const int);
template dimensionedScalar operator-(const dimensionedScalar&, const int);
template dimensionedScalar operator*(const dimensionedScalar&, const int);
template dimensionedScalar operator/(const dimensionedScalar&, const int);

/*
dimensionedScalar operator+(const dimensionedScalar& ds1, const scalar s2)
{
    return ds1 + dimensionedScalar(s2);
}


dimensionedScalar operator+(const scalar s1, const dimensionedScalar& ds2)
{
    return dimensionedScalar(s1) + ds2;
}


dimensionedScalar operator-(const dimensionedScalar& ds1, const scalar s2)
{
    return ds1 - dimensionedScalar(s2);
}


dimensionedScalar operator-(const scalar s1, const dimensionedScalar& ds2)
{
    return dimensionedScalar(s1) - ds2;
}


dimensionedScalar operator*(const dimensionedScalar& ds1, const scalar s2)
{
    return ds1 * dimensionedScalar(s2);
}


dimensionedScalar operator/(const scalar s1, const dimensionedScalar& ds2)
{
    return dimensionedScalar(s1)/ds2;
}
*/

dimensionedScalar pow
(
    const dimensionedScalar& ds,
    const dimensionedScalar& expt
)
{
    return dimensionedScalar
    (
        "pow(" + ds.name() + ',' + expt.name() + ')',
        pow(ds.dimensions(), expt),
        pow(ds.value(), expt.value())
    );
}

// codi:
dimensionedScalar pow
(
    const dimensionedScalar& ds,
    const double& expt
)
{
    return dimensionedScalar
    (
        "pow(" + ds.name() + ',' + Foam::name(scalar(expt)) + ')',
        pow(ds.dimensions(), dimensionedScalar(scalar(expt))),
        pow(ds.value(), scalar(expt))
    );
}


dimensionedScalar pow3(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pow3(" + ds.name() + ')',
        pow3(ds.dimensions()),
        pow3(ds.value())
    );
}


dimensionedScalar pow4(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pow4(" + ds.name() + ')',
        pow4(ds.dimensions()),
        pow4(ds.value())
    );
}


dimensionedScalar pow5(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pow5(" + ds.name() + ')',
        pow5(ds.dimensions()),
        pow5(ds.value())
    );
}


dimensionedScalar pow6(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pow6(" + ds.name() + ')',
        pow6(ds.dimensions()),
        pow6(ds.value())
    );
}


dimensionedScalar pow025(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pow025(" + ds.name() + ')',
        pow025(ds.dimensions()),
        pow025(ds.value())
    );
}


dimensionedScalar sqrt(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "sqrt(" + ds.name() + ')',
        pow(ds.dimensions(), dimensionedScalar("0.5", dimless, 0.5)),
        sqrt(ds.value())
    );
}


dimensionedScalar cbrt(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "cbrt(" + ds.name() + ')',
        pow(ds.dimensions(), dimensionedScalar("(1|3)", dimless, 1.0/3.0)),
        cbrt(ds.value())
    );
}


dimensionedScalar sign(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "sign(" + ds.name() + ')',
        sign(ds.dimensions()),
        ::Foam::sign(ds.value())
    );
}


dimensionedScalar pos(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pos(" + ds.name() + ')',
        pos(ds.dimensions()),
        ::Foam::pos(ds.value())
    );
}


dimensionedScalar pos0(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "pos0(" + ds.name() + ')',
        pos0(ds.dimensions()),
        ::Foam::pos0(ds.value())
    );
}


dimensionedScalar neg(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "neg(" + ds.name() + ')',
        neg(ds.dimensions()),
        ::Foam::neg(ds.value())
    );
}


dimensionedScalar neg0(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "neg0(" + ds.name() + ')',
        neg0(ds.dimensions()),
        ::Foam::neg0(ds.value())
    );
}


dimensionedScalar posPart(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "posPart(" + ds.name() + ')',
        posPart(ds.dimensions()),
        ::Foam::pos0(ds.value())
    );
}


dimensionedScalar negPart(const dimensionedScalar& ds)
{
    return dimensionedScalar
    (
        "negPart(" + ds.name() + ')',
        negPart(ds.dimensions()),
        ::Foam::neg(ds.value())
    );
}


#define transFunc(func)                                                        \
dimensionedScalar func(const dimensionedScalar& ds)                            \
{                                                                              \
    if (dimensionSet::checking() && !ds.dimensions().dimensionless())          \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "scalar is not dimensionless: " << ds.dimensions() << nl        \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    return dimensionedScalar                                                   \
    (                                                                          \
        #func "(" + ds.name() + ')',                                           \
        dimless,                                                               \
        func(ds.value())                                                     \
    );                                                                         \
}

transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)
transFunc(erf)
transFunc(erfc)

// codi: comment out these functions
/*
transFunc(lgamma)
transFunc(j0)
transFunc(j1)
transFunc(y0)
transFunc(y1)
*/

#undef transFunc


#define transFunc(func)                                                        \
dimensionedScalar func(const int n, const dimensionedScalar& ds)               \
{                                                                              \
    if (dimensionSet::checking() && !ds.dimensions().dimensionless())          \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "scalar is not dimensionless: " << ds.dimensions() << nl        \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    return dimensionedScalar                                                   \
    (                                                                          \
        #func "(" + name(n) + ',' + ds.name() + ')',                           \
        dimless,                                                               \
        func(n, ds.value())                                                  \
    );                                                                         \
}

// codi: comment out these functions
//transFunc(jn)
//transFunc(yn)

#undef transFunc


dimensionedScalar atan2
(
    const dimensionedScalar& x,
    const dimensionedScalar& y
)
{
    return dimensionedScalar
    (
        "atan2(" + x.name() + ',' + y.name() + ')',
        atan2(x.dimensions(), y.dimensions()),
        atan2(x.value(), y.value())
    );
}


dimensionedScalar hypot
(
    const dimensionedScalar& x,
    const dimensionedScalar& y
)
{
    return dimensionedScalar
    (
        "hypot(" + x.name() + ',' + y.name() + ')',
        hypot(x.dimensions(), y.dimensions()),
        hypot(x.value(), y.value())
    );
}


dimensionedScalar stabilise
(
    const dimensionedScalar& x,
    const dimensionedScalar& y
)
{
    return dimensionedScalar
    (
        "stabilise(" + x.name() + ',' + y.name() + ')',
        stabilise(x.dimensions(), y.dimensions()),
        stabilise(x.value(), y.value())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
