/* 
   codi: This is a new file to explicit instantiations for min/max with CoDiPack AD scalar.
   This avoids missing symbols when scalar aliases to an active type. 
*/

#include "Field.H"
#include "FieldField.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"

namespace Foam
{

#undef min
#undef max

// Template definitions for the template-template min/max functions
// This must be defined before instantiation
template<template<class> class FieldType, class Type>
void min
(
    FieldType<Type>& result,
    const FieldType<Type>& f1,
    const FieldType<Type>& f2
)
{
    const label n = result.size();

    for (label i = 0; i < n; ++i)
    {
        result[i] = Foam::min(f1[i], f2[i]);
    }
}

template<template<class> class FieldType, class Type>
void max
(
    FieldType<Type>& result,
    const FieldType<Type>& f1,
    const FieldType<Type>& f2
)
{
    const label n = result.size();

    for (label i = 0; i < n; ++i)
    {
        result[i] = Foam::max(f1[i], f2[i]);
    }
}

// Macro to explicitly instantiate min/max for a given type
#define INSTANTIATE_MIN_MAX(Type) \
    \
    /* FieldField on fvPatchField */ \
    template void min<fvPatchField, Type> \
    ( \
        FieldField<fvPatchField, Type>&, \
        const FieldField<fvPatchField, Type>&, \
        const FieldField<fvPatchField, Type>& \
    ); \
    template void max<fvPatchField, Type> \
    ( \
        FieldField<fvPatchField, Type>&, \
        const FieldField<fvPatchField, Type>&, \
        const FieldField<fvPatchField, Type>& \
    ); \
    \
    /* FieldField on fvsPatchField */ \
    template void min<fvsPatchField, Type> \
    ( \
        FieldField<fvsPatchField, Type>&, \
        const FieldField<fvsPatchField, Type>&, \
        const FieldField<fvsPatchField, Type>& \
    ); \
    template void max<fvsPatchField, Type> \
    ( \
        FieldField<fvsPatchField, Type>&, \
        const FieldField<fvsPatchField, Type>&, \
        const FieldField<fvsPatchField, Type>& \
    ); \
    \
    /* fvPatchField min/max (single field) */ \
    template void min \
    ( \
        fvPatchField<Type>&, \
        const fvPatchField<Type>&, \
        const fvPatchField<Type>& \
    ); \
    template void max \
    ( \
        fvPatchField<Type>&, \
        const fvPatchField<Type>&, \
        const fvPatchField<Type>& \
    ); \
    \
    /* fvsPatchField min/max (single field) */ \
    template void min \
    ( \
        fvsPatchField<Type>&, \
        const fvsPatchField<Type>&, \
        const fvsPatchField<Type>& \
    ); \
    template void max \
    ( \
        fvsPatchField<Type>&, \
        const fvsPatchField<Type>&, \
        const fvsPatchField<Type>& \
    ); \
    \
    /* Field min/max (direct Field usage) */ \
    template void min<Field, Type> \
    ( \
        Field<Type>&, \
        const Field<Type>&, \
        const Field<Type>& \
    ); \
    template void max<Field, Type> \
    ( \
        Field<Type>&, \
        const Field<Type>&, \
        const Field<Type>& \
    );

// Instantiate for all required types
INSTANTIATE_MIN_MAX(scalar);
INSTANTIATE_MIN_MAX(vector);
INSTANTIATE_MIN_MAX(tensor);
INSTANTIATE_MIN_MAX(symmTensor);
INSTANTIATE_MIN_MAX(sphericalTensor);

#undef INSTANTIATE_MIN_MAX

} // namespace Foam
