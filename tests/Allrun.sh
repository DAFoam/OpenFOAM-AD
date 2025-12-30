#!/bin/sh
set -e  # exit if any command fails

# Test if --oversubscribe is supported
if mpirun --help 2>&1 | grep -q "oversubscribe"; then
    OVERSUBSCRIBE="--oversubscribe"
else
    OVERSUBSCRIBE=""
fi

cd Channel

# test simpleFoam
cp 0.incompressible/* 0/
cp system.incompressible/* system/
rm -rf processor*
decomposePar
simpleFoam${WM_AD_MODE} -dvName U0
simpleFoam${WM_AD_MODE} -dvName Xv
mpirun $OVERSUBSCRIBE -np 4 simpleFoam${WM_AD_MODE} -dvName U0 -parallel
mpirun $OVERSUBSCRIBE -np 4 simpleFoam${WM_AD_MODE} -dvName Xv -parallel

# test rhoSimpleFoam
cp 0.compressible/* 0/
cp system.compressible/* system/
rm -rf processor*
decomposePar
rhoSimpleFoam${WM_AD_MODE} -dvName U0
rhoSimpleFoam${WM_AD_MODE} -dvName Xv
mpirun $OVERSUBSCRIBE -np 4 rhoSimpleFoam${WM_AD_MODE} -dvName U0 -parallel
mpirun $OVERSUBSCRIBE -np 4 rhoSimpleFoam${WM_AD_MODE} -dvName Xv -parallel

# unit tests for AMPI
if [ "$WM_AD_MODE" = "ADR" ]; then
    mpirun $OVERSUBSCRIBE -np 4 MPIUnitTests -parallel
fi

