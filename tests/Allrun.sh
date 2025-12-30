#!/bin/sh
set -e  # exit if any command fails

OVERSUBSCRIBE=""

mpi_version=$(mpirun -V --help 2>&1 | head -1)

if echo "$mpi_version" | grep -qi "open.*mpi"; then
    OVERSUBSCRIBE="-oversubscribe"
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

cd ../ChannelAMI
rm -rf processor*
decomposePar
simpleFoam${WM_AD_MODE} -dvName XvAMI
mpirun $OVERSUBSCRIBE -np 4 simpleFoam${WM_AD_MODE} -dvName XvAMI -parallel

# test rhoSimpleFoam
cd ../Channel
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

