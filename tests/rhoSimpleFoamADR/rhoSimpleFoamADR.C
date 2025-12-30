/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    rhoSimpleFoam

Group
    grpCompressibleSolvers

Description
    Steady-state solver for compressible turbulent flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for compressible turbulent flow."
    );

    argList::addOption(
        "dvName",
        "U0",
        "name of the design variable");

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    word dvName = "None";
    if (args.found("dvName"))
    {
        dvName = word(args.lookup("dvName")());
    }
    else
    {
        Info << "dvName not set!" << endl;
    }

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    scalar U0 = 10.0;
    scalar pWall = 0.0;
    label patchIWalls = mesh.boundaryMesh().findPatchID("walls");
    pointField meshPoints = mesh.points();
    
    if (dvName == "U0")
    {
        tape.setActive();
        label patchIInlet = mesh.boundaryMesh().findPatchID("inlet");
        
        tape.registerInput(U0);
        forAll(U.boundaryField()[patchIInlet], faceI)
        {
            U.boundaryFieldRef()[patchIInlet][faceI][0] = U0;
        }
        U.correctBoundaryConditions();
    }
    else if (dvName == "Xv")
    {
        tape.setActive();
        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                label pointI = 69;
                label comp = 1;
                Info << "Seed mesh coords " << meshPoints[pointI] << endl;
                tape.registerInput(meshPoints[pointI][comp]);
            }
        }
        else
        {
            label pointI = 195;
            label comp = 1;
            Info << "Seed mesh coords " << meshPoints[pointI] << endl;
            tape.registerInput(meshPoints[pointI][comp]);
        }
        mesh.movePoints(meshPoints);
    }

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "EEqn.H"

        if (simple.consistent())
        {
            #include "pcEqn.H"
        }
        else
        {
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);

        pWall = 0.0;
        forAll(p.boundaryField()[patchIWalls], faceI)
        {
            pWall += p.boundaryField()[patchIWalls][faceI];
        }
        reduce(pWall, sumOp<scalar>());
        Info << "pWall: " << pWall << endl;
    }

    if (dvName != "None")
    {
        tape.registerOutput(pWall);
        tape.setPassive();
    
        if (Pstream::master())
        {
            pWall.setGradient(1.0);
        }
        tape.evaluate();
    
        if (dvName == "U0")
        {
            scalar total = U0.getGradient();
            reduce(total, sumOp<scalar>());
            scalar ref = 1265.293612277306;
            Info << "dpWall/dU0 ADR: " << total << endl;
            Info << "dpWall/dU0 REF: " << ref << endl;
            if (mag(total - ref) / ref < 1e-8)
            {
                Info << "dpWall/dU0 test passed!" << endl;
            }
            else
            {
                Info << "dpWall/dU0 test failed!" << endl;
                return 1;
            }
        }
        else if (dvName == "Xv")
        {
            scalar total = 0.0; 
            if (Pstream::parRun())
            {
                if (Pstream::master())
                {
                    label pointI = 69;
                    label comp = 1;
                    total = meshPoints[pointI][comp].getGradient();
                }
            }
            else
            {
                label pointI = 195;
                label comp = 1;
                total = meshPoints[pointI][comp].getGradient();
            }
            // in parallel, only master has the total value, other proces have total=0
            // so we need to reduce() the total value to all procs
            reduce(total, sumOp<scalar>());

            scalar ref = -5120.116613136466;
            Info << "dpWall/dXv ADR: " << total << endl;
            Info << "dpWall/dXv REF: " << ref << endl;
            if (mag(total - ref) / mag(ref) < 1e-7)
            {
                Info << "dpWall/dXv test passed!" << endl;
            }
            else
            {
                Info << "dpWall/dXv test failed!" << endl;
                return 1;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
