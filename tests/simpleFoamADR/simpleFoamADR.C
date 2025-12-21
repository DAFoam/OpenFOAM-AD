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
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    argList::addOption(
        "dvName",
        "U0",
        "name of the design variable");

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
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
    label pointI = 196;
    label comp = 1;
    
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
        Info << "Seed mesh coords " << meshPoints[pointI] << endl;
        tape.registerInput(meshPoints[pointI][comp]);
        mesh.movePoints(meshPoints);
    }
    
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.controlledUpdate();

        if (mesh.changing())
        {
            MRF.update();
        }

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);

        pWall = 0.0;
        forAll(p.boundaryField()[patchIWalls], faceI)
        {
            pWall += p.boundaryField()[patchIWalls][faceI];
        }
        Info << "pWall: " << pWall << endl;
    }

    if (dvName != "None")
    {
        tape.registerOutput(pWall);
        tape.setPassive();
    
        pWall.setGradient(1.0);
        tape.evaluate();
    
        if (dvName == "U0")
        {
            scalar total = U0.getGradient();
            scalar ref = 1068.670037013205;
            Info << "dpWall/dU0 ADR: " << total << endl;
            Info << "dpWall/dU0 REF: " << ref << endl;
            if (mag(total - ref) / ref < 1e-10)
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
            scalar total = meshPoints[pointI][comp].getGradient();
            scalar ref = -5922.669045154439;
            Info << "dpWall/dXv ADR: " << total << endl;
            Info << "dpWall/dXv REF: " << ref << endl;
            if (mag(total - ref) / ref < 1e-10)
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
