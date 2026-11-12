/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Modified by VTT
Copyright (C) 2026 VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "setDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::setDeltaT(Time& runTime, const PtrList<solver>& solvers)
{
    if
    (
        runTime.timeIndex() == 0
     && runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
    )
    {
        bool transient = false;
        scalar deltaT = vGreat;

        forAll(solvers, i)
        {
            if (solvers[i].transient())
            {
                transient = true;
                deltaT = min(deltaT, solvers[i].maxDeltaT());
            }
        }

        if (transient)
        {
            deltaT = min(deltaT, runTime.functionObjects().maxDeltaT());
        }

        if (transient && deltaT < rootVGreat)
        {
            runTime.setDeltaT(min(runTime.deltaTValue(), deltaT));
        }
    }
}


void Foam::adjustDeltaT(Time& runTime, const PtrList<solver>& solvers)
{
    // Update the time-step limited by the solvers maxDeltaT
    if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        bool transient = false;
        scalar deltaT = vGreat;

        forAll(solvers, i)
        {
            if (solvers[i].transient())
            {
                transient = true;
                deltaT = min(deltaT, solvers[i].maxDeltaT());
            }
        }

        if (transient)
        {
            deltaT = min(deltaT, runTime.functionObjects().maxDeltaT());
        }

        if (transient && deltaT < rootVGreat)
        {
            // Adjust original dt set by Cerberus to remain below max dt, and
            // also match it with Cerberus time interval. Note deltaT = max dt
            if (runTime.deltaTValue() > deltaT)
            {
                const scalar remaining =
                    runTime.endTime().value() - runTime.value();
                const long nStep = long(remaining/deltaT) + 1;
                runTime.setDeltaTNoAdjust(remaining/nStep);
            }

            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }
    }
}


// ************************************************************************* //
