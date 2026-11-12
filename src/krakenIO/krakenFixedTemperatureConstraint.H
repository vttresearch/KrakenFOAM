/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

Class
    Foam::fv::krakenFixedTemperatureConstraint

Description
    Fixed temperature constraint. Sets a temperature distribution provided by
    Cerberus (associated with a set of points) to a specified set of cells.
    The point set is read from the dictionary and each cell is associated
    with the nearest point. Cerberus provides the temperature distribution.
    The routine integrates the heat flow rate and returns it back to
    Cerberus. Intended to be used with a fixed solid phase in multiphase
    simulations. In addition, volume averaged temperature and density fields
    (averaged over outputPhases) are also provided as output.

    Input and output field names are appended with the name of the
    constraint, e.g. core in the example below.

    Input fields:
        foam_if_temperature_<name>  Temperature [K] for each point

    Output fields:
        foam_of_density_<name>      Density [kg/m^3] for each point, averaged
                                    over outputPhases
        foam_of_temperature_<name>  Temperature [K] for each point, averaged
                                    over outputPhases
        foam_of_power_<name>        Heat transfer rate [W] for each point

Usage
    \table
        Property      | Description           | Required | Default value
        cellZone      | cell zone selection   | yes      |
        inputPhase    | constrained phase     | yes      |
        outputPhases  | list of output phases | yes      |
        points        | list of points        | yes      |
    \endtable

    Currently supports cell selection from:
    \verbatim
        // Select all cells
        cellZone    all;

        // Select the cells within the given cellZone
        cellZone    core;
    \endverbatim

Example
    \verbatim
    core
    {
        type            krakenFixedTemperatureConstraint;
        libs            ("libkrakenIO.so");

        cellZone        fluid;

        inputPhase      solid;
        outputPhases    (gas liquid);

        points
        (
            (0.2150364 -0.4300728 0.05)
            (0.2150364 -0.4300728 0.15)
            .
            .
            .
        );
    }
    \endverbatim

SourceFiles
    krakenFixedTemperatureConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef krakenFixedTemperatureConstraint_H
#define krakenFixedTemperatureConstraint_H

#include "fvConstraint.H"
#include "krakenModelBase.H"
#include "fvCellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
             Class krakenFixedTemperatureConstraint Declaration
\*---------------------------------------------------------------------------*/

class krakenFixedTemperatureConstraint
:
    public fvConstraint,
    public krakenModelBase
{
    // Private Data

        //- The cellZone the fvConstraint applies to
        fvCellZone zone_;

        //- The constraint is applied for this phase
        word inputPhase_;

        //- Fluid temperature and density are averaged over these phases
        wordList outputPhases_;

        //- Number of points in the transferred data
        label nSamplePoints_;

        //- Index of the nearest sample point for each cell
        labelList samplePointIndices_;

        //- Volume of the zone cells associated with each sample point
        scalarList samplePointVolumes_;

        //- Volume of the output phases associated with each sample point
        scalarList samplePointSumAlphas_;


public:

    //- Runtime type information
    TypeName("krakenFixedTemperatureConstraint");


    // Constructors

        //- Construct from components
        krakenFixedTemperatureConstraint
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        krakenFixedTemperatureConstraint
        (
            const krakenFixedTemperatureConstraint&
        ) = delete;


    //- Destructor
    virtual ~krakenFixedTemperatureConstraint()
    {}


    // Member Functions

        // Constrain

            //- Return the list of fields constrained by the fvConstraint
            virtual wordList constrainedFields() const;

            //- Constrain energy equation to fix the temperature
            virtual bool constrain
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // Kraken model base class functions

            //- Define krakenIO field templates
            virtual void createFieldTemplates
            (
                List<krakenIO::fieldTemplate>& inputFields,
                List<krakenIO::fieldTemplate>& outputFields
            ) const;

            //- Update/gather values of krakenIO arrays to master
            virtual void updateOutputFields();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const krakenFixedTemperatureConstraint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
