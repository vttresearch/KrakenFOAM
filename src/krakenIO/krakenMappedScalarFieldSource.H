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
    Foam::fv::krakenMappedScalarFieldSource

Description
    Mapped power source for coupling with Cerberus. Sets a power distribution
    provided by Cerberus (associated with a set of points) to a specified set
    of cells. The point set is read from the dictionary and each cell is
    associated with the nearest point. Cerberus provides the power
    distribution. Volume averaged temperature and density fields (averaged
    over outputPhases) are provided as output. Supports only multiphase
    simulations.

    Input and output field names are appended with the name of the model,
    e.g. core in the example below.

    Input fields:
        foam_if_power_<name>        Power or heat transfer rate [W] for each
                                    point

    Output fields:
        foam_of_density_<name>      Density [kg/m^3] for each point, averaged
                                    over outputPhases
        foam_of_temperature_<name>  Temperature [K] for each point, averaged
                                    over outputPhases

Usage
    \table
        Property      | Description           | Required | Default value
        cellZone      | cell zone selection   | yes      |
        inputPhase    | phase of the source   | yes      |
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
        type            krakenMappedScalarFieldSource;
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
    krakenMappedScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef krakenMappedScalarFieldSource_H
#define krakenMappedScalarFieldSource_H

#include "fvModel.H"
#include "krakenModelBase.H"
#include "fvCellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class krakenMappedScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class krakenMappedScalarFieldSource
:
    public fvModel,
    public krakenModelBase
{
    // Private Data

        //- The cellZone the fvModel applies to
        fvCellZone zone_;

        //- The source is put into this phase
        word inputPhase_;

        //- Fluid temperature and density are averaged over these phases
        wordList outputPhases_;

        //- Number of points in the transferred data
        label nSamplePoints_;

        //- Volumetric heat source field [W/m^3]
        mutable volScalarField pVol_;

        //- Index of the nearest sample point for each cell
        labelList samplePointIndices_;

        //- Volume of the zone cells associated with each sample point
        scalarList samplePointVolumes_;

        //- Volume of the output phases associated with each sample point
        scalarList samplePointSumAlphas_;


public:

    //- Runtime type information
    TypeName("krakenMappedScalarFieldSource");


    // Constructors

        //- Construct from components
        krakenMappedScalarFieldSource
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        krakenMappedScalarFieldSource
        (
            const krakenMappedScalarFieldSource&
        ) = delete;


    //- Destructor
    virtual ~krakenMappedScalarFieldSource()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to the transport equation (not implemented)
            virtual void addSup(fvMatrix<scalar>& eqn) const
            {
                NotImplemented;
            }

            //- Source term to the compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to the phase compressible energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const
            {
                return this->addSup(rho, he, eqn);
            }


        // Correction

            //- Correct the fvModel, i.e. update the source from the
            //  Cerberus input
            virtual void correct();


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
        void operator=(const krakenMappedScalarFieldSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
