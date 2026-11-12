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
    Foam::fv::krakenScalarFieldSource

Description
    Power source for coupling with Cerberus. Sets a power distribution
    provided by Cerberus. Cerberus provides the power [W] for each cell in
    the mesh, indexed as in the undecomposed mesh. Cell temperatures and
    densities are provided as output. Supports only single phase
    simulations.

    Input and output field names are appended with the name of the model,
    e.g. core in the example below.

    Input fields:
        foam_if_power_<name>        Power [W] for each cell

    Output fields:
        foam_of_density_<name>      Density [kg/m^3] for each cell
        foam_of_temperature_<name>  Temperature [K] for each cell

Usage
    Example usage:
    \verbatim
    core
    {
        type            krakenScalarFieldSource;
        libs            ("libkrakenIO.so");
    }
    \endverbatim

SourceFiles
    krakenScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef krakenScalarFieldSource_H
#define krakenScalarFieldSource_H

#include "fvModel.H"
#include "krakenModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class krakenScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class krakenScalarFieldSource
:
    public fvModel,
    public krakenModelBase
{
    // Private Data

        //- Volumetric heat source field [W/m^3]
        mutable volScalarField pVol_;

        //- Number of cells before decomposition
        label nGlobalCells_;

        //- Local to global mesh indexing (needed on all procs)
        labelList localToGlobalIdxs_;

        //- Local to global mesh indexing for all procs (needed on master)
        labelListList localToGlobalIdxsAllProcs_;


    // Private Member Functions

        //- Gather the given field to the output array on master using the
        //  global (undecomposed) cell indexing
        void updateOutputArray
        (
            List<scalar>& outputArray,
            const volScalarField& field
        ) const;

public:

    //- Runtime type information
    TypeName("krakenScalarFieldSource");


    // Constructors

        //- Construct from components
        krakenScalarFieldSource
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        krakenScalarFieldSource(const krakenScalarFieldSource&) = delete;


    //- Destructor
    virtual ~krakenScalarFieldSource()
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
        void operator=(const krakenScalarFieldSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
