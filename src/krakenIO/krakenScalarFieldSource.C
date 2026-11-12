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

\*---------------------------------------------------------------------------*/

#include "krakenScalarFieldSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(krakenScalarFieldSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        krakenScalarFieldSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::krakenScalarFieldSource::updateOutputArray
(
    List<scalar>& outputArray,
    const volScalarField& field
) const
{
    outputArray = 0;

    if (Pstream::parRun())
    {
        List<List<scalar>> gatheredValues(Pstream::nProcs());
        gatheredValues[Pstream::myProcNo()] = field;
        Pstream::gatherList(gatheredValues);

        if (Pstream::master())
        {
            forAll(gatheredValues, j)
            {
                const labelList& globalIdxs = localToGlobalIdxsAllProcs_[j];

                forAll(globalIdxs, i)
                {
                    outputArray[globalIdxs[i]] = gatheredValues[j][i];
                }
            }
        }
    }
    else
    {
        outputArray = field;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::krakenScalarFieldSource::krakenScalarFieldSource
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    krakenModelBase(),
    pVol_
    (
        IOobject
        (
            "pVol." + sourceName,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),
    nGlobalCells_(0),
    localToGlobalIdxs_(),
    localToGlobalIdxsAllProcs_(Pstream::nProcs())
{
    // Calculate the number of cells in the mesh before decomposition
    nGlobalCells_ = returnReduce(mesh.nCells(), sumOp<label>());

    // Read local to global mesh indexing if this a parallel run and gather
    // all index mappings to the master
    if (Pstream::parRun())
    {
        localToGlobalIdxs_ =
            labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        localToGlobalIdxsAllProcs_[Pstream::myProcNo()] = localToGlobalIdxs_;
        Pstream::gatherList(localToGlobalIdxsAllProcs_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::krakenScalarFieldSource::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::krakenScalarFieldSource::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    if (eqn.psi().name() == thermo.he().name())
    {
        eqn += pVol_;
    }
}


void Foam::fv::krakenScalarFieldSource::correct()
{
    // Reference to the IO object (constness removed)
    krakenIO& IO =
        const_cast<krakenIO&>
        (
            mesh().time().lookupObject<krakenIO>("krakenIO")
        );

    // Reference to the krakenIO input array
    const List<scalar>& power =
        IO.getFieldData("foam_if_power_" + this->name());

    // Convert powers to power densities and update the source term field
    const bool isParallel = Pstream::parRun();
    const scalarField& V = mesh().V();

    forAll(V, i)
    {
        const label globali = isParallel ? localToGlobalIdxs_[i] : i;
        pVol_[i] = power[globali]/V[i];
    }

    Info<< "power"
        << ": min = " << gMin(power)
        << ", mean = " << gAverage(power)
        << ", max = " << gMax(power)
        << endl;

    Info<< "pVol_"
        << ": min = " << gMin(pVol_)
        << ", mean = " << gAverage(pVol_)
        << ", max = " << gMax(pVol_)
        << endl;

    // Reset flag to update output fields
    setOutputUpdateRequired();
}


bool Foam::fv::krakenScalarFieldSource::movePoints()
{
    return true;
}


void Foam::fv::krakenScalarFieldSource::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::krakenScalarFieldSource::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::krakenScalarFieldSource::distribute
(
    const polyDistributionMap& map
)
{}


void Foam::fv::krakenScalarFieldSource::createFieldTemplates
(
    List<krakenIO::fieldTemplate>& inputFields,
    List<krakenIO::fieldTemplate>& outputFields
) const
{
    const word& name = this->name();
    const word& regionName = this->mesh().name();
    const label nCells = nGlobalCells_;

    krakenIO::fieldTemplate powerField;
    powerField.fieldName = "foam_if_power_" + name;
    powerField.meshName = "foam_mesh_" + name;
    powerField.modelName = name;
    powerField.regionName = regionName;
    powerField.nCells = nCells;
    powerField.dataType = 1;
    powerField.unitVector[0] = 1;
    powerField.unitVector[1] = 2;
    powerField.unitVector[2] = -3;
    powerField.unitVector[3] = 0;
    powerField.unitVector[4] = 0;
    powerField.multiplier = 1;
    powerField.offset = 0;
    powerField.fieldData.setSize(nCells, 0);
    inputFields.append(powerField);

    krakenIO::fieldTemplate densityField;
    densityField.fieldName = "foam_of_density_" + name;
    densityField.meshName = "foam_mesh_" + name;
    densityField.modelName = name;
    densityField.regionName = regionName;
    densityField.nCells = nCells;
    densityField.dataType = 1;
    densityField.unitVector[0] = 1;
    densityField.unitVector[1] = -3;
    densityField.unitVector[2] = 0;
    densityField.unitVector[3] = 0;
    densityField.unitVector[4] = 0;
    densityField.multiplier = 1;
    densityField.offset = 0;
    densityField.fieldData.setSize(nCells, 0);
    outputFields.append(densityField);

    krakenIO::fieldTemplate temperatureField;
    temperatureField.fieldName = "foam_of_temperature_" + name;
    temperatureField.meshName = "foam_mesh_" + name;
    temperatureField.modelName = name;
    temperatureField.regionName = regionName;
    temperatureField.nCells = nCells;
    temperatureField.dataType = 1;
    temperatureField.unitVector[0] = 0;
    temperatureField.unitVector[1] = 0;
    temperatureField.unitVector[2] = 0;
    temperatureField.unitVector[3] = 1;
    temperatureField.unitVector[4] = 0;
    temperatureField.multiplier = 1;
    temperatureField.offset = 0;
    temperatureField.fieldData.setSize(nCells, 0);
    outputFields.append(temperatureField);
}


void Foam::fv::krakenScalarFieldSource::updateOutputFields()
{
    if (!outputUpdateRequired())
    {
        return;
    }

    clearOutputUpdateRequired();

    // Reference to the IO object (constness removed)
    krakenIO& IO =
        const_cast<krakenIO&>
        (
            mesh().time().lookupObject<krakenIO>("krakenIO")
        );

    // References to the krakenIO output arrays
    List<scalar>& temperature =
        IO.setFieldData("foam_of_temperature_" + this->name());
    List<scalar>& density =
        IO.setFieldData("foam_of_density_" + this->name());

    // Temperatures and densities for output
    const volScalarField& T = mesh().lookupObject<volScalarField>("T");
    const volScalarField& rho = mesh().lookupObject<volScalarField>("rho");

    updateOutputArray(temperature, T);
    updateOutputArray(density, rho);
}


// ************************************************************************* //
