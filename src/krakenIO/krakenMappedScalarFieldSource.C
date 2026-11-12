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

#include "krakenMappedScalarFieldSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(krakenMappedScalarFieldSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        krakenMappedScalarFieldSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::krakenMappedScalarFieldSource::krakenMappedScalarFieldSource
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    krakenModelBase(),
    zone_(mesh, dict),
    inputPhase_(dict.lookup<word>("inputPhase")),
    outputPhases_(dict.lookup<wordList>("outputPhases")),
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
    )
{
    // Interpolation (to nearest)

    // Get the list of sample points
    const List<vector> samplePoints(dict.lookup<List<vector>>("points"));
    nSamplePoints_ = samplePoints.size();

    Info<< "    Interpolating " << mesh.globalData().nTotalCells()
        << " cells against " << samplePoints.size() << " points... ";

    // Save samplePoint indices for each cell
    const volVectorField& cellCentres = mesh.C();
    samplePointIndices_.setSize(cellCentres.internalField().size());

    forAll(cellCentres, i)
    {
        label nearesti = -1;
        scalar dist = great;
        const vector& C = cellCentres[i];

        forAll(samplePoints, j)
        {
            const scalar magV = mag(C - samplePoints[j]);

            if (magV <= dist)
            {
                dist = magV;
                nearesti = j;
            }
        }

        if (nearesti < 0)
        {
            FatalErrorInFunction
                << "Nearest samplePoint not found!"
                << abort(FatalError);
        }
        else
        {
            samplePointIndices_[i] = nearesti;
        }
    }

    Info<< "Complete!" << endl;


    // Sum cell volumes addressed to each sample point

    // Constant reference to cell volumes
    const scalarField& V = mesh.V();

    // Constant reference to the selected set of cells
    const labelUList& sCells = zone_.zone();

    samplePointVolumes_.setSize(nSamplePoints_);
    samplePointVolumes_ = 0;

    Info<< "    Summing sample point volumes over "
        << returnReduce(sCells.size(), sumOp<label>()) << " cells... ";

    forAll(sCells, i)
    {
        samplePointVolumes_[samplePointIndices_[sCells[i]]] += V[sCells[i]];
    }

    reduce(samplePointVolumes_, sumOp<List<scalar>>());

    Info<< "Complete!" << endl;

    Info<< "    Cell volumes"
        << ": min = " << gMin(samplePointVolumes_)
        << ", mean = " << gAverage(samplePointVolumes_)
        << ", max = " << gMax(samplePointVolumes_)
        << endl;


    // Sum fluid volumes addressed to each sample point

    if (outputPhases_.size() < 1)
    {
        FatalErrorInFunction
            << "Invalid outputPhases list!"
            << abort(FatalError);
    }

    samplePointSumAlphas_.setSize(nSamplePoints_);
    samplePointSumAlphas_ = 0;

    Info<< "    Summing sample point output phase volumes over "
        << returnReduce(sCells.size(), sumOp<label>()) << " cells... ";

    forAll(outputPhases_, phasei)
    {
        const volScalarField& alphai = mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", outputPhases_[phasei])
        );

        forAll(sCells, i)
        {
            samplePointSumAlphas_[samplePointIndices_[sCells[i]]] +=
                alphai[sCells[i]]*V[sCells[i]];
        }
    }

    reduce(samplePointSumAlphas_, sumOp<List<scalar>>());

    Info<< "Complete!" << endl;

    Info<< "    Output phase volumes"
        << ": min = " << gMin(samplePointSumAlphas_)
        << ", mean = " << gAverage(samplePointSumAlphas_)
        << ", max = " << gMax(samplePointSumAlphas_)
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::krakenMappedScalarFieldSource::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, inputPhase_)
        );

    return wordList(1, thermo.he().name());
}


void Foam::fv::krakenMappedScalarFieldSource::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, inputPhase_)
        );

    if (eqn.psi().name() == thermo.he().name())
    {
        eqn += pVol_;
    }
}


void Foam::fv::krakenMappedScalarFieldSource::correct()
{
    // Constant reference to the selected set of cells
    const labelUList& sCells = zone_.zone();

    // Reference to the IO object (constness removed)
    krakenIO& IO =
        const_cast<krakenIO&>
        (
            mesh().time().lookupObject<krakenIO>("krakenIO")
        );

    // Reference to the input point field
    const List<scalar>& power =
        IO.getFieldData("foam_if_power_" + this->name());

    // Update the input field, i.e. the source term
    forAll(sCells, i)
    {
        const label celli = sCells[i];
        const label pointi = samplePointIndices_[celli];
        const scalar V = samplePointVolumes_[pointi];

        if (V >= 1e-16)
        {
            // Convert the power [W] to power density [W/m^3]
            pVol_[celli] = power[pointi]/V;
        }
        else
        {
            pVol_[celli] = 0;
        }
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


bool Foam::fv::krakenMappedScalarFieldSource::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::krakenMappedScalarFieldSource::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::krakenMappedScalarFieldSource::mapMesh
(
    const polyMeshMap& map
)
{
    zone_.mapMesh(map);
}


void Foam::fv::krakenMappedScalarFieldSource::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


void Foam::fv::krakenMappedScalarFieldSource::createFieldTemplates
(
    List<krakenIO::fieldTemplate>& inputFields,
    List<krakenIO::fieldTemplate>& outputFields
) const
{
    const word& name = this->name();
    const word& regionName = this->mesh().name();
    const label nPoints = nSamplePoints_;

    krakenIO::fieldTemplate powerField;
    powerField.fieldName = "foam_if_power_" + name;
    powerField.meshName = "foam_mesh_" + name;
    powerField.modelName = name;
    powerField.regionName = regionName;
    powerField.nCells = nPoints;
    powerField.dataType = 1;
    powerField.unitVector[0] = 1;
    powerField.unitVector[1] = 2;
    powerField.unitVector[2] = -3;
    powerField.unitVector[3] = 0;
    powerField.unitVector[4] = 0;
    powerField.multiplier = 1;
    powerField.offset = 0;
    powerField.fieldData.setSize(nPoints, 0);
    inputFields.append(powerField);

    krakenIO::fieldTemplate densityField;
    densityField.fieldName = "foam_of_density_" + name;
    densityField.meshName = "foam_mesh_" + name;
    densityField.modelName = name;
    densityField.regionName = regionName;
    densityField.nCells = nPoints;
    densityField.dataType = 1;
    densityField.unitVector[0] = 1;
    densityField.unitVector[1] = -3;
    densityField.unitVector[2] = 0;
    densityField.unitVector[3] = 0;
    densityField.unitVector[4] = 0;
    densityField.multiplier = 1;
    densityField.offset = 0;
    densityField.fieldData.setSize(nPoints, 0);
    outputFields.append(densityField);

    krakenIO::fieldTemplate temperatureField;
    temperatureField.fieldName = "foam_of_temperature_" + name;
    temperatureField.meshName = "foam_mesh_" + name;
    temperatureField.modelName = name;
    temperatureField.regionName = regionName;
    temperatureField.nCells = nPoints;
    temperatureField.dataType = 1;
    temperatureField.unitVector[0] = 0;
    temperatureField.unitVector[1] = 0;
    temperatureField.unitVector[2] = 0;
    temperatureField.unitVector[3] = 1;
    temperatureField.unitVector[4] = 0;
    temperatureField.multiplier = 1;
    temperatureField.offset = 0;
    temperatureField.fieldData.setSize(nPoints, 0);
    outputFields.append(temperatureField);
}


void Foam::fv::krakenMappedScalarFieldSource::updateOutputFields()
{
    if (!outputUpdateRequired())
    {
        return;
    }

    clearOutputUpdateRequired();

    // Constant reference to the selected set of cells
    const labelUList& sCells = zone_.zone();

    // Constant reference to cell volumes [m^3]
    const scalarField& V = mesh().V();

    // Reference to the IO object (constness removed)
    krakenIO& IO =
        const_cast<krakenIO&>
        (
            mesh().time().lookupObject<krakenIO>("krakenIO")
        );

    // References to the output point fields
    List<scalar>& temperature =
        IO.setFieldData("foam_of_temperature_" + this->name());
    List<scalar>& density =
        IO.setFieldData("foam_of_density_" + this->name());

    // Update the output fields
    temperature = 0;
    density = 0;

    forAll(outputPhases_, phasei)
    {
        const word& phaseName = outputPhases_[phasei];

        const volScalarField& alphai =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("alpha", phaseName)
            );

        const volScalarField& Ti =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("T", phaseName)
            );

        const volScalarField& rhoi =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("rho", phaseName)
            );

        forAll(sCells, i)
        {
            const label celli = sCells[i];
            const label pointi = samplePointIndices_[celli];

            // Volume averaged values
            temperature[pointi] += alphai[celli]*V[celli]*Ti[celli];
            density[pointi] += alphai[celli]*V[celli]*rhoi[celli];
        }
    }

    reduce(temperature, sumOp<List<scalar>>());
    reduce(density, sumOp<List<scalar>>());

    forAll(samplePointSumAlphas_, i)
    {
        const scalar alphaV = samplePointSumAlphas_[i];

        if (alphaV >= 1e-16)
        {
            temperature[i] /= alphaV;
            density[i] /= alphaV;
        }
        else
        {
            temperature[i] = 0;
            density[i] = 0;
        }
    }
}


// ************************************************************************* //
