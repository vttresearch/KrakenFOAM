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

#include "krakenFixedTemperatureConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "heatTransferSystem.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(krakenFixedTemperatureConstraint, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        krakenFixedTemperatureConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::krakenFixedTemperatureConstraint::krakenFixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    krakenModelBase(),
    zone_(mesh, dict),
    inputPhase_(dict.lookup<word>("inputPhase")),
    outputPhases_(dict.lookup<wordList>("outputPhases"))
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

Foam::wordList
Foam::fv::krakenFixedTemperatureConstraint::constrainedFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, inputPhase_)
        );

    return wordList(1, thermo.he().name());
}


bool Foam::fv::krakenFixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, inputPhase_)
        );

    // Constant reference to the selected set of cells
    const labelUList& sCells = zone_.zone();

    // Reference to the IO object (constness removed)
    krakenIO& IO =
        const_cast<krakenIO&>
        (
            mesh().time().lookupObject<krakenIO>("krakenIO")
        );

    // Get the data related to the point field
    const List<scalar>& T =
        IO.getFieldData("foam_if_temperature_" + this->name());

    // Cell temperature for the constrained cells
    scalarField Tc(sCells.size(), scalar(0));

    // Map temperature
    forAll(sCells, i)
    {
        Tc[i] = T[samplePointIndices_[sCells[i]]];
    }

    // Temperature converted to enthalpy
    const scalarField hec(thermo.he(Tc, sCells));

    // Constrain
    eqn.setValues(sCells, hec);

    DebugInfo<< "Setting Tc " << gMin(Tc) << " - " << gMax(Tc) << endl;

    // Reset flag to update output fields
    setOutputUpdateRequired();

    return sCells.size();
}


bool Foam::fv::krakenFixedTemperatureConstraint::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::krakenFixedTemperatureConstraint::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::krakenFixedTemperatureConstraint::mapMesh
(
    const polyMeshMap& map
)
{
    zone_.mapMesh(map);
}


void Foam::fv::krakenFixedTemperatureConstraint::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


void Foam::fv::krakenFixedTemperatureConstraint::createFieldTemplates
(
    List<krakenIO::fieldTemplate>& inputFields,
    List<krakenIO::fieldTemplate>& outputFields
) const
{
    const word& name = this->name();
    const word& regionName = this->mesh().name();
    const label nPoints = nSamplePoints_;

    krakenIO::fieldTemplate TField;
    TField.fieldName = "foam_if_temperature_" + name;
    TField.meshName = "foam_mesh_" + name;
    TField.modelName = name;
    TField.regionName = regionName;
    TField.nCells = nPoints;
    TField.dataType = 1;
    TField.unitVector[0] = 0;
    TField.unitVector[1] = 0;
    TField.unitVector[2] = 0;
    TField.unitVector[3] = 1;
    TField.unitVector[4] = 0;
    TField.multiplier = 1;
    TField.offset = 0;
    TField.fieldData.setSize(nPoints, 0);
    inputFields.append(TField);

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

    krakenIO::fieldTemplate powerField;
    powerField.fieldName = "foam_of_power_" + name;
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
    outputFields.append(powerField);
}


void Foam::fv::krakenFixedTemperatureConstraint::updateOutputFields()
{
    if (!outputUpdateRequired())
    {
        return;
    }

    clearOutputUpdateRequired();

    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, inputPhase_)
        );

    const heatTransferSystem& hts =
        mesh().lookupObject<heatTransferSystem>("heatTransfer");

    const phaseModel& phase =
        mesh().lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", inputPhase_)
        );

    // Heat transfer table (as in EEqns.H of the solver)
    autoPtr<HashPtrTable<fvScalarMatrix>> heatTransferPtr =
        hts.heatTransfer();
    HashPtrTable<fvScalarMatrix>& heatTransfer = heatTransferPtr();

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
    List<scalar>& P =
        IO.setFieldData("foam_of_power_" + this->name());

    // Update the output fields
    temperature = 0;
    density = 0;
    P = 0;

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

    // Power density [W/m^3]
    const volScalarField Qv(*heatTransfer[phase.name()] & thermo.he());

    Info<< " integral(Qv) = " << fvc::domainIntegrate(Qv).value()
        << ", min(Qv) = " << gMin(Qv)
        << ", max(Qv) = " << gMax(Qv) << endl;

    // Power [W]
    forAll(sCells, i)
    {
        P[samplePointIndices_[sCells[i]]] -= Qv[sCells[i]]*V[sCells[i]];
    }

    reduce(P, sumOp<List<scalar>>());

    Info<< " sum(P) = " << sum(P)
        << ", min(P) = " << min(P)
        << ", max(P) = " << max(P) << endl;
}


// ************************************************************************* //
