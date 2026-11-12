/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
#include "krakenIO.H"
#include "Time.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "krakenModelBase.H"

#include "clientIO.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::krakenIO::findField
(
    const List<fieldTemplate>& fields,
    const word& fieldName
)
{
    forAll(fields, i)
    {
        if (fields[i].fieldName == fieldName)
        {
            return i;
        }
    }

    return -1;
}


Foam::string Foam::krakenIO::receiveString(const char* name) const
{
    char charMessage[MAX_STRING_SIZE];
    int64_t stringLength;

    // Global C function, not this member function
    const label status =
        ::receiveString(charMessage, MAX_STRING_SIZE, &stringLength);

    if (status != 0)
    {
        FatalErrorInFunction
            << "receiveString(" << name << ") failed! status: " << status
            << abort(FatalError);
    }

    return string(charMessage, stringLength);
}


void Foam::krakenIO::createFieldTemplates(const PtrList<solver>& solvers)
{
    forAll(solvers, i)
    {
        // List of all fvModels
        const fvModels& models = solvers[i].fvModels();

        forAll(models, j)
        {
            const fvModel& model = models[j];

            if (isA<krakenModelBase>(model))
            {
                const krakenModelBase& krakenModel =
                    dynamic_cast<const krakenModelBase&>(model);

                krakenModel.createFieldTemplates(inputFields_, outputFields_);
            }
        }

        // List of all fvConstraints
        const fvConstraints& constraints = solvers[i].fvConstraints();

        forAll(constraints, j)
        {
            const fvConstraint& constraint = constraints[j];

            if (isA<krakenModelBase>(constraint))
            {
                const krakenModelBase& krakenModel =
                    dynamic_cast<const krakenModelBase&>(constraint);

                krakenModel.createFieldTemplates(inputFields_, outputFields_);
            }
        }
    }
}


void Foam::krakenIO::terminate()
{
    if (Pstream::master())
    {
        Info<< "Disconnecting from Cerberus..." << endl;

        const label status = closeConnection();
        if (status != 0)
        {
            FatalErrorInFunction
                << "closeConnection() failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::returnCurrentTime()
{
    // This is used only in depletion calculations and since depletion is not
    // simulated in OpenFOAM give an error
    FatalErrorInFunction
        << "OpenFOAM does not support depletion specific signals"
        << abort(FatalError);
}


void Foam::krakenIO::returnCurrentTimeInterval()
{
    const Time& runTime = this->db().time();
    const double startTime = runTime.value();
    const double endTime = runTime.endTime().value();

    if (Pstream::master())
    {
        Info<< "Return current time interval to Cerberus: "
            << startTime << " " << endTime << endl;

        label status = sendDouble(startTime);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendDouble(startTime) failed! status: " << status
                << abort(FatalError);
        }

        status = sendDouble(endTime);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendDouble(endTime) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::setCurrentTime()
{
    // This is used only in depletion calculations and since depletion is not
    // simulated in OpenFOAM give an error
    FatalErrorInFunction
        << "OpenFOAM does not support depletion specific signals"
        << abort(FatalError);
}


void Foam::krakenIO::setCurrentTimeInterval()
{
    advanceToTimeInterval();
}


void Foam::krakenIO::suggestNextTime()
{
    // This is used only in depletion calculations and since depletion is not
    // simulated in OpenFOAM give an error
    FatalErrorInFunction
        << "OpenFOAM does not support depletion specific signals"
        << abort(FatalError);
}


void Foam::krakenIO::suggestNextTimeInterval()
{
    FatalErrorInFunction
        << "OpenFOAM does not support suggesting a time step size"
        << abort(FatalError);
}


void Foam::krakenIO::advanceToTime()
{
    // This is used only in depletion calculations and since depletion is not
    // simulated in OpenFOAM give an error
    FatalErrorInFunction
        << "OpenFOAM does not support depletion specific signals"
        << abort(FatalError);
}


void Foam::krakenIO::advanceToTimeInterval()
{
    const Time& runTime = this->db().time();
    double startTime, endTime;

    if (Pstream::master())
    {
        Info<< "Receive interval start and end times from Cerberus" << endl;

        label status = receiveDouble(&startTime);
        if (status != 0)
        {
            FatalErrorInFunction
                << "receiveDouble(&startTime) failed! status: " << status
                << abort(FatalError);
        }

        status = receiveDouble(&endTime);
        if (status != 0)
        {
            FatalErrorInFunction
                << "receiveDouble(&endTime) failed! status: " << status
                << abort(FatalError);
        }

        const double currentTime = runTime.value();

        if (mag(currentTime - startTime) > 1e-6)
        {
            FatalErrorInFunction
                << "Times do not match: " << currentTime
                << " vs " << startTime
                << abort(FatalError);
        }
    }

    // Calculate the time step
    double delta = endTime - runTime.value();

    Pstream::scatter(endTime);
    Pstream::scatter(delta);

    // Remove the constness from runTime to modify it during the simulation
    Time& modTime = const_cast<Time&>(runTime);

    // Set the end time and the time step
    modTime.setEndTime(scalar(endTime));
    modTime.setDeltaTNoAdjust(scalar(delta));

    Info<< "Set end time to: " << runTime.endTime().value() << endl;
    Info<< "Adjusted deltaT to: " << runTime.deltaTValue() << endl;
}


void Foam::krakenIO::giveInputFieldNames()
{
    if (Pstream::master())
    {
        label status = sendInt64(inputFields_.size());
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(nInputFields) failed! status: " << status
                << abort(FatalError);
        }

        forAll(inputFields_, i)
        {
            const word& fieldName = inputFields_[i].fieldName;

            status = sendString(fieldName.c_str(), fieldName.length());
            if (status != 0)
            {
                FatalErrorInFunction
                    << "sendString(fieldName) failed! status: " << status
                    << abort(FatalError);
            }
        }
    }
}


void Foam::krakenIO::giveOutputFieldNames()
{
    if (Pstream::master())
    {
        label status = sendInt64(outputFields_.size());
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(nOutputFields) failed! status: " << status
                << abort(FatalError);
        }

        forAll(outputFields_, i)
        {
            const word& fieldName = outputFields_[i].fieldName;

            status = sendString(fieldName.c_str(), fieldName.length());
            if (status != 0)
            {
                FatalErrorInFunction
                    << "sendString(fieldName) failed! status: " << status
                    << abort(FatalError);
            }
        }
    }
}


void Foam::krakenIO::giveFieldTemplate()
{
    if (Pstream::master())
    {
        const word fieldName = receiveString("fieldName");

        Info<< "Sending field template for " << fieldName << endl;

        // Find the matching field
        const fieldTemplate* fieldPtr = nullptr;

        label fieldi = findField(inputFields_, fieldName);

        if (fieldi >= 0)
        {
            fieldPtr = &inputFields_[fieldi];
        }
        else
        {
            fieldi = findField(outputFields_, fieldName);

            if (fieldi < 0)
            {
                FatalErrorInFunction
                    << "Matching field not found for " << fieldName
                    << abort(FatalError);
            }

            fieldPtr = &outputFields_[fieldi];
        }

        const fieldTemplate& field = *fieldPtr;

        // Send the field template
        label status = sendInt64(field.dataType);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(dataType) failed! status: " << status
                << abort(FatalError);
        }

        status = sendInt64Array(field.unitVector, 5);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64Array(unitVector) failed! status: " << status
                << abort(FatalError);
        }

        status = sendDouble(field.multiplier);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendDouble(multiplier) failed! status: " << status
                << abort(FatalError);
        }

        status = sendDouble(field.offset);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendDouble(offset) failed! status: " << status
                << abort(FatalError);
        }

        status = sendString(field.meshName.c_str(), field.meshName.length());
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendString(meshName) failed! status: " << status
                << abort(FatalError);
        }

        status = sendInt64(0);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(remappingFlag) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::giveMeshData()
{
    if (Pstream::master())
    {
        const word meshName = receiveString("meshName");

        label fieldi = -1;

        forAll(inputFields_, i)
        {
            if (inputFields_[i].meshName == meshName)
            {
                fieldi = i;
                break;
            }
        }

        if (fieldi < 0)
        {
            FatalErrorInFunction
                << "Matching mesh not found for " << meshName
                << abort(FatalError);
        }

        Info<< "Sending mesh data for " << meshName << endl;

        label status = sendInt64(MESH_TYPE_LIST);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(meshType) failed! status: " << status
                << abort(FatalError);
        }

        status = sendInt64(inputFields_[fieldi].nCells);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(meshSize) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::giveIndexingArray()
{
    FatalErrorInFunction
        << "Nothing implemented! Should not enter here."
        << abort(FatalError);
}


void Foam::krakenIO::takeFieldData()
{
    label fieldi = -1;

    // Receive the data
    if (Pstream::master())
    {
        const word fieldName = receiveString("fieldName");

        fieldi = findField(inputFields_, fieldName);

        if (fieldi < 0)
        {
            FatalErrorInFunction
                << "No match for field: " << fieldName
                << abort(FatalError);
        }

        fieldTemplate& field = inputFields_[fieldi];

        int64_t nValues;
        label status = receiveInt64(&nValues);
        if (status != 0)
        {
            FatalErrorInFunction
                << "receiveInt64(nValues) failed! status: " << status
                << abort(FatalError);
        }

        if (nValues != field.nCells)
        {
            FatalErrorInFunction
                << "Size mismatch: nValues != nCells " << nl
                << nValues << " != " << field.nCells
                << abort(FatalError);
        }

        // Receive into a double buffer since scalar may not be double
        List<double> dataBuffer(nValues);

        status = receiveDoubleArray(dataBuffer.begin(), nValues);
        if (status != 0)
        {
            FatalErrorInFunction
                << "receiveDoubleArray(dataBuffer) failed! status: " << status
                << abort(FatalError);
        }

        forAll(field.fieldData, i)
        {
            field.fieldData[i] = dataBuffer[i];
        }
    }

    // Scatter the data
    Pstream::scatter(fieldi);
    Pstream::scatter(inputFields_[fieldi].fieldData);
}


void Foam::krakenIO::giveFieldData()
{
    // Receive and scatter the field name
    word fieldName;

    if (Pstream::master())
    {
        fieldName = receiveString("fieldName");
    }

    Pstream::scatter(fieldName);

    const label fieldi = findField(outputFields_, fieldName);

    if (fieldi < 0)
    {
        FatalErrorInFunction
            << "No match for field: " << fieldName
            << abort(FatalError);
    }

    const fieldTemplate& field = outputFields_[fieldi];

    // Update/gather the data to master if needed
    const objectRegistry& region =
        this->db().lookupObject<objectRegistry>(field.regionName);

    fvModels& fvModelsObj = region.lookupObjectRef<fvModels>("fvModels");

    fvConstraints& fvConstraintsObj =
        region.lookupObjectRef<fvConstraints>("fvConstraints");

    const word& modelName = field.modelName;

    if (fvModelsObj.dictionary::found(modelName))
    {
        krakenModelBase& krakenModel =
            dynamic_cast<krakenModelBase&>(fvModelsObj[modelName]);

        krakenModel.updateOutputFields();
    }
    else if (fvConstraintsObj.dictionary::found(modelName))
    {
        krakenModelBase& krakenModel =
            dynamic_cast<krakenModelBase&>(fvConstraintsObj[modelName]);

        krakenModel.updateOutputFields();
    }
    else
    {
        FatalErrorInFunction
            << "No fvModel or fvConstraint with name: " << modelName
            << abort(FatalError);
    }

    // Send the data
    if (Pstream::master())
    {
        label status = sendInt64(field.nCells);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(nValues) failed! status: " << status
                << abort(FatalError);
        }

        // Send from a double buffer since scalar may not be double
        List<double> dataBuffer(field.nCells);

        forAll(field.fieldData, i)
        {
            dataBuffer[i] = field.fieldData[i];
        }

        status = sendDoubleArray(dataBuffer.begin(), field.nCells);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendDoubleArray(dataBuffer) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::giveInputVariableDefs()
{
    if (Pstream::master())
    {
        const label status = sendInt64(nInputVariables_);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(nInputVariables) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::giveOutputVariableDefs()
{
    if (Pstream::master())
    {
        const label status = sendInt64(nOutputVariables_);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(nOutputVariables) failed! status: " << status
                << abort(FatalError);
        }
    }
}


void Foam::krakenIO::takeVariableData()
{
    FatalErrorInFunction
        << "Nothing implemented! Should not enter here."
        << abort(FatalError);
}


void Foam::krakenIO::giveVariableData()
{
    FatalErrorInFunction
        << "Nothing implemented! Should not enter here."
        << abort(FatalError);
}


void Foam::krakenIO::writeRestart()
{
    // Remove the constness from runTime to force writing
    Time& modTime = const_cast<Time&>(this->db().time());

    Info<< "Forcing restart write..." << endl;

    modTime.writeNow();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::krakenIO::krakenIO(const PtrList<solver>& solvers)
:
    IOdictionary
    (
        IOobject
        (
            "krakenIO",
            solvers[0].time().system(),
            solvers[0].time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hostName_(lookupOrDefault<word>("hostName", "localhost")),
    portNo_(lookupOrDefault<label>("port", 2230)),
    isUnixSocket_(lookupOrDefault<bool>("unixSocket", true))
{
    // Create the field templates from the Kraken fvModels and fvConstraints
    createFieldTemplates(solvers);

    if (Pstream::master())
    {
        label status = -100;

        Info<< "Trying to connect to Cerberus... @" << hostName_ << endl;

        // Use negative port number for unix sockets
        const label portSign = isUnixSocket_ ? -1 : 1;

        // Cerberus currently supports port numbers between 2211 - 2230
        for (label port = portNo_; port > 2210; port--)
        {
            status = openConnection(hostName_.c_str(), portSign*port);

            if (status == 0)
            {
                break;
            }
        }
        if (status != 0)
        {
            FatalErrorInFunction
                << "Connection initialization failed! status: " << status
                << abort(FatalError);
        }

        Info<< "Connection established!" << endl;

        // Send a test message
        const string message("OpenFOAM activated!");

        status = sendString(message.c_str(), message.length());
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendString(message) failed! status: " << status
                << abort(FatalError);
        }

        // Receive a test message
        Info<< "Cerberus says: " << receiveString("message") << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::krakenIO::~krakenIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::krakenIO::control(label state)
{
    // Return states
    //  -1  Keep looping = false, Solve flow = false (terminates the execution)
    //   0  Keep looping = true,  Solve flow = false (skips solution phase)
    //   1  Keep looping = true,  Solve flow = true  (proceeds with solution)
    //
    // Input state
    //  The return state of the previous call, used to decide when to send
    //  the solver status to Cerberus
    //
    // Solver status (needs to be implemented properly)
    //   3  Diverged (or major continuity error)
    //   2  Converged
    //   1  Not fully converged (but not diverged either)

    // Send the solver status
    if (state == 1 && Pstream::master())
    {
        const int64_t solverStatus = 2;

        const label status = sendInt64(solverStatus);
        if (status != 0)
        {
            FatalErrorInFunction
                << "sendInt64(solverStatus) failed! status: " << status
                << abort(FatalError);
        }

        Info<< "Sent solver status to Cerberus!" << endl;
    }

    // Get the signal
    label krakenSignal = -1;

    if (Pstream::master())
    {
        int64_t sig;

        const label status = receiveInt64(&sig);
        if (status != 0)
        {
            FatalErrorInFunction
                << "receiveInt64(&sig) failed! status: " << status
                << abort(FatalError);
        }

        krakenSignal = sig;
    }

    Pstream::scatter(krakenSignal);

    Info<< nl << "krakenSignal = " << krakenSignal << endl;

    switch (krakenSignal)
    {
        case KRAKEN_TERMINATE:
        {
            Info<< "KRAKEN_TERMINATE = " << KRAKEN_TERMINATE << endl;
            terminate();
            return -1;
        }
        case KRAKEN_RETURN_CURRENT_TIME:
        {
            Info<< "KRAKEN_RETURN_CURRENT_TIME = "
                << KRAKEN_RETURN_CURRENT_TIME << endl;
            returnCurrentTime();
            return 0;
        }
        case KRAKEN_SET_CURRENT_TIME:
        {
            Info<< "KRAKEN_SET_CURRENT_TIME = "
                << KRAKEN_SET_CURRENT_TIME << endl;
            setCurrentTime();
            return 0;
        }
        case KRAKEN_SUGGEST_NEXT_TIME:
        {
            Info<< "KRAKEN_SUGGEST_NEXT_TIME = "
                << KRAKEN_SUGGEST_NEXT_TIME << endl;
            suggestNextTime();
            return 0;
        }
        case KRAKEN_ADVANCE_TO_TIME:
        {
            Info<< "KRAKEN_ADVANCE_TO_TIME = "
                << KRAKEN_ADVANCE_TO_TIME << endl;
            advanceToTime();
            return 0;
        }
        case KRAKEN_SOLVE_CURRENT_TIME:
        {
            Info<< "KRAKEN_SOLVE_CURRENT_TIME = "
                << KRAKEN_SOLVE_CURRENT_TIME << endl;
            return 1;
        }
        case KRAKEN_GIVE_INPUT_FIELD_NAMES:
        {
            Info<< "KRAKEN_GIVE_INPUT_FIELD_NAMES = "
                << KRAKEN_GIVE_INPUT_FIELD_NAMES << endl;
            giveInputFieldNames();
            return 0;
        }
        case KRAKEN_GIVE_OUTPUT_FIELD_NAMES:
        {
            Info<< "KRAKEN_GIVE_OUTPUT_FIELD_NAMES = "
                << KRAKEN_GIVE_OUTPUT_FIELD_NAMES << endl;
            giveOutputFieldNames();
            return 0;
        }
        case KRAKEN_GIVE_FIELD_TEMPLATE:
        {
            Info<< "KRAKEN_GIVE_FIELD_TEMPLATE = "
                << KRAKEN_GIVE_FIELD_TEMPLATE << endl;
            giveFieldTemplate();
            return 0;
        }
        case KRAKEN_GIVE_MESH_DATA:
        {
            Info<< "KRAKEN_GIVE_MESH_DATA = " << KRAKEN_GIVE_MESH_DATA << endl;
            giveMeshData();
            return 0;
        }
        case KRAKEN_GIVE_INDEXING_ARRAY:
        {
            Info<< "KRAKEN_GIVE_INDEXING_ARRAY = "
                << KRAKEN_GIVE_INDEXING_ARRAY << endl;
            giveIndexingArray();
            return 0;
        }
        case KRAKEN_GIVE_FIELD_DATA:
        {
            Info<< "KRAKEN_GIVE_FIELD_DATA = "
                << KRAKEN_GIVE_FIELD_DATA << endl;
            giveFieldData();
            return 0;
        }
        case KRAKEN_TAKE_FIELD_DATA:
        {
            Info<< "KRAKEN_TAKE_FIELD_DATA = "
                << KRAKEN_TAKE_FIELD_DATA << endl;
            takeFieldData();
            return 0;
        }
        case KRAKEN_GIVE_OUTPUT_VARIABLE_DEFS:
        {
            Info<< "KRAKEN_GIVE_OUTPUT_VARIABLE_DEFS = "
                << KRAKEN_GIVE_OUTPUT_VARIABLE_DEFS << endl;
            giveOutputVariableDefs();
            return 0;
        }
        case KRAKEN_GIVE_VARIABLE_DATA:
        {
            Info<< "KRAKEN_GIVE_VARIABLE_DATA = "
                << KRAKEN_GIVE_VARIABLE_DATA << endl;
            giveVariableData();
            return 0;
        }
        case KRAKEN_GIVE_INPUT_VARIABLE_DEFS:
        {
            Info<< "KRAKEN_GIVE_INPUT_VARIABLE_DEFS = "
                << KRAKEN_GIVE_INPUT_VARIABLE_DEFS << endl;
            giveInputVariableDefs();
            return 0;
        }
        case KRAKEN_TAKE_VARIABLE_DATA:
        {
            Info<< "KRAKEN_TAKE_VARIABLE_DATA = "
                << KRAKEN_TAKE_VARIABLE_DATA << endl;
            takeVariableData();
            return 0;
        }
        case KRAKEN_RETURN_CURRENT_TIME_INTERVAL:
        {
            Info<< "KRAKEN_RETURN_CURRENT_TIME_INTERVAL = "
                << KRAKEN_RETURN_CURRENT_TIME_INTERVAL << endl;
            returnCurrentTimeInterval();
            return 0;
        }
        case KRAKEN_SET_CURRENT_TIME_INTERVAL:
        {
            Info<< "KRAKEN_SET_CURRENT_TIME_INTERVAL = "
                << KRAKEN_SET_CURRENT_TIME_INTERVAL << endl;
            setCurrentTimeInterval();
            return 0;
        }
        case KRAKEN_SUGGEST_NEXT_TIME_INTERVAL:
        {
            Info<< "KRAKEN_SUGGEST_NEXT_TIME_INTERVAL = "
                << KRAKEN_SUGGEST_NEXT_TIME_INTERVAL << endl;
            suggestNextTimeInterval();
            return 0;
        }
        case KRAKEN_ADVANCE_TO_TIME_INTERVAL:
        {
            Info<< "KRAKEN_ADVANCE_TO_TIME_INTERVAL = "
                << KRAKEN_ADVANCE_TO_TIME_INTERVAL << endl;
            advanceToTimeInterval();
            return 0;
        }
        case KRAKEN_WRITE_RESTART:
        {
            Info<< "KRAKEN_WRITE_RESTART = " << KRAKEN_WRITE_RESTART << endl;
            writeRestart();
            return 0;
        }
        default:
        {
            FatalErrorInFunction
                << "Unrecognised signal: " << krakenSignal
                << abort(FatalError);
        }
    }

    return -100;
}


const Foam::List<Foam::scalar>& Foam::krakenIO::getFieldData
(
    const word& fieldName
) const
{
    const label fieldi = findField(inputFields_, fieldName);

    if (fieldi < 0)
    {
        wordList fieldNames(inputFields_.size());

        forAll(inputFields_, i)
        {
            fieldNames[i] = inputFields_[i].fieldName;
        }

        FatalErrorInFunction
            << "No match for field: " << fieldName << nl
            << "Available input fields: " << fieldNames
            << abort(FatalError);
    }

    return inputFields_[fieldi].fieldData;
}


Foam::List<Foam::scalar>& Foam::krakenIO::setFieldData
(
    const word& fieldName
)
{
    const label fieldi = findField(outputFields_, fieldName);

    if (fieldi < 0)
    {
        wordList fieldNames(outputFields_.size());

        forAll(outputFields_, i)
        {
            fieldNames[i] = outputFields_[i].fieldName;
        }

        FatalErrorInFunction
            << "No match for field: " << fieldName << nl
            << "Available output fields: " << fieldNames
            << abort(FatalError);
    }

    return outputFields_[fieldi].fieldData;
}


// ************************************************************************* //
