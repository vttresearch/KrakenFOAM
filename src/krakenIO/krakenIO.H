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

Class
    Foam::krakenIO

Description
    Socket communication with the Cerberus multiphysics driver in Kraken
    simulations.

    Reads the connection settings from the \c system/krakenIO dictionary,
    collects the input and output field templates from all Kraken fvModels
    and fvConstraints of the given solvers and responds to the signals sent
    by Cerberus.

Usage
    Example of the \c system/krakenIO dictionary:
    \verbatim
        hostName    localhost;
        port        2230;
        unixSocket  true;
    \endverbatim

    \table
        Property    | Description                         | Required | Default
        hostName    | Name of the Cerberus host           | no       | localhost
        port        | Port number                         | no       | 2230
        unixSocket  | Use a unix socket instead of inet   | no       | true
    \endtable

SourceFiles
    krakenIO.C

\*---------------------------------------------------------------------------*/

#ifndef krakenIO_H
#define krakenIO_H

#include "IOdictionary.H"
#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class krakenIO Declaration
\*---------------------------------------------------------------------------*/

class krakenIO
:
    public IOdictionary
{
public:

    // Public Classes

        //- Template describing a field exchanged with Cerberus
        struct fieldTemplate
        {
            //- Name of the field
            word fieldName;

            //- Name of the mesh the field is associated with
            word meshName;

            //- Name of the fvModel or fvConstraint providing the field
            word modelName;

            //- Name of the mesh region
            word regionName;

            //- Number of values in the field
            int64_t nCells;

            //- Cerberus data type
            int64_t dataType;

            //- Exponents of the SI base units [kg m s K mol]
            int64_t unitVector[5];

            //- Unit conversion multiplier
            double multiplier;

            //- Unit conversion offset
            double offset;

            //- Field values
            List<scalar> fieldData;
        };


private:

    // Private Enumerations

        //- Signals sent by Cerberus
        enum KRAKEN_SIGNALS
        {
            KRAKEN_TERMINATE,
            KRAKEN_RETURN_CURRENT_TIME,
            KRAKEN_SET_CURRENT_TIME,
            KRAKEN_SUGGEST_NEXT_TIME,
            KRAKEN_ADVANCE_TO_TIME,
            KRAKEN_SOLVE_CURRENT_TIME,
            KRAKEN_GIVE_INPUT_FIELD_NAMES,
            KRAKEN_GIVE_OUTPUT_FIELD_NAMES,
            KRAKEN_GIVE_FIELD_TEMPLATE,
            KRAKEN_GIVE_MESH_DATA,
            KRAKEN_GIVE_INDEXING_ARRAY,
            KRAKEN_GIVE_FIELD_DATA,
            KRAKEN_TAKE_FIELD_DATA,
            KRAKEN_GIVE_OUTPUT_VARIABLE_DEFS,
            KRAKEN_GIVE_VARIABLE_DATA,
            KRAKEN_GIVE_INPUT_VARIABLE_DEFS,
            KRAKEN_TAKE_VARIABLE_DATA,
            KRAKEN_WRITE_RESTART,
            KRAKEN_READ_RESTART,
            KRAKEN_RETURN_CURRENT_TIME_INTERVAL,
            KRAKEN_SET_CURRENT_TIME_INTERVAL,
            KRAKEN_SUGGEST_NEXT_TIME_INTERVAL,
            KRAKEN_ADVANCE_TO_TIME_INTERVAL,
            KRAKEN_N_SIGNALS
        };


    // Private Data

        //- Cerberus mesh types
        const label MESH_TYPE_REGULAR_CARTESIAN = 1;
        const label MESH_TYPE_IRREGULAR_CARTESIAN = 2;
        const label MESH_TYPE_NESTED = 8;
        const label MESH_TYPE_LIST = 9;

        //- Maximum length of the strings received from Cerberus
        static const label MAX_STRING_SIZE = 100;

        //- Host name
        const word hostName_;

        //- Port number
        const label portNo_;

        //- True if a unix socket is used, false if inet
        const bool isUnixSocket_;

        //- Input field templates
        List<fieldTemplate> inputFields_;

        //- Output field templates
        List<fieldTemplate> outputFields_;

        //- Number of input variables
        const label nInputVariables_ = 0;

        //- Number of output variables
        const label nOutputVariables_ = 0;


    // Private Member Functions

        //- Return the index of the named field in the given list
        //  or -1 if not found
        static label findField
        (
            const List<fieldTemplate>& fields,
            const word& fieldName
        );

        //- Receive a string from Cerberus
        string receiveString(const char* name) const;

        //- Create the field templates from the Kraken fvModels and
        //  fvConstraints of the given solvers
        void createFieldTemplates(const PtrList<solver>& solvers);


        // Cerberus communication

            void terminate();
            void returnCurrentTime();
            void returnCurrentTimeInterval();
            void setCurrentTime();
            void setCurrentTimeInterval();
            void suggestNextTime();
            void suggestNextTimeInterval();
            void advanceToTime();
            void advanceToTimeInterval();
            void giveInputFieldNames();
            void giveOutputFieldNames();
            void giveFieldTemplate();
            void giveMeshData();
            void giveIndexingArray();
            void giveFieldData();
            void takeFieldData();
            void giveOutputVariableDefs();
            void giveVariableData();
            void giveInputVariableDefs();
            void takeVariableData();
            void writeRestart();


public:

    // Constructors

        //- Construct from the list of solvers and connect to Cerberus
        krakenIO(const PtrList<solver>& solvers);

        //- Disallow default bitwise copy construction
        krakenIO(const krakenIO&) = delete;


    //- Destructor
    ~krakenIO();


    // Member Functions

        // Control

            //- Receive and execute the next signal from Cerberus
            //  Returns:
            //      -1: terminate the execution
            //       0: keep looping without solving
            //       1: solve the current time interval
            label control(label state);


        // Access

            //- Return the data of the named input field
            const List<scalar>& getFieldData(const word& fieldName) const;

            //- Return the data of the named output field for modification
            List<scalar>& setFieldData(const word& fieldName);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const krakenIO&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
