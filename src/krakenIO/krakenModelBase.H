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
    Foam::krakenModelBase

Description
    Base class for fvModels and fvConstraints in Kraken simulations.

    Derived classes define the fields exchanged with Cerberus via
    createFieldTemplates and update the output fields on request via
    updateOutputFields.

SourceFiles
    krakenModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef krakenModelBase_H
#define krakenModelBase_H

#include "krakenIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class krakenModelBase Declaration
\*---------------------------------------------------------------------------*/

class krakenModelBase
{
    // Private Data

        //- Flag indicating that the output fields need to be updated
        mutable bool outputUpdateRequired_;


public:

    // Constructors

        //- Construct null
        krakenModelBase();

        //- Disallow default bitwise copy construction
        krakenModelBase(const krakenModelBase&) = delete;


    //- Destructor
    virtual ~krakenModelBase();


    // Member Functions

        //- Create krakenIO field templates
        virtual void createFieldTemplates
        (
            List<krakenIO::fieldTemplate>& inputFields,
            List<krakenIO::fieldTemplate>& outputFields
        ) const = 0;

        //- Update krakenIO output fields to master if needed
        virtual void updateOutputFields() = 0;

        //- Set output update flag
        void setOutputUpdateRequired() const
        {
            outputUpdateRequired_ = true;
        }

        //- Clear output update flag
        void clearOutputUpdateRequired() const
        {
            outputUpdateRequired_ = false;
        }

        //- Return whether output fields need to be updated
        bool outputUpdateRequired() const
        {
            return outputUpdateRequired_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const krakenModelBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
