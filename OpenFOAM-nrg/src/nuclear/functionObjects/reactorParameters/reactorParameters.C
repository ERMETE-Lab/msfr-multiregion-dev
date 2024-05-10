/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "reactorParameters.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reactorParameters, 0);
    addToRunTimeSelectionTable(functionObject, reactorParameters, dictionary);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::reactorParameters::modeType,
    2
>::names[] = {"magnitude", "component"};

const Foam::NamedEnum
<
    Foam::functionObjects::reactorParameters::modeType,
    2
> Foam::functionObjects::reactorParameters::modeTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reactorParameters::writeFileHeader(const label i)
{
    OFstream& file = this->file();

    writeHeader(file, "Reactor parameters");
    writeCommented(file, "Time");

    forAll(fieldSet_, fieldi)
    {
        writeTabbed(file, fieldSet_[fieldi]);
    }

    file<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reactorParameters::reactorParameters
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    mode_(mdMag),
    fieldSet_()
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::reactorParameters::~reactorParameters()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reactorParameters::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::reactorParameters::execute()
{
    return true;
}


bool Foam::functionObjects::reactorParameters::write()
{
    logFiles::write();

    if (Pstream::master()) writeTime(file());
    Log << type() << " " << name() <<  " write:" << nl;

    forAll(fieldSet_, fieldi)
    {
        calcReactorParameters<scalar>(fieldSet_[fieldi], mdCmpt);
        calcReactorParameters<vector>(fieldSet_[fieldi], mode_);
        calcReactorParameters<sphericalTensor>(fieldSet_[fieldi], mode_);
        calcReactorParameters<symmTensor>(fieldSet_[fieldi], mode_);
        calcReactorParameters<tensor>(fieldSet_[fieldi], mode_);
    }

    if (Pstream::master()) file() << endl;
    Log << endl;

    return true;
}



// ************************************************************************* //
