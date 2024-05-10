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

#include "albedoOldFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::albedoOldFvPatchField<Type>::albedoOldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gamma_(p.size()),
    diffCoeff_(p.size())
    //diffCoeffName_("diffCoeffName")
{}


template<class Type>
Foam::albedoOldFvPatchField<Type>::albedoOldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    gamma_("gamma", dict, p.size()),
    diffCoeff_("diffCoeff", dict, p.size())
    //diffCoeffName_(dict.lookupOrDefault<word>("diffCoeffName", "diffCoeffName"))
{
    evaluate();
}


template<class Type>
Foam::albedoOldFvPatchField<Type>::albedoOldFvPatchField
(
    const albedoOldFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }

    mapper(gamma_, ptf.gamma_);
    mapper(diffCoeff_, ptf.diffCoeff_);
}


template<class Type>
Foam::albedoOldFvPatchField<Type>::albedoOldFvPatchField
(
    const albedoOldFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gamma_(ptf.gamma_),
    diffCoeff_(ptf.diffCoeff_)
    //diffCoeffName_(ptf.diffCoeffName_)
{}


template<class Type>
Foam::albedoOldFvPatchField<Type>::albedoOldFvPatchField
(
    const albedoOldFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gamma_(ptf.gamma_),
    diffCoeff_(ptf.diffCoeff_)
    //diffCoeffName_(ptf.diffCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::albedoOldFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    m(gamma_, gamma_);
    m(diffCoeff_, diffCoeff_); // andrà aggiunta??? nella versione di Stefano et al. non c'è
}


template<class Type>
void Foam::albedoOldFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const albedoOldFvPatchField<Type>& fgptf =
        refCast<const albedoOldFvPatchField<Type>>(ptf);

    gamma_.rmap(fgptf.gamma_, addr);
    diffCoeff_.rmap(fgptf.diffCoeff_, addr);
}


template<class Type>
void Foam::albedoOldFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() - (this->patchInternalField()*gamma_/diffCoeff_)/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoOldFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>(new Field<Type>(this->size(), pTraits<Type>::zero));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoOldFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>(new Field<Type>(this->size(), pTraits<Type>::zero));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoOldFvPatchField<Type>::gradientInternalCoeffs() const
{
    return Type(pTraits<Type>::one)*(- gamma_/diffCoeff_);
    /*return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );*/
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoOldFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoOldFvPatchField<Type>::snGrad() const
{
    return (- (this->patchInternalField()*gamma_/diffCoeff_)/this->patch().deltaCoeffs());
}

template<class Type>
void Foam::albedoOldFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "diffCoeff", diffCoeff_);

    //writeEntry("value", os);
    os.writeKeyword("value uniform 0") << token::END_STATEMENT << nl;
}


// ************************************************************************* //
