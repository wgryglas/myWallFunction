/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "myEpsilonWallFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#define ASSERT_WALL_PATCH\
        if (!isA<wallFvPatch>(patch())) \
        { \
            FatalErrorIn("constructor of myEpsilonWallFunction") \
                << "Invalid wall function specification" << nl \
                << "    Patch type for patch " << patch().name() \
                << " must be wall" << nl \
                << "    Current patch type is " << patch().type() << nl << endl \
                << abort(FatalError);\
        }


#define ASSERT_MASTER\
    if(!master_) return;

#define ASSERT_INIT\
    if(!init_) \
    { \
       ASSERT_INIT \
       calacNBFaces(); \
       initG(); \
       init_=true;\
    }


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * * * //
void Foam::myEpsilonWallFunction::calcNBFaces()
{
    const volScalarField& epsilon = static_cast<const volScalarField&>( this->dimensionedInternalField() );
    const volScalarField::GeometricBoundaryField& epsBF = epsilon.boundaryField();

    DynamicList<label> dNb;

    forAll(epsBF, patchId)
    {
        const fvPatchField& pf = epsBF[patchId];
        const fvPatch & p = pf.patch();

        if(isA<myEpsilonWallFunction>(pf))
        {
            const labelUList& faceCells = p.faceCells();

            forAll(faceCells, i)
            {
                label cellI = faceCells[i];

                Map<label>::iterator itr = cellId2Local.find(cellI);
                if( itr != cellId2Local.end())
                {
                        dNb[*itr]++;
                }
                else
                {
                    cellId2Local.insert(cellI,dNb.size());
                    dNb.append(1);
                }
            }
        }
    }

    nbBFaces=dNb;
}

void Foam::myEpsilonWallFunction::initG()
{
    G_.setSize(nbBFaces.size(),0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myEpsilonWallFunction::myEpsilonWallFunction
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    master_(false)
{
    ASSERT_WALL_PATCH
}


Foam::myEpsilonWallFunction::myEpsilonWallFunction
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    master_(false)
{
    ASSERT_WALL_PATCH
}


Foam::myEpsilonWallFunction::myEpsilonWallFunction
(
    const myEpsilonWallFunction& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    master_(false)
{
    ASSERT_WALL_PATCH
}


Foam::myEpsilonWallFunction::myEpsilonWallFunction
(
    const myEpsilonWallFunction& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    master_(false)
{
    ASSERT_WALL_PATCH
}


Foam::myEpsilonWallFunction::myEpsilonWallFunction
(
    const myEpsilonWallFunction& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    master_(false)
{
    ASSERT_WALL_PATCH
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::myEpsilonWallFunction::updateCoeffs()
{
    ASSERT_INIT
    // Apply epsilon value to bc
    db().



    ASSERT_MASTER
    //Modify G field


}


void Foam::myEpsilonWallFunction::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "G", "G", GName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myEpsilonWallFunction
    );
}

// ************************************************************************* //
