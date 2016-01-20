/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::myEpsilonWallFunction

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    !todo

Note
    The default boundary behaviour is for subsonic, incompressible flow.

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    myEpsilonWallFunction.C

\*---------------------------------------------------------------------------*/

#ifndef myEpsilonWallFunction_H
#define myEpsilonWallFunction_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class myEpsilonWallFunction Declaration
\*---------------------------------------------------------------------------*/

class myEpsilonWallFunction
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Weigths defining boundary influence into cell
        //  (if corner more then on cell face influence
        //   the value in the wall boundary region)
        labelList nbBFaces;
        Map<label> cellId2Local;

        //Is it the main patch which does wall cells specific calcualtions
        bool master_;

        //  Turbulecne parameters
        scalar Cmu_;
        scalar E_;
        scalar kappa_;

        bool init_;

    // Member helper functions
        void calcNBFaces();

        void initG();


public:

    //- Runtime type information
    TypeName("myEpsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        myEpsilonWallFunction
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myEpsilonWallFunction
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given myEpsilonWallFunction
        //  onto a new patch
        myEpsilonWallFunction
        (
            const myEpsilonWallFunction&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myEpsilonWallFunction
        (
            const myEpsilonWallFunction&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallFunction(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myEpsilonWallFunction
        (
            const myEpsilonWallFunction&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myEpsilonWallFunction(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& GName() const
            {
                return GName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            virtual void manipulateMatrix(fvMatrix<scalar>& matrix, const scalarField& weights);



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
