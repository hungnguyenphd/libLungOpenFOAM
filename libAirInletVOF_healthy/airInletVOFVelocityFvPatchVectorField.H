/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::airInletVOFVelocityFvPatchVectorField

Group
    grpOutletBoundaryConditions Hung

Description
    Velocity outlet boundary condition which corrects the extrapolated velocity
    to match the specified flow rate.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m^3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the outlet density using the \c rhoOutlet entry

    For a volumetric-based flux:
    - the flow rate is in m^3/s

Usage
    \table
        Property     | Description             | Required    | Default value
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m^3/s]| no |
        rho          | density field name      | no          | rho
        rhoOutlet    | outlet density          | no          |
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                airInletVOFVelocity;
        volumetricFlowRate  0.2;
        value               uniform (0 0 0);
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                airInletVOFVelocity;
        massFlowRate        0.2;
        rhoOutlet           1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoOutlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive out of the domain (as an outlet)
    - May not work correctly for transonic outlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::flowRateInletVelocityFvPatchVectorField

SourceFiles
    airInletVOFVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef airInletVOFVelocityFvPatchVectorField_H
#define airInletVOFVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "coupleGroupIdentifier.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class airInletVOFVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class airInletVOFVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Outlet integral flow rate
        autoPtr<Function1<scalar>> flowRateFraction_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoOutlet_;

        //- Is volumetric?
        bool volumetric_;

        //- Set true to extrapolate the velocity profile from the interior
        Switch extrapolateProfile_;
        
        //- Is parabolic?
        bool parabolic_;
        
        //- Normalised distance to the wall
        scalarField y_;

        //- Total area of patch
        scalar area_;

    // Private Member Functions

        //- Set the wall distance field y_
        void setWallDist();

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("airInletVOFVelocity");


   // Constructors

        //- Construct from patch and internal field
        airInletVOFVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        airInletVOFVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  airInletVOFVelocityFvPatchVectorField
        //  onto a new patch
        airInletVOFVelocityFvPatchVectorField
        (
            const airInletVOFVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        airInletVOFVelocityFvPatchVectorField
        (
            const airInletVOFVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        airInletVOFVelocityFvPatchVectorField
        (
            const airInletVOFVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
