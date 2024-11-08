/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "flowRateInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"
#include "patchPatchDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletOutletVelocityFvPatchVectorField::
flowRateInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(nullptr),
    rhoName_("rho"),
    rhoOutlet_(0),
    volumetric_(false),
    extrapolateProfile_(false),
    parabolic_(false),
    y_(p.size(), one())
{}


Foam::flowRateInletOutletVelocityFvPatchVectorField::
flowRateInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    flowRate_(nullptr),
    rhoName_("rho"),
    rhoOutlet_(dict.getOrDefault<scalar>("rhoOutlet", -VGREAT)),
    volumetric_(false),
    extrapolateProfile_
    (
        dict.getOrDefault<Switch>("extrapolateProfile", false)
    ),
    parabolic_(false),
    y_(p.size(), one())
{
    flowRate_ =
        Function1<scalar>::NewIfPresent("volumetricFlowRate", dict, &db());

    if (flowRate_)
    {
        volumetric_ = true;
    }
    else
    {
        dict.readIfPresent("rho", rhoName_);
        flowRate_ =
            Function1<scalar>::NewIfPresent("massFlowRate", dict, &db());
    }

    dict.readIfPresent("parabolic", parabolic_);
    if (parabolic_)
    {
        setWallDist();
    }

    if (!flowRate_)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' (optional: with 'rho')" << nl
            << exit(FatalIOError);
    }

    // Value field required if mass based
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::buffered);
    }
}


Foam::flowRateInletOutletVelocityFvPatchVectorField::
flowRateInletOutletVelocityFvPatchVectorField
(
    const flowRateInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


Foam::flowRateInletOutletVelocityFvPatchVectorField::
flowRateInletOutletVelocityFvPatchVectorField
(
    const flowRateInletOutletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


Foam::flowRateInletOutletVelocityFvPatchVectorField::
flowRateInletOutletVelocityFvPatchVectorField
(
    const flowRateInletOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::flowRateInletOutletVelocityFvPatchVectorField::setWallDist()
{
    const labelHashSet otherPatchIDs
    (
        patch().patch().boundaryMesh().findPatchIDs<polyPatch>()
    ); 

    Foam::Info << "para_mucus" << endl;
    // Foam::Info << patch().patch().name() << " otherPatchIDs " << otherPatchIDs.sortedToc() << endl;
    // //Foam::Info << patch().patch().name() << " " << patch().patch().boundaryMesh().mesh().boundaryConnections() << endl;
    const patchPatchDist pwd(patch().patch(), otherPatchIDs);

    const scalarField r_(pwd/gMax(pwd));
    // //Foam::Info << "max(r) " << gMax(pwd) << endl;

    // boundBox bb_(patch().patch().localPoints(), true);
    // vector ctr_ = 0.5*(bb_.max() + bb_.min());
    // const vectorField& c_ = patch().Cf();
    // scalarField rp_ = mag(c_ - ctr_);

    // const scalarField rr_(rp_/gMax(rp_));

    //Foam::Info << "Center of patch inlet" << ctr_ << endl;
    //Foam::Info << rp_ << endl;

    //y_ = 2*(1 - sqr(rr_));
    y_ = 2*(1 - sqr(1-mag(r_)));

    area_ = gSum(patch().magSf());
}

// Foam::tmp<Foam::scalarField>
// Foam::flowRateInletOutletVelocityFvPatchVectorField::profile()
// {
//     if (profile_.valid())
//     {
//         return profile_->value(y_);
//     }
//     else
//     {
//         return tmp<scalarField>(new scalarField(size(), scalar(1)));
//     }
// }

template<class RhoType>
void Foam::flowRateInletOutletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{

    const scalarField profile(this->y_);

    //Foam::Info << patch().patch().name() << " " << profile << endl;

    const scalar t = db().time().timeOutputValue();

    const scalar flowRate = flowRate_->value(t);

    const scalar avgU = flowRate/gSum(rho*patch().magSf());

    //Normal vector to patch
    const vectorField n(patch().nf());

    //vectorField Up(avgU*profile*n);
    vectorField Up(avgU*profile*n);

    scalarField nUp(n & Up);

    const scalar estimatedFlowRate = gSum(rho*(this->patch().magSf()*nUp));
    //Foam::Info << patch().patch().name() << " estimatedFlowRate " << estimatedFlowRate << endl;

    const scalar ratio = mag(estimatedFlowRate)/mag(flowRate);

    if (ratio > 0.5)
    {
        nUp *= (mag(flowRate)/mag(estimatedFlowRate));
    }
    else
    {
        nUp += ((flowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
    }

    // Add the corrected normal component of velocity to the patch velocity
    Up = nUp*n;
    
    const scalar correctedFlowRate = gSum(rho*(this->patch().magSf()*(n & Up)));
    //Foam::Info << "correctedFlowRate " << correctedFlowRate << endl;

    // Correct the patch velocity
    this->operator==(Up);
}


void Foam::flowRateInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one{});
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const auto& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoOutlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoOutlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoOutlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    os.writeEntry<bool>("parabolic", parabolic_);

    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoOutlet", -VGREAT, rhoOutlet_);
    }
    if (extrapolateProfile_)
    {
        os.writeEntry("extrapolateProfile", extrapolateProfile_);
    }
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletOutletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
