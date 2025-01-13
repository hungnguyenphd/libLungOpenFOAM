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

#include "airInletVOFVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"
#include "patchPatchDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::airInletVOFVelocityFvPatchVectorField::
airInletVOFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRateFraction_(nullptr),
    rhoName_("rho"),
    rhoOutlet_(0),
    volumetric_(false),
    extrapolateProfile_(false),
    parabolic_(false),
    y_(p.size(), one())
{}


Foam::airInletVOFVelocityFvPatchVectorField::
airInletVOFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    flowRateFraction_(nullptr),
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
    flowRateFraction_ =
        Function1<scalar>::NewIfPresent("flowRateFraction", dict, &db());

    if (flowRateFraction_)
    {
        volumetric_ = true;
    }
    else
    {
        dict.readIfPresent("rho", rhoName_);
        flowRateFraction_ =
            Function1<scalar>::NewIfPresent("massFlowRate", dict, &db());
    }

    dict.readIfPresent("parabolic", parabolic_);
    if (parabolic_)
    {
        setWallDist();
    }

    if (!flowRateFraction_)
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


Foam::airInletVOFVelocityFvPatchVectorField::
airInletVOFVelocityFvPatchVectorField
(
    const airInletVOFVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRateFraction_(ptf.flowRateFraction_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


Foam::airInletVOFVelocityFvPatchVectorField::
airInletVOFVelocityFvPatchVectorField
(
    const airInletVOFVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRateFraction_(ptf.flowRateFraction_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


Foam::airInletVOFVelocityFvPatchVectorField::
airInletVOFVelocityFvPatchVectorField
(
    const airInletVOFVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRateFraction_(ptf.flowRateFraction_.clone()),
    rhoName_(ptf.rhoName_),
    rhoOutlet_(ptf.rhoOutlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_),
    parabolic_(ptf.parabolic_),
    y_(ptf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::airInletVOFVelocityFvPatchVectorField::setWallDist()
{
    boundBox bb_(patch().patch().localPoints(), true);
    //Foam::Info << patch().patch().localPoints() << endl;
    //Foam::Info << bb_.max() << endl;
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& c_ = patch().Cf();
    scalarField rp_ = mag(c_ - ctr_);

    const scalarField rr_(rp_/gMax(rp_));

    //Foam::Info << "Center of patch inlet" << ctr_ << endl;
    //Foam::Info << rp_ << endl;

    y_ = 2*(1 - sqr(rr_));

    area_ = gSum(patch().magSf());
}

// Foam::tmp<Foam::scalarField>
// Foam::airInletVOFVelocityFvPatchVectorField::profile()
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
void Foam::airInletVOFVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{

    const List<double> time_data{0, 0.002251363, 0.005598246, 0.009234403, 0.013235621, 0.016101601, 0.019443466, 
                0.02325706, 0.027070209, 0.03135574, 0.035638334, 0.039917724, 0.044666115, 
                0.049888727, 0.055576482, 0.061260678, 0.066938049, 0.072611267, 0.078279738, 
                0.083941089, 0.089600363, 0.095252516, 0.100904669, 0.106550889, 0.112195328, 
                0.117836801, 0.123475604, 0.129113219, 0.134748758, 0.140384, 0.146017759, 
                0.151651518, 0.157285574, 0.162917849, 0.168551311, 0.17418596, 0.179820312, 
                0.185454862, 0.191091538, 0.196729006, 0.202366622, 0.208006017, 0.213646007, 
                0.220227082, 0.225868851, 0.231511214, 0.237155357, 0.24279861, 0.248443346, 
                0.254088379, 0.259733411, 0.265379038, 0.271024961, 0.27667029, 0.282315916, 
                0.287961543, 0.293606872, 0.299251905, 0.304896641, 0.310541377, 0.31618552, 
                0.321829663, 0.327473509, 0.333117355, 0.338760905, 0.344404454, 0.350048597, 
                0.355692443, 0.361336883, 0.366981322, 0.372626652, 0.378271685, 0.383918201, 
                0.389565014, 0.395212717, 0.400860717, 0.406510497, 0.412160276, 0.417810946, 
                0.4234631, 0.42911644, 0.434768593, 0.44042223, 0.44607735, 0.451731283, 
                0.45738492, 0.463038557, 0.468692194, 0.47434583, 0.479998577, 0.48564895, 
                0.491298137, 0.496475, 0.502595175, 0.507769813};

    const List<double> flowrate_data_Ls{0, 0.14480054, 0.437884497, 0.770442972, 1.032092584, 1.266687654, 1.529724963, 
                1.79633832, 2.060287442, 2.331720156, 2.585568921, 2.820235195, 3.042136877, 
                3.282534322, 3.48707714, 3.670306079, 3.812683419, 3.930194567, 4.019287211, 
                4.065752099, 4.099783891, 4.091187926, 4.082591962, 4.038472866, 3.983696832, 
                3.911159233, 3.822636225, 3.727008591, 3.618947861, 3.509110975, 3.390393306, 
                3.271675637, 3.154734124, 3.027135673, 2.906641847, 2.793252648, 2.678087292, 
                2.56410604, 2.462853911, 2.366338198, 2.270710564, 2.185739869, 2.104321487, 
                2.015845949, 1.945084506, 1.877875376, 1.821323185, 1.759442525, 1.706442647, 
                1.655218926, 1.603995205, 1.556323797, 1.510428545, 1.460980981, 1.413309573, 
                1.365638165, 1.3161906, 1.264966879, 1.211967001, 1.158967123, 1.102414933, 
                1.045862742, 0.987534395, 0.929206048, 0.869101544, 0.80899704, 0.752444849, 
                0.694116502, 0.639340468, 0.584564434, 0.535116869, 0.483893148, 0.44155021, 
                0.400983428, 0.365745115, 0.332282959, 0.309477742, 0.286672526, 0.269195778, 
                0.260599814, 0.259108475, 0.250512511, 0.250797328, 0.259962929, 0.262023904, 
                0.262308722, 0.26259354, 0.262878357, 0.263163175, 0.258119524, 0.23886662, 
                0.21250909, 0.1790232, 0.151801326, 0.104994262};

    const List<double> flowrate_data = flowrate_data_Ls/1000;

    const scalar t = db().time().timeOutputValue();

    int lower_id = 0;
    int upper_id = 0;
    //int i;
    forAll(time_data, i)
    {
        if (t >= time_data[i])
        {
            lower_id = i;
        }
    }

    forAll(time_data, i)
    {
        if (t <= time_data[i])
        {
            upper_id = i;
            break;
        }
    }

    const scalar time_weight = (t - time_data[lower_id])
                    /(time_data[upper_id] - time_data[lower_id]);

    const scalar flowrate = time_weight*(flowrate_data[upper_id]-flowrate_data[lower_id]) + flowrate_data[lower_id];

    const scalarField profile(y_);

    //Foam::Info << patch().patch().name() << " " << profile << endl;

    const scalar outletFlowRate = flowrate*flowRateFraction_->value(t);

    const scalar avgU = outletFlowRate/gSum(rho*patch().magSf());

    //Normal vector to patch
    const vectorField n(patch().nf());

    //vectorField Up(avgU*profile*n);
    vectorField Up(avgU*profile*n);

    scalarField nUp(n & Up);

    const scalar estimatedFlowRate = gSum(rho*(this->patch().magSf()*nUp));
    //Foam::Info << patch().patch().name() << " estimatedFlowRate " << estimatedFlowRate << endl;

    const scalar ratio = mag(estimatedFlowRate)/mag(outletFlowRate);
    //Foam::Info << "ratio = " << ratio << endl;
    if (ratio > 0.5)
    {
        nUp *= (mag(outletFlowRate)/mag(estimatedFlowRate));
    }
    else
    {
        nUp += ((outletFlowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
    }

    // Add the corrected normal component of velocity to the patch velocity
    Up = nUp*n;

    //const scalar correctedFlowRate = gSum(rho*(this->patch().magSf()*(n & Up)));
    //Foam::Info << "correctedFlowRate " << correctedFlowRate << endl;

    // Correct the patch velocity
    this->operator==(Up);
}


void Foam::airInletVOFVelocityFvPatchVectorField::updateCoeffs()
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


void Foam::airInletVOFVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRateFraction_->writeData(os);
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
       airInletVOFVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
