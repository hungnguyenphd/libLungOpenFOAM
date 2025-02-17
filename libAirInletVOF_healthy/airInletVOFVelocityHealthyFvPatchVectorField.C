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

#include "airInletVOFVelocityHealthyFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"
#include "patchPatchDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::airInletVOFVelocityHealthyFvPatchVectorField::
airInletVOFVelocityHealthyFvPatchVectorField
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


Foam::airInletVOFVelocityHealthyFvPatchVectorField::
airInletVOFVelocityHealthyFvPatchVectorField
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


Foam::airInletVOFVelocityHealthyFvPatchVectorField::
airInletVOFVelocityHealthyFvPatchVectorField
(
    const airInletVOFVelocityHealthyFvPatchVectorField& ptf,
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


Foam::airInletVOFVelocityHealthyFvPatchVectorField::
airInletVOFVelocityHealthyFvPatchVectorField
(
    const airInletVOFVelocityHealthyFvPatchVectorField& ptf
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


Foam::airInletVOFVelocityHealthyFvPatchVectorField::
airInletVOFVelocityHealthyFvPatchVectorField
(
    const airInletVOFVelocityHealthyFvPatchVectorField& ptf,
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
void Foam::airInletVOFVelocityHealthyFvPatchVectorField::setWallDist()
{
    const labelHashSet otherPatchIDs
    (
        patch().patch().boundaryMesh().findPatchIDs<wallPolyPatch>()
    );

    const patchPatchDist pwd(patch().patch(), otherPatchIDs);

    const scalarField r_(pwd/gMax(pwd));

    boundBox bb_(patch().patch().localPoints(), true);
    //Foam::Info << patch().patch().localPoints() << endl;
    //Foam::Info << bb_.max() << endl;
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& c_ = patch().Cf();
    scalarField rp_ = mag(c_ - ctr_);

    const scalarField rr_(rp_/gMax(rp_));

    //Foam::Info << "Center of patch inlet" << ctr_ << endl;
    //Foam::Info << rp_ << endl;

    //y_ = 2*(1 - sqr(rr_));
    y_ = 2*(1 - sqr(r_));

    area_ = gSum(patch().magSf());
}

// Foam::tmp<Foam::scalarField>
// Foam::airInletVOFVelocityHealthyFvPatchVectorField::profile()
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
void Foam::airInletVOFVelocityHealthyFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{

    const List<double> time_data{0.        , 0.0021978 , 0.00439559, 0.00659339, 0.00879118,
       0.01098898, 0.01318678, 0.01538457, 0.01758237, 0.01978016,
       0.02197796, 0.02417576, 0.02637355, 0.02857135, 0.03076914,
       0.03296694, 0.03516474, 0.03736253, 0.03956033, 0.04175812,
       0.04395592, 0.04615372, 0.04835151, 0.05054931, 0.0527471 ,
       0.0549449 , 0.0571427 , 0.05934049, 0.06153829, 0.06373608,
       0.06593388, 0.06813168, 0.07032947, 0.07252727, 0.07472506,
       0.07692286, 0.07912065, 0.08131845, 0.08351625, 0.08571404,
       0.08791184, 0.09010963, 0.09230743, 0.09450523, 0.09670302,
       0.09890082, 0.10109861, 0.10329641, 0.10549421, 0.107692  ,
       0.1098898 , 0.11208759, 0.11428539, 0.11648319, 0.11868098,
       0.12087878, 0.12307657, 0.12527437, 0.12747217, 0.12966996,
       0.13186776, 0.13406555, 0.13626335, 0.13846115, 0.14065894,
       0.14285674, 0.14505453, 0.14725233, 0.14945013, 0.15164792,
       0.15384572, 0.15604351, 0.15824131, 0.16043911, 0.1626369 ,
       0.1648347 , 0.16703249, 0.16923029, 0.17142809, 0.17362588,
       0.17582368, 0.17802147, 0.18021927, 0.18241707, 0.18461486,
       0.18681266, 0.18901045, 0.19120825, 0.19340605, 0.19560384,
       0.19780164, 0.19999943, 0.20219723, 0.20439503, 0.20659282,
       0.20879062, 0.21098841, 0.21318621, 0.21538401, 0.2175818 ,
       0.2197796 , 0.22197739, 0.22417519, 0.22637299, 0.22857078,
       0.23076858, 0.23296637, 0.23516417, 0.23736196, 0.23955976,
       0.24175756, 0.24395535, 0.24615315, 0.24835094, 0.25054874,
       0.25274654, 0.25494433, 0.25714213, 0.25933992, 0.26153772,
       0.26373552, 0.26593331, 0.26813111, 0.2703289 , 0.2725267 ,
       0.2747245 , 0.27692229, 0.27912009, 0.28131788, 0.28351568,
       0.28571348, 0.28791127, 0.29010907, 0.29230686, 0.29450466,
       0.29670246, 0.29890025, 0.30109805, 0.30329584, 0.30549364,
       0.30769144, 0.30988923, 0.31208703, 0.31428482, 0.31648262,
       0.31868042, 0.32087821, 0.32307601, 0.3252738 , 0.3274716 ,
       0.3296694 , 0.33186719, 0.33406499, 0.33626278, 0.33846058,
       0.34065838, 0.34285617, 0.34505397, 0.34725176, 0.34944956,
       0.35164736, 0.35384515, 0.35604295, 0.35824074, 0.36043854,
       0.36263634, 0.36483413, 0.36703193, 0.36922972, 0.37142752,
       0.37362532, 0.37582311, 0.37802091, 0.3802187 , 0.3824165 ,
       0.3846143 , 0.38681209, 0.38900989, 0.39120768, 0.39340548,
       0.39560327, 0.39780107, 0.39999887, 0.40219666, 0.40439446,
       0.40659225, 0.40879005, 0.41098785, 0.41318564, 0.41538344,
       0.41758123, 0.41977903, 0.42197683, 0.42417462, 0.42637242,
       0.42857021, 0.43076801, 0.43296581, 0.4351636 , 0.4373614};

    const List<double> flowrate_data_Ls{0.00000000e+00, 1.37186674e-02, 6.46436522e-02, 1.55465159e-01,
       2.83887107e-01, 4.45786700e-01, 6.36295259e-01, 8.50327273e-01,
       1.08287196e+00, 1.32916082e+00, 1.58476318e+00, 1.84563704e+00,
       2.10815132e+00, 2.36908919e+00, 2.62563908e+00, 2.87537751e+00,
       3.11624687e+00, 3.34653014e+00, 3.56482411e+00, 3.77001203e+00,
       3.96123651e+00, 4.13787324e+00, 4.29950567e+00, 4.44590118e+00,
       4.57698865e+00, 4.69283754e+00, 4.79363866e+00, 4.87968642e+00,
       4.95136267e+00, 5.00912197e+00, 5.05347834e+00, 5.08499327e+00,
       5.10426506e+00, 5.11191929e+00, 5.10860038e+00, 5.09496422e+00,
       5.07167160e+00, 5.03938264e+00, 4.99875189e+00, 4.95042415e+00,
       4.89503100e+00, 4.83318780e+00, 4.76549131e+00, 4.69251767e+00,
       4.61482092e+00, 4.53293174e+00, 4.44735666e+00, 4.35857746e+00,
       4.26705083e+00, 4.17320833e+00, 4.07745642e+00, 3.98017673e+00,
       3.88172647e+00, 3.78243891e+00, 3.68262398e+00, 3.58256898e+00,
       3.48253931e+00, 3.38277927e+00, 3.28351288e+00, 3.18494477e+00,
       3.08726104e+00, 2.99063017e+00, 2.89520389e+00, 2.80111807e+00,
       2.70849362e+00, 2.61743733e+00, 2.52804275e+00, 2.44039096e+00,
       2.35455145e+00, 2.27058282e+00, 2.18853360e+00, 2.10844291e+00,
       2.03034122e+00, 1.95425094e+00, 1.88018710e+00, 1.80815797e+00,
       1.73816557e+00, 1.67020625e+00, 1.60427124e+00, 1.54034706e+00,
       1.47841602e+00, 1.41845667e+00, 1.36044414e+00, 1.30435057e+00,
       1.25014547e+00, 1.19779599e+00, 1.14726730e+00, 1.09852280e+00,
       1.05152445e+00, 1.00623299e+00, 9.62608141e-01, 9.20608850e-01,
       8.80193461e-01, 8.41319893e-01, 8.03945805e-01, 7.68028733e-01,
       7.33526231e-01, 7.00395985e-01, 6.68595925e-01, 6.38084319e-01,
       6.08819865e-01, 5.80761761e-01, 5.53869780e-01, 5.28104329e-01,
       5.03426497e-01, 4.79798108e-01, 4.57181755e-01, 4.35540832e-01,
       4.14839567e-01, 3.95043036e-01, 3.76117191e-01, 3.58028861e-01,
       3.40745773e-01, 3.24236549e-01, 3.08470713e-01, 2.93418689e-01,
       2.79051799e-01, 2.65342258e-01, 2.52263162e-01, 2.39788486e-01,
       2.27893065e-01, 2.16552589e-01, 2.05743581e-01, 1.95443391e-01,
       1.85630171e-01, 1.76282866e-01, 1.67381191e-01, 1.58905616e-01,
       1.50837346e-01, 1.43158304e-01, 1.35851108e-01, 1.28899058e-01,
       1.22286111e-01, 1.15996864e-01, 1.10016533e-01, 1.04330939e-01,
       9.89264806e-02, 9.37901228e-02, 8.89093734e-02, 8.42722667e-02,
       7.98673449e-02, 7.56836402e-02, 7.17106573e-02, 6.79383565e-02,
       6.43571368e-02, 6.09578196e-02, 5.77316326e-02, 5.46701944e-02,
       5.17654993e-02, 4.90099026e-02, 4.63961063e-02, 4.39171452e-02,
       4.15663733e-02, 3.93374512e-02, 3.72243330e-02, 3.52212542e-02,
       3.33227201e-02, 3.15234940e-02, 2.98185868e-02, 2.82032458e-02,
       2.66729445e-02, 2.52233733e-02, 2.38504293e-02, 2.25502075e-02,
       2.13189920e-02, 2.01532476e-02, 1.90496111e-02, 1.80048846e-02,
       1.70160268e-02, 1.60801465e-02, 1.51944957e-02, 1.43564627e-02,
       1.35635657e-02, 1.28134472e-02, 1.21038676e-02, 1.14327000e-02,
       1.07979247e-02, 1.01976241e-02, 9.62997803e-03, 9.09325865e-03,
       8.58582638e-03, 8.10612544e-03, 7.65267984e-03, 7.22408944e-03,
       6.81902627e-03, 6.43623098e-03, 6.07450944e-03, 5.73272951e-03,
       5.40981799e-03, 5.10475763e-03, 4.81658441e-03, 4.54438477e-03,
       4.28729318e-03, 4.04448963e-03, 3.81519736e-03, 3.59868069e-03,
       3.39424286e-03, 3.20122414e-03, 3.01899983e-03, 2.84697854e-03};

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

    //Foam::Info << "Current flowrate at " << patch().patch().name() << " = "
    //            << outletFlowRate << endl;

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


void Foam::airInletVOFVelocityHealthyFvPatchVectorField::updateCoeffs()
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


void Foam::airInletVOFVelocityHealthyFvPatchVectorField::write(Ostream& os) const
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
       airInletVOFVelocityHealthyFvPatchVectorField
   );
}


// ************************************************************************* //
