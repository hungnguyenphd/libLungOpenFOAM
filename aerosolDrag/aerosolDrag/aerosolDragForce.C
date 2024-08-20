/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "aerosolDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::aerosolDragForce<CloudType>::CdRe(const scalar Re) const
{
    return 24.0;
}

template<class CloudType>
Foam::scalar Foam::aerosolDragForce<CloudType>::Cc
(
    const typename CloudType::parcelType& p
) const
{
    //Cunningham correction factor
    //(wiki: https://en.wikipedia.org/wiki/Cunningham_correction_factor)

    const scalar lambda = 6.14e-08; //mean free path (Jimin's draft)

    const scalar A1 = 1.257;

    const scalar A2 = 0.400;

    const scalar A3 = 0.55;

    return 1 + (2*lambda/p.d())*(A1 + A2*exp(-A3*(p.d()/(lambda))));
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::aerosolDragForce<CloudType>::aerosolDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::aerosolDragForce<CloudType>::aerosolDragForce
(
    const aerosolDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::aerosolDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar alpha = 1; //particle-particle interaction factor
                            // set to 1 due to the negligible particle
                            // volume per unit gas volume of 1.30*10^(-7)
                            // for 2.5um particles (Miyawaki et al., 2012,
                            // Annals of Biomedical Engineering)

    return forceSuSp(Zero, mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d())*Cc(p)*pow(alpha,3.17)));
}


// ************************************************************************* //
