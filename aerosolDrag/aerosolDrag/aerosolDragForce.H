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

Class
    Foam::aerosolDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid
    velocity) are dynamically computed based on empirical expressions.

    \f[
        \mathrm{F}_\mathrm{D} =
            \frac{3}{4}
            \frac{\mu_c\,\mathrm{C}_\mathrm{D}\,\mathrm{Re}_p}{\rho_p \, d_p^2}
    \f]
    with

    \f[
        \mathrm{C}_\mathrm{D} =
            \frac{24}{\mathrm{Re}_p}
            \left(1 + \frac{1}{6}\mathrm{Re}_p^{2/3} \right)
            \quad \mathrm{if} \quad \mathrm{Re}_p \leq 1000
    \f]
    \f[
        \mathrm{C}_\mathrm{D} =
            0.424 \quad \mathrm{if} \quad \mathrm{Re}_p > 1000
    \f]
    and
    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    where
    \vartable
        \mathrm{F}_\mathrm{D} | Drag force per carrier-fluid velocity [kg/s]
        \mathrm{C}_\mathrm{D} | Particle drag coefficient
        \mathrm{Re}_p         | Particle Reynolds number
        \rho_p                | Particle mass density
        \mu_c    | Dynamic viscosity of carrier at the cell occupying particle
        d_p                   | Particle diameter
        \rho_c   | Density of carrier at the cell occupying particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
    \endvartable

    Constraints:
    - Particles remain spherical throughout the force
    computation, hence no particle distortion.
    - Applicable to particles with a spatially homogeneous distribution.

    References:
    \verbatim
        Standard model:
            Putnam, A. (1961).
            Integratable form of droplet drag coefficient.
            ARS Journal, 31(10), 1467-1468.

        Expressions (tag:AOB), (Eq. 34-35):
            Amsden, A. A., O'Rourke, P. J., & Butler, T. D. (1989).
            KIVA-II: A computer program for chemically
            reactive flows with sprays (No. LA-11560-MS).
            Los Alamos National Lab.(LANL), Los Alamos, NM (United States).
            DOI:10.2172/6228444
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            aerosolDrag;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: aerosolDrag             | word | yes  | -
    \endtable

Note
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass/density
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

SourceFiles
    aerosolDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef aerosolDragForce_H
#define aerosolDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class aerosolDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class aerosolDragForce
:
    public ParticleForce<CloudType>
{
    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;

        scalar Cc(const typename CloudType::parcelType& p) const;

public:

    //- Runtime type information
    TypeName("aerosolDrag");


    // Constructors

        //- Construct from mesh
        aerosolDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        aerosolDragForce(const aerosolDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new aerosolDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const aerosolDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~aerosolDragForce() = default;


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "aerosolDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
