/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "odePyjac.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::odePyjac<ChemistryModel>::odePyjac(typename ChemistryModel::reactionThermo& thermo)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("odeCoeffs")),
    odeSolver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())  // cTp is initialized with (Nsp - 1) + temp + pressure
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::odePyjac<ChemistryModel>::~odePyjac()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::odePyjac<ChemistryModel>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cTp_);
    }

    const label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector (N+1)
    // First element of solver vector will be pressure
    cTp_[0] = p;
    // Second element of solver vector will be temperature
    cTp_[1] = T;
    // Now we update the rest of the vector for Nsp - 1 species:

    for (int i=0; i<nSpecie - 1; i++)
    {
        cTp_[i+2] = c[i];
    }

    // Here we call the ODE solver. This could be the solver we select in chemistryProperties
    odeSolver_->solve(0, deltaT, cTp_, subDeltaT);
    
    // Now we back substitute, now that cTp_ is edited based on ODE solver
    p = cTp_[0];
    T = cTp_[1];
    scalar csum = 0;
    // Also for Yi
    for (int i=0; i<nSpecie - 1; i++)
    {
        c[i] = max(0.0, cTp_[i+2]);
	csum += c[i];
    }
    // Update last species as well
    c[nSpecie - 1] = 1 - csum;
}


// ************************************************************************* //
