
rho = thermo.rho();

// Reciprocal of momentum equation diagonal coefficient
volScalarField rAU(1.0/UEqn.A());

surfaceScalarField rhorAUf("rhorAuf", fvc::interpolate(rho*rAU));

// Compute intermediate velocity field (HbyA)
//volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
volVectorField HbyA("HbyA",U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        (fvc::interpolate(rho*HbyA) & mesh.Sf())
     // + rhorAUf*fvc::ddtCorr(rho, U, phi)
     + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

//fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

// Handle subsonic (Low Mach) flow regime
/*
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA) // Flux based on predicted velocity and density
     // + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi))
    )
);
*/

//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Add density-weighted pressure gradient terms
//phiHbyA += fvc::interpolate(rho*rAU)*fvc::snGrad(pPrime)*mesh.magSf();

// Define a reference thermodynamic pressure
// const scalar pReff = 1.01325e5; // Atmospheric pressure (example)

// Solve for perturbation pressure (pPrime)
while (pimple.correctNonOrthogonal())
{
    /*
    fvScalarMatrix pPrimeEqn
    (
        fvm::laplacian(rhorAUf, pPrime)
        == fvc::div(phiHbyA)
        + fvc::ddt(psi) * pReff // Include the time derivative of reference pressure
    );
   */	
	
    
    fvScalarMatrix pPrimeEqn
    (
        fvc::ddt(psi) * pReff // Quasi-incompressible term
      + fvc::div(phiHbyA) // Flux divergence
     // fvm::laplacian(rhorAUf, pPrime) // Density-weighted Laplacian
      -fvm::laplacian(rho*rAU, pPrime)
     ==
        fvOptions(psi, pPrime, rho.name())
    );
   
    fvOptions.constrain(pPrimeEqn); 

    pPrimeEqn.solve(mesh.solver(pPrime.select(pimple.finalInnerIter())));
    //pPrime.boundaryField().updateCoeffs(pPrimeEqn.snGrad());
    //pPrime.correctBoundaryConditions();
    
    // Print out pPrime max and min values for convergence check
    Info << "min/max(pPrime) = "
    << min(pPrime()).value() << ", " << max(pPrime()).value() << endl;

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pPrimeEqn.flux(); // Update flux with solved pPrime
    }
}


#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Update the total pressure field
//p.internalField() = pReff + pPrime.internalField();
//p.correctBoundaryConditions();
//p = pReff + pPrime;


// Update velocity field based on solved pPrime
U = HbyA - rAU*fvc::grad(pPrime);
U.correctBoundaryConditions();

// Apply fvOptions corrections
fvOptions.correct(U);

// Calculate kinetic energy for post-processing
K = 0.5*magSqr(U);

// Ensure pressure consistency and update density if needed
/*
if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
    rho = thermo.rho();
}
*/
// Compute time derivative of pPrime for energy or other equations
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

