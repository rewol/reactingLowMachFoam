
rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

// Subsonic (Low Mach) Flow Handling
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi))
    )
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Add pressure perturbation gradient terms
phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(pPrime)*mesh.magSf();
HbyA -= (rAU - rAtU)*fvc::grad(pPrime); // Use pPrime instead of p for low Mach number

volScalarField rhorAtU("rhorAtU", rho*rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(pPrime, rho, U, phiHbyA, rhorAtU, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pPrimeEqn
    (
        fvc::ddt(psi) * pReff // Quasi-incompressible term
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAtU, pPrime)
     ==
        fvOptions(psi, pPrime, rho.name())
    );

    pPrimeEqn.solve(mesh.solver(pPrime.select(pimple.finalInnerIter())));
    pPrime.correctBoundaryConditions();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pPrimeEqn.flux(); // Update flux with solved pPrime
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Relax pPrime for momentum corrector
pPrime.relax();

// Update velocity field based on solved pPrime
U = HbyA - rAtU*fvc::grad(pPrime);
U.correctBoundaryConditions();

// Apply fvOptions corrections
fvOptions.correct(U);
K = 0.5*magSqr(U);

// Ensure pressure consistency and update density if needed
if (pressureControl.limit(pPrime))
{
    //p.internalField() = pReff + pPrime.internalField();
    pPrime.correctBoundaryConditions();
    rho = thermo.rho();
}

// Time derivative of pPrime for energy equation
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

