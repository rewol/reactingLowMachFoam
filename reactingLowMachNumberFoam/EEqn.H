/*
{
volScalarField& he = thermo.he(); // Access enthalpy or internal energy field

    // Construct the energy equation
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)                // Time derivative of rho * he
      + mvConvection->fvmDiv(phi, he)    // Convection term
      //+ fvc::ddt(rho, K)                 // Kinetic energy contribution
      //+ fvc::div(phi, K)                 // Divergence of kinetic energy flux
      - fvm::laplacian(turbulence->alphaEff(), he) // Diffusion term
     ==
        Qdot                              // Heat source term
      + fvOptions(rho, he)                // User-defined energy sources
    );

    // Remove compressibility effects (low-Mach)
    // Eliminate dp/dt term as thermodynamic pressure (p0) is constant
    if (he.name() == "e")
    {
        EEqn -= fvc::div
        (
            fvc::absolute(phi / fvc::interpolate(rho), U),
            p,
            "div(phiv,p)"
        );
    }

    // Relaxation for numerical stability
    EEqn.relax();

    // Apply additional constraints (e.g., boundary conditions)
    fvOptions.constrain(EEqn);

    // Solve the energy equation
    EEqn.solve();

    // Correct energy-related fields after solving
    fvOptions.correct(he);

    // Update thermodynamic properties
    thermo.correct();

    // Debug output for monitoring temperature range
    Info << "min/max(T) = "
         << min(thermo.T()).value() << ", " << max(thermo.T()).value() << endl;
}
*/


{
    
    volScalarField& he = thermo.he();	
	
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
//      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : 0 * -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}

