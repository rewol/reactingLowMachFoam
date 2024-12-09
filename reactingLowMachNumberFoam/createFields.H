#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

// THIS IS NEW. WE ADD DYNAMIC PRESSURE PERTURBATION pPrime
// This will be solved and drive pressure differences !
volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //dimensionedScalar("pPrime", p.dimensions(), 0.0)
);

// THIS IS NEW. Define reference pressure
/* Define Yourself */
//const scalar pReff = 1.01325e5; // Atmospheric pressure in Pascals
/* Define from Dictionary */
IOdictionary flameProperties
(
	IOobject
	(
	 "flameProperties",
	 runTime.constant(),
	 mesh,
	 IOobject::MUST_READ,
	 IOobject::NO_WRITE
	)
);

// THIS IS NEW. Assign pReff to p. p wont chaange. only pPrime will change
dimensionedScalar pReff(flameProperties.lookup("pReff"));
Info << "pReff is " << pReff << endl;

p = pReff;

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(pPrime.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"
