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
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);

// volScalarField alpha ("alpha", 0.*mu);  // kg/m/s   (kg/m3 * m2/s), already rho*alpha
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.alpha()  // kg/m/s   already rho*alpha (kg/m3 * m2/s)
);
//const int size = 174400;
//i
//int he[size]; // Create an array of size elements

volScalarField he
(
    IOobject
    (
        "he",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, //MUST_READ
        IOobject::AUTO_WRITE
    ),
  //  mesh
    thermo.he()

); 




/*
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);
*/
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
   // mesh
    thermo.psi()
);

volScalarField& T =thermo.T();

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

// --------------------- using NN inference to initialize fields ----------//
#include "python_rho_cp_he_psi_T_alpha_mu__create.H"
//#include "python_rho_cp_he_psi_alpha_mu__get_input_data.H"
//#include "python_T__get_input_data.H"
Info<< "" << nl << endl;
Info<< "overwriting initial fields with NN" << endl;
Info<< "----------------------------------" << endl;

Info<< "min/max(he)    before = " << min(he).value() << ",  " << max(he).value() << endl;
#include "python__calc_he.H"  
Info<< "min/max(he)    after  = " << min(he).value() << ",  " << max(he).value() << endl;

Info<< "min/max(T)     before = " << min(T).value() << ",   " << max(T).value() << endl;
#include "python__calc_T.H"  
Info<< "min/max(T)     after  = " << min(T).value() << ",   " << max(T).value() << endl;

Info<< "min/max(he)    before = " << min(he).value() << ",  " << max(he).value() << endl;
#include "python__calc_he.H"  
Info<< "min/max(he)    after  = " << min(he).value() << ",  " << max(he).value() << endl;

Info<< "min/max(T)     before = " << min(T).value() << ",   " << max(T).value() << endl;
#include "python__calc_T.H"  
Info<< "min/max(T)     after  = " << min(T).value() << ",   " << max(T).value() << endl;

Info<< "min/max(rho)   before = " << min(rho).value() << ", " << max(rho).value() << endl;
#include "python__calc_rho.H"  
Info<< "min/max(rho)   after  = " << min(rho).value() << ", " << max(rho).value() << endl;

Info<< "min/max(psi)   before = " << min(psi).value() << ",  " << max(psi).value() << endl;
#include "python__calc_psi.H"
Info<< "min/max(psi)   after  = " << min(psi).value() << ",  " << max(psi).value() << endl;

Info<< "min/max(mu)    before = " << min(mu).value() << ",  " << max(mu).value() << endl;
//#include "python__calc_mu.H"
Info<< "min/max(mu)    after  = " << min(mu).value() << ",  " << max(mu).value() << endl;

Info<< "min/max(alpha) before = " << min(alpha).value() << ",  " << max(alpha).value() << endl;
//#include "python__calc_alpha.H"
Info<< "min/max(alpha) after  = " << min(alpha).value() << ",  " << max(alpha).value() << endl;

Info<< "min/max(kappa)        = " << min(kappa).value() << ",  " << max(kappa).value() << endl;


/*const dimensionedScalar min_psi
    (
        "min_psi",
        dimensionSet(0, -2, 2, 0, 0, 0, 0),
        scalar(1e-7)
    );

psi = Foam::max(psi, min_psi);

Info<< "min/max(psi) clipped  = " << min(psi).value() << ",  " << max(psi).value() << endl;

const dimensionedScalar min_mu
    (
        "min_mu",
        dimensionSet(1, -1, -1, 0, 0, 0, 0),
        scalar(1e-7)
    );

mu = Foam::max(mu, min_mu);

Info<< "min/max(mu)  clipped  = " << min(mu).value() << ",  " << max(mu).value() << endl;


Info<< "----------------------------------" << endl;
Info<< "done with NN fields initialization" << endl;
Info<< "" << nl << endl;*/

// --------------------- end of NN inference for initializations ----------//

#include "compressibleCreatePhi.H"    // here phi is constructed from rho and U

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

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

