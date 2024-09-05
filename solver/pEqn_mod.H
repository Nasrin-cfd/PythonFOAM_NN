//rho = thermo.rho();
Info<< "min/max(rho) before 1st pEqn inference = "
    << min(rho).value() << ", " << max(rho).value() << endl;
#include "python__calc_rho.H"
Info<< "min/max(rho) after  1st pEqn inference = "
    << min(rho).value() << ", " << max(rho).value() << endl;

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

// MB
const volScalarField p0(p);     // check
const volScalarField psi0(psi); //
const volScalarField rho0(rho); //
//

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi0, p)     // MB check
          + fvc::ddt(rho0)        // MB check
          - fvc::ddt(psi0,p0)     // MB check
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
        Info<< "I am in pimple.transonic - pEqn_mod.H" << endl;
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
        )
    );

     MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi0, p)     // MB check
          + fvc::ddt(rho0)        // MB check
          - fvc::ddt(psi0,p0)     // MB check             
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
        Info<< "I am in pimple std - pEqn_mod.H" << endl;
    }
}

#include "rhoEqn.H"
//#include "NNcompressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

Info<< "min/max(p) after solving pEqn_mod = "
    << min(p).value() << ", " << max(p).value() << endl;

//rho = thermo.rho();
Info<< "min/max(rho) before 2nd pEqn inference = "
    << min(rho).value() << ", " << max(rho).value() << endl;
#include "python__calc_rho.H"
Info<< "min/max(rho) after  2nd pEqn inference = "
    << min(rho).value() << ", " << max(rho).value() << endl;

/*if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
} */

dpdt = fvc::ddt(p);
