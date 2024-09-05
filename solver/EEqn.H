{
    // volScalarField& he = thermo.he();
    // Info<< "min/max(he) before enthalpy NN module= "
    //    << min(he).value() << ", " << max(he).value() << endl;
   // #include "calc_he.H"  // uses NN to get he from (Y1,T,p)
    // Info<< "min/max(he) after enthalpy NN module= "
    //    << min(he).value() << ", " << max(he).value() << endl;

    // #include "calc_rho.H"
    // #include "calc_alpha.H"    


       fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
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
    Info<< "min/max(he) after EEqn is solved = "
        << min(he).value() << ", " << max(he).value() << endl;

    Info<< "min/max(T) before T inference = "
        << min(T).value() << ", " << max(T).value() << endl;
    #include "python__calc_T.H""
    Info<< "min/max(T) after  T inference = "
        << min(T).value() << ", " << max(T).value() << endl;
    

    #include "python__calc_rho.H""
    #include "python__calc_psi.H""
    //#include "python__calc_mu.H""
    //#include "python__calc_alpha.H""
    //rho=thermo.rho();
    //psi=thermo.psi();
    mu=thermo.mu();
    alpha=thermo.alpha();
    //psi = Foam::max(psi, min_psi);
    //mu = Foam::max(mu, min_mu);
    
    //thermo.correct();
}
