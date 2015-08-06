{
  //alphat = turbulence->nut()/Prt;
  //alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
      - fvc::div(qwall)
    );

    TEqn.relax();
    TEqn.solve();

    rhok = 1.0 - ( (T - TRef)/TRef );
}
