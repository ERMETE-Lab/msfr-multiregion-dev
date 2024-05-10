{
    volScalarField Dprec("Dprec", turbulence->nu()/Sc + turbulence->nut()/Sct);

    forAll(prec, i)
    {
        volScalarField& preci = prec[i];

        fvScalarMatrix precEqni
        (
            fvm::div(phi, preci)
          - fvm::laplacian(Dprec, preci)
          + fvm::Sp(precLambda[i], preci)
          // + fvm::Sp( C*(1-liquidFraction)*(1-liquidFraction)/(liquidFraction*liquidFraction*liquidFraction+1e-3), preci)
         ==
            precBeta[i]*NuRf
          + fvOptions(preci)
        );

        precEqni.relax();
        precEqni.solve();
    }
}
