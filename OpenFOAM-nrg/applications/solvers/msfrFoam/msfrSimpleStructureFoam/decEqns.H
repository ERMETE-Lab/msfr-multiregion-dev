{
    volScalarField Ddec("Ddec", turbulence->nu()/Sc + turbulence->nut()/Sct);

    forAll(dec, i)
    {
        volScalarField& deci = dec[i];

        fvScalarMatrix decEqni
        (
            fvm::div(phi, deci)
          - fvm::laplacian(Ddec, deci)
          + fvm::Sp(decLambda[i], deci)
          // + fvm::Sp(C*(1-liquidFraction)*(1-liquidFraction)/(liquidFraction*liquidFraction*liquidFraction+1e-3), deci)
         ==
            decBeta[i]*EfRf
          + fvOptions(deci)
        );

        decEqni.relax();
        decEqni.solve();
    }
}
