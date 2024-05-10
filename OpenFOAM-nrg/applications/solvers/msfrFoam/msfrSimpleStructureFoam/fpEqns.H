{
    volScalarField Dfp("Dfp", turbulence->nu()/Sc + turbulence->nut()/Sct);

    forAll(fp, i)
    {
        volScalarField& fpi = fp[i];

        fvScalarMatrix fpEqni
        (
            fvm::div(phi, fpi)
          - fvm::laplacian(Dfp, fpi)
          + fvm::Sp(FPLambda[i], fpi)
         ==
            FPYield[i]*Rf
          + fvOptions(fpi)
        );

        fpEqni.relax();
        fpEqni.solve();

        dimensionedScalar fpiTot(fvc::domainIntegrate(fpi));

        volScalarField& fpdi = fpd[i];
        volScalarField lambdai = unitField*FPLambda[i];

        volScalarField::Boundary& fpdiBf = fpdi.boundaryFieldRef();

        const volScalarField::Boundary& fpiBf = fpi.boundaryField();
        const volScalarField::Boundary& DfpBf = Dfp.boundaryField();
        const volScalarField::Boundary& lambdaiBf = lambdai.boundaryField();

        dimensionedScalar fpdiTot("fpdiTot", dimless, 0.0);

        forAll(fpdiBf, patchi)
        {
            if (!fpdiBf[patchi].coupled())
            {
                fpdiBf[patchi] = - DfpBf[patchi]*fpiBf[patchi].snGrad()/lambdaiBf[patchi];
                fpdiTot += gSum(fpdiBf[patchi]*fpdiBf[patchi].patch().magSf());
            }
        }

        Info<< "Total specie: " << fpiTot.value() << " atoms (fluid), " << fpdiTot.value() << " atoms (walls)" << endl;
        Info<< "Total balance: " << fvc::domainIntegrate(FPYield[i]*Rf).value() << " atoms/s (source), " << (fpiTot.value() + fpdiTot.value())*FPLambda[i].value() << " atoms/s (decay)" << endl;
    }
}
