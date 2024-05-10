{
    forAll(flux, i)
    {
        volScalarField& fluxi = flux[i];

        diffCoeffAlbedo = D[i];

        fluxi.correctBoundaryConditions();

        fvScalarMatrix fluxEqni
        (
          - fvm::laplacian(D[i], fluxi)
          - fvm::Sp(M[i], fluxi)
         ==
            neutronSource[i]
        );

        fluxEqni.relax();
        fluxEqni.solve();
    }

    // Update total flux
    fluxTot *= 0;
    forAll(flux, i)
    {
        fluxTot += flux[i];
    }
}
