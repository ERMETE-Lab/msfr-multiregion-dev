{
    // Evaluate normalisation factor
    volScalarField normFactor = Q/Qnominal;

    Info<< endl;
    Info<< "Q (MW)       : " << min(Q).value()/1E+06 << endl;
    Info<< "Q nom. (MW)  : " << Qnominal.value()/1E+06 << nl << endl;
    Info<< "Norm. factor : " << min(normFactor).value() << nl << endl;

    // Normalise fluxes
    forAll(flux, i)
    {
        volScalarField& fluxi = flux[i];
        fluxi = fluxi/normFactor;
    }
}
