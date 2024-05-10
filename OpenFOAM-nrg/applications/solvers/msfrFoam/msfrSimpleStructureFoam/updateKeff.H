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

    Info<< "Keff old     : " << min(Keff).value() << endl;

    // Update multiplication factor
    Keff = Keff*normFactor;

    // Bound multiplication factor to limit power excursions between iterations
    dimensionedScalar minKeff ("minKeff", dimless, 0.75);
    dimensionedScalar maxKeff ("maxKeff", dimless, 1.25);

    if (min(Keff).value() > maxKeff.value())
    {
        Keff = maxKeff*unitField;
    }
    else if (min(Keff).value() < minKeff.value())
    {
        Keff = minKeff*unitField;
    }

    // Evaluate reactivity
    reactivity = (Keff - 1)/Keff;

    Info<< "Keff         : " << min(Keff).value() << nl << endl;
    Info<< "React. (pcm) : " << min(reactivity).value()*1e5 << nl << endl;
}
