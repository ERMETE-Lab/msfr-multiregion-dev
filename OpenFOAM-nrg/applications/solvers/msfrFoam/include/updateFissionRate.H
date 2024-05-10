{
    // Reset fields
    Rf *= 0;
    NuRf *= 0;
    EfRf *= 0;

    // Evaluate new fission rates based on current fluxes and cross-sections
    forAll(flux, i)
    {
        Rf += Sf[i]*flux[i];
        NuRf += Nu[i]*Sf[i]*flux[i];
        EfRf += Ef[i]*Sf[i]*flux[i];
    }
}
