// Read nuclear properties
    Info<< "\n\nReading nuclear properties" << endl;

    IOdictionary nuclearPropertiesDict
    (
        IOobject
        (
            "nuclearProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

// Read neutron transport constants
    Info<< "\nReading neutron transport constants" << endl;

    dictionary neutronTransportDict(nuclearPropertiesDict.subDict("neutronTransport"));

    label nEnergyGroups(readInt(neutronTransportDict.lookup("energyGroups")));

    // Enable Doppler effect switch for cross sections update
    Switch useDoppler(neutronTransportDict.lookup("useDoppler"));
    dimensionedScalar DopplerControl("DopplerControl", dimless, 1);
    if (!useDoppler)
    {
        DopplerControl = 0;
    }

    // Read reference state
    dimensionedScalar rhoRefXS("rhoRefXS", dimDensity, neutronTransportDict.lookup("rhoRef"));
    dimensionedScalar TRefXS("TRefXS", dimTemperature, neutronTransportDict.lookup("TRef"));

    // For each property, read reference values (and temp. coefficients) from
    // dictionary into a List<scalar>, and declare a PtrList<dimensionedScalar>
    // to be populated
/*    PtrList<dimensionedScalar> D_ref(nEnergyGroups);
    PtrList<dimensionedScalar> alphaD_ref(nEnergyGroups);
//    List<scalar> D_ref_(neutronTransportDict.lookup("DRef"));
    List<scalar> alphaD_ref_(neutronTransportDict.lookup("alphaDRef"));


    PtrList<dimensionedScalar> Sa_ref(nEnergyGroups);
    PtrList<dimensionedScalar> alphaSa_ref(nEnergyGroups);
    List<scalar> Sa_ref_(neutronTransportDict.lookup("SaRef"));
    List<scalar> alphaSa_ref_(neutronTransportDict.lookup("alphaSaRef"));

    PtrList<dimensionedScalar> Sf_ref(nEnergyGroups);
    PtrList<dimensionedScalar> alphaSf_ref(nEnergyGroups);
    List<scalar> Sf_ref_(neutronTransportDict.lookup("SfRef"));
    List<scalar> alphaSf_ref_(neutronTransportDict.lookup("alphaSfRef"));

    PtrList<dimensionedScalar> Ss_ref(nEnergyGroups*nEnergyGroups);
    PtrList<dimensionedScalar> alphaSs_ref(nEnergyGroups*nEnergyGroups);
    List<scalar> Ss_ref_(neutronTransportDict.lookup("SsRef"));
    List<scalar> alphaSs_ref_(neutronTransportDict.lookup("alphaSsRef"));

    PtrList<dimensionedScalar> Nu(nEnergyGroups);
    List<scalar> Nu_(neutronTransportDict.lookup("Nu"));
    PtrList<dimensionedScalar> Ef(nEnergyGroups);
    List<scalar> Ef_(neutronTransportDict.lookup("Ef"));
    PtrList<dimensionedScalar> invVel(nEnergyGroups);
    List<scalar> invVel_(neutronTransportDict.lookup("inverseVelocity"));
    PtrList<dimensionedScalar> chiPrompt(nEnergyGroups);
    List<scalar> chiPrompt_(neutronTransportDict.lookup("promptYeld"));
    PtrList<dimensionedScalar> chiDelayed(nEnergyGroups);
    List<scalar> chiDelayed_(neutronTransportDict.lookup("delayedYeld"));

    // Check input data size
    if( D_ref_.size() == nEnergyGroups )
    {
        Info<< "    Found " << nEnergyGroups << " energy groups" << endl;
    }
    else
    {
        FatalError
            << "in " << neutronTransportDict.name() << " , invalid number of neutron transport constants" << exit(FatalError);
    }

    // Populate pointer lists with dimensioned values
    forAll(Sa_ref_, i)
    {
        word iName = name(i + 1);

  //      D_ref.set(i, new dimensionedScalar("D" + iName + "_ref", dimLength, D_ref_[i]));
        Sa_ref.set(i, new dimensionedScalar("Sa" + iName + "_ref", dimless/dimLength, Sa_ref_[i]));
        Sf_ref.set(i, new dimensionedScalar("Sf" + iName + "_ref", dimless/dimLength, Sf_ref_[i]));

        alphaD_ref.set(i, new dimensionedScalar("alphaD" + iName + "_ref", dimLength, alphaD_ref_[i]));
        alphaSa_ref.set(i, new dimensionedScalar("alphaSa" + iName + "_ref", dimless/dimLength, alphaSa_ref_[i]));
        alphaSf_ref.set(i, new dimensionedScalar("alphaSf" + iName + "_ref", dimless/dimLength, alphaSf_ref_[i]));

        forAll(Sa_ref_, j)
        {
            label ijInd = nEnergyGroups*i + j;
            word ijName = name(i + 1) + name(j + 1);

            Ss_ref.set(ijInd, new dimensionedScalar("Ss" + ijName + "_ref", dimless/dimLength, Ss_ref_[ijInd]));
            alphaSs_ref.set(ijInd, new dimensionedScalar("alphaSs" + ijName + "_ref", dimless/dimLength, alphaSs_ref_[ijInd]));
        }

        Nu.set(i, new dimensionedScalar("Nu" + name(i + 1), dimless, Nu_[i]));
        Ef.set(i, new dimensionedScalar("Ef" + name(i + 1), dimEnergy, Ef_[i]));
        invVel.set(i, new dimensionedScalar("invVel" + name(i + 1), dimTime/dimLength, invVel_[i]));
        chiPrompt.set(i, new dimensionedScalar("chiPrompt" + name(i + 1), dimless, chiPrompt_[i]));
        chiDelayed.set(i, new dimensionedScalar("chiDelayed" + name(i + 1), dimless, chiDelayed_[i]));
    }
// */
// Read delayed neutron precursors constants
   Info<< "\nReading delayed neutron precursors constants" << endl;

    dictionary neutronPrecursorsDict(nuclearPropertiesDict.subDict("delayedNeutronPrecursors"));

    label nPrecGroups(readInt(neutronPrecursorsDict.lookup("groups")));
   // PtrList<dimensionedScalar> precLambda(nPrecGroups);
   // PtrList<dimensionedScalar> precBeta(nPrecGroups);
   // dimensionedScalar precBetaTot("precBetaTot", dimless, 0);

    // Scalar lists to lookup dictionary
   // List<scalar> precLambda_(neutronPrecursorsDict.lookup("lambda"));
   // List<scalar> precBeta_(neutronPrecursorsDict.lookup("beta"));

/*    // Check input data size
    if( precLambda.size() == nPrecGroups && precBeta_.size() == nPrecGroups )
    {
        Info<< "    Found " << nPrecGroups << " groups of delayed neutron precursors" << endl;
    }
    else
    {
        FatalError
            << "in " << neutronPrecursorsDict.name() << " , invalid number of delayed neutron precursors constants" << exit(FatalError);
    }

    // Assign dimensioned values to pointer list
   /* forAll(precLambda_, i)
    {
        precLambda.set(i, new dimensionedScalar("precLambda" + name(i + 1), dimless/dimTime, precLambda_[i]));
        precBeta.set(i, new dimensionedScalar("precBeta" + name(i + 1), dimless, precBeta_[i]));
        precBetaTot += precBeta[i];
    }
*/
// Read decay heat precursors constants
    Info<< "\nReading decay heat precursors constants" << endl;

    dictionary decayHeatPrecursorsDict(nuclearPropertiesDict.subDict("decayHeatPrecursors"));

    label nDecGroups(readInt(decayHeatPrecursorsDict.lookup("groups")));
  //  PtrList<dimensionedScalar> decLambda(nDecGroups);
  //  PtrList<dimensionedScalar> decBeta(nDecGroups);
  //  dimensionedScalar decBetaTot("decBetaTot", dimless, 0);

  //  Scalar lists to lookup dictionary
  //  List<scalar> decLambda_(decayHeatPrecursorsDict.lookup("lambda"));
  //  List<scalar> decBeta_(decayHeatPrecursorsDict.lookup("beta"));

    // Check input data size
/*    if( decLambda.size() == nDecGroups && decBeta_.size() == nDecGroups )
    {
        Info<< "    Found " << nDecGroups << " groups of decay heat precursors" << endl;
    }
    else
    {
        FatalError
            << "in " << decayHeatPrecursorsDict.name() << " , invalid number of decay heat precursors constants" << exit(FatalError);
    }

  /*  // Assign dimensioned values to pointer list
    forAll(decLambda_, i)
    {
        decLambda.set(i, new dimensionedScalar("decLambda" + name(i + 1), dimless/dimTime, decLambda_[i]));
        decBeta.set(i, new dimensionedScalar("decBeta" + name(i + 1), dimless, decBeta_[i]));
        decBetaTot += decBeta[i];
    }
*/
// Read fission products constants (if present, otherwise leave empty)
    bool fpTransport(nuclearPropertiesDict.found("fissionProducts"));
    label nFPSpecies(0);
    PtrList<dimensionedScalar> FPYield;
    PtrList<dimensionedScalar> FPLambda;

    if (fpTransport)
    {
        Info<< "\nReading fission products constants" << endl;

        dictionary fissionProductsDict(nuclearPropertiesDict.subDict("fissionProducts"));

        nFPSpecies = readInt(fissionProductsDict.lookup("species"));
        FPYield.setSize(nFPSpecies);
        FPLambda.setSize(nFPSpecies);

        // Scalar lists to lookup dictionary
        List<scalar> FPYield_(fissionProductsDict.lookup("yield"));
        List<scalar> FPLambda_(fissionProductsDict.lookup("lambda"));

        // Check input data size
        if( FPLambda_.size() == nFPSpecies && FPYield_.size() == nFPSpecies )
        {
            Info<< "    Found " << nFPSpecies << " species of fission products" << endl;
        }
        else
        {
            FatalError
                << "in " << fissionProductsDict.name() << " , invalid number of fission product species constants" << exit(FatalError);
        }

        // Assign dimensioned values to pointer list
        forAll(FPYield_, i)
        {
            FPYield.set(i, new dimensionedScalar("FPYield" + name(i + 1), dimless, FPYield_[i]));
            FPLambda.set(i, new dimensionedScalar("FPLambda" + name(i + 1), dimless/dimTime, FPLambda_[i]));
        }
    }
