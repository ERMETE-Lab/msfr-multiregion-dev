    PtrList<volScalarField> D_ref(nEnergyGroups);
    PtrList<volScalarField> alphaD_ref(nEnergyGroups);
    PtrList<volScalarField> Sa_ref(nEnergyGroups);
    PtrList<volScalarField> alphaSa_ref(nEnergyGroups);
    PtrList<volScalarField> Sf_ref(nEnergyGroups);
    PtrList<volScalarField> alphaSf_ref(nEnergyGroups);
    PtrList<volScalarField> Ss_ref(nEnergyGroups*nEnergyGroups);
    PtrList<volScalarField> alphaSs_ref(nEnergyGroups*nEnergyGroups);
    PtrList<volScalarField> Nu(nEnergyGroups);
    PtrList<volScalarField> Ef(nEnergyGroups);
    PtrList<volScalarField> invVel(nEnergyGroups);
    PtrList<volScalarField> chiPrompt(nEnergyGroups);
    PtrList<volScalarField> chiDelayed(nEnergyGroups);


    // Check input data size
    if( D_ref.size() == nEnergyGroups )
    {
        Info<< "    Found " << nEnergyGroups << " energy groups" << endl;
    }
    else
    {
        FatalError
    << "in " << neutronTransportDict.name() << " , invalid number of neutron transport constants" << exit(FatalError);
    }



   forAll(D_ref, i)
    {
        word iName = name(i + 1);

        D_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "D" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            //   dimensionedScalar("", dimLength, 0)//   D_ref[i]
                )
        );

        Sa_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sa" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
               // dimensionedScalar("", dimless/dimLength, 0)//Sa_ref[i]
            )
        );

       Sf_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sf" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
               // dimensionedScalar("", dimless/dimLength, 0)//Sf_ref[i]
            )
        );


        alphaD_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "alphaD" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
          //      dimensionedScalar("", dimLength, 0)//   D_ref[i]
            )
        );


       alphaSa_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "alphaSa" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
               // dimensionedScalar("", dimless/dimLength, 0)//Sa_ref[i]
            )
        );

        alphaSf_ref.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "alphaSf" + iName + "_ref",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            //    dimensionedScalar("", dimless/dimLength, 0)//Sf_ref[i]
            )
        );

      forAll(D_ref, j)
        {
            label ijInd = nEnergyGroups*i + j;
            word ijName = name(i + 1) + name(j + 1);

            Ss_ref.set
            (
                ijInd,
                new volScalarField
                (
                    IOobject
                    (
                        "Ss" + ijName + "_ref",
                        runTime.timeName(),
                        mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                    ),
                    mesh
                   //  dimensionedScalar("", dimless/dimLength, 0)//Ss_ref[ijInd]
                )
            );
      alphaSs_ref.set
            (
                ijInd,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaSs" + ijName + "_ref",
                        runTime.timeName(),
                        mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                    ),
                    mesh
              //       dimensionedScalar("", dimless/dimLength, 0)//Ss_ref[ijInd]
                )
            );
        }


       Nu.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Nu" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
                //dimensionedScalar("", dimless, 0)//Sf_ref[i]
            )
        );

        Ef.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Ef" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
                //dimensionedScalar("", dimEnergy, 0)//Sf_ref[i]
            )
        );

        invVel.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "invVel" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
               // dimensionedScalar("", dimTime/dimLength, 0)  //Sf_ref[i]
            )
        );

     chiPrompt.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "chiPrompt" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
                //dimensionedScalar("", dimless, 0)   //Sf_ref[i]
            )
        );

       chiDelayed.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "chiDelayed" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
     //           dimensionedScalar("", dimless, 0) //Sf_ref[i]
            )
        );

    }
