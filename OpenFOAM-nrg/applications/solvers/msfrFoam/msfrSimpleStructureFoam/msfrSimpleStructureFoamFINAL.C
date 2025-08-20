/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    msfrSimpleFoam

Description
    Steady-state solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

// Specialised class for incompressible flows, using kinematic viscosity.
// To be changed if we want to model full compressibility.
#include "kinematicMomentumTransportModel.H" 

// This implements the radiation model interface but without radiation effects.
// To be changed once the radiation module will be fully completed.
#include "noRadiation.H"

#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    // Add to the default fields (T, p_rgh, p, U, turbulence quantities, 
    // and thermo-dynamics properties) the neutronics-specific fields 
    // (total flux, group fluxes, delayed neutron precursors fields, decay heat
    // precursors fields, neutronics properties, volumetric heat sources)
    #include "createNuclearFieldsStructure.H"

    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "normFlux.H" // Evaluate the normalisation factor and normalise the fluxes

    // Here we initialise the temperature and density-dependent nuclear properties (diffusion
    // coefficient and cross sections) with the initial temperature field.
    #include "updateCrossSections.H"

    // Initialise the fission rates based on the initial fluxes and cross-sections
    #include "updateFissionRate.H"

    // Initialise the neutron sources based on the initial fluxes and cross-sections
    #include "updateNeutronSource.H"

    // Initialise the prompt and delayed heat sources from initial fluxes, cross-sections and precursors
    #include "updatePowerSource.H"

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"      // Solve the velocity equation
            #include "TEqn.H"      // Solve the temperature equation

            #include "fluxEqns.H"  // Solve the SP3 multigroup diffusion equations
            #include "precEqns.H"  // Solve the delayed neutron precursor equations
            #include "decEqns.H"   // Solve the decay heat precursors equations

            // Flag to determine whether we consider or not the transport of the 
            // solid fission products (defined in constant/nuclearProperties, default False)
            if (fpTransport)      
            {
                #include "fpEqns.H"  // Solve the fission product transport equation
            }

            #include "pEqn.H"

            // Update the nuclear properties, fission rates, neutron sources and
            // heat sources based on the new values of temperature and density
            #include "updateCrossSections.H"
            #include "updateFissionRate.H"
            #include "updateNeutronSource.H"
            #include "updatePowerSource.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        #include "updateKeff.H"     // For the steady-state solver, compute the effective
                                    // multiplication factor and the reactivity

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
