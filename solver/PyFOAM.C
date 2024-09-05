/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    NNICELIB-OpenFOAM_Unipg_keras_MB

Group
    grpCombustionSolvers

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
//#include "nnice_single.h"
//...................................................pyhton.....................................................
// Following code to link with Python
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <omp.h>
// #include </home/michele/.local/lib/python3.10/site-packages/numpy/core/include/numpy/arrayobject.h>
// # include </usr/local/lib/python3.10/dist-packages/numpy/core/include/numpy/arrayobject.h>
/*Done importing Python functionality*/


// This initializes numpy for use in OpenFOAM
void init_numpy() {
  import_array1();
} 


/*void distributeData(const std::vector<double>& pValueHe, std::vector<double>& he, int rank, int num_ranks) {
//void distributeData(double pValueHe, double he, int rank, int num_ranks) {
                 //  int num_ranks=3;
                   int size=PyArray_DIM(pValueHe, 0);
                   int num_elements_per_rank = size / num_ranks;
                   int start_index = rank * num_elements_per_rank;
                   int end_index = (rank + 1) * num_elements_per_rank;
                   if (rank == num_ranks - 1) {
                   end_index = total_elements; // Last rank may have fewer elements
                     }

    // Copy data from pValueHe to he for the assigned range
                   for (int i = start_index; i < end_index; ++i) {
                      he[i - start_index] = pValueHe[i];
                          }
                    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <fstream>
#include <vector>
#include <stdexcept>
#include <iomanip>
//void writeScalarField(const std::vector<double>& values, const std::string& fileName);

int main(int argc, char *argv[])
{
 

    (
        "Solver for combustion with chemical reactions"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

/////////////////////////////////////////////////////////////////////////    





///////////////////////////////////////////////////////////////////////    
//    double* he1; // Assuming he is a std::vector
//    double* pValueHe; // Assuming pValueHe is a pointer to the NumPy array
   /* #include "python_rho_cp_he_psi_alpha_mu__get_input_data.H"
     pValueHe = reinterpret_cast<PyArrayObject*>(
                // PyObject_CallObject(pFunc_rho, pArgs_rho)
                PyObject_CallObject(pFunc_he, my_func_args)
                );*/
    
    
    
    
     /*forAll(U.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                         //std::array<int,1> mu;

                         U[id]= *((double*)PyArray_GETPTR1(pValueHe, id));

                 }*/
////////////////////////////////////////////////////////////////////////////////

 
 // Assuming pValueHe is your array of double values and id is the index

     
     
   turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        while (pimple.loop())
        {
            #include "UEqn.H"
            // #include "YEqn.H"  // original location
            // #include "EEqn.H"  // original location

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "YEqn.H"    // alternative location 
                #include "EEqn.H"    // alternative location
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    // #include "pEqn.H"
                    #include "pEqn_mod.H"  // modified pressure equation
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
	        Info<< "min/max(p) = "
            << min(p).value() << ", " << max(p).value() << endl;
            Info<< "min/max(rho) = "
            << min(rho).value() << ", " << max(rho).value() << endl;
        }

        // rho = thermo.rho();
//        #include "python__calc_rho.H"

        // Additional output fields
        // he = thermo.he();
        // Cp = thermo.Cp();
        // mu = thermo.mu();
        // kappa = thermo.kappa();

        #include "Infos.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

/*void writeScalarField(const std::vector<double>& values, const std::string& fileName) {
    std::ofstream outFile(fileName);
    if (!outFile.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    outFile << "| =========                 |\n";
    outFile << "| \\\\      /  F ield         |\n";
    outFile << "|  \\\\    /   O peration     |\n";
    outFile << "|   \\\\  /    A nd           |\n";
    outFile << "|    \\\\/     M anipulation  |\n";
    outFile << "\\*---------------------------------------------------------------------------*///\n";
   /* outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    version     2.0;\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       volScalarField;\n";
    outFile << "    location    \"0\";\n";
    outFile << "    object      T;\n"; // change to appropriate field name
    outFile << "}\n";
    outFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    outFile << "dimensions      [0 2 -2 0 0 0 0];\n"; // adjust dimensions if necessary
    outFile << "internalField   nonuniform List<scalar> \n";
    outFile << values.size() << "\n";
    outFile << "(\n";
    
    outFile << std::setprecision(10);
    for (const auto& value : values) {
        outFile << value << "\n";
    }
    outFile << ");\n\n";
    
    outFile << "boundaryField\n";
    outFile << "{\n";
    // Define boundary conditions appropriately
    outFile << "    inlet\n";
    outFile << "    {\n";
    outFile << "        type fixedValue;\n";
    outFile << "        value uniform 0;\n"; // adjust as necessary
    outFile << "    }\n";
    outFile << "    outlet\n";
    outFile << "    {\n";
    outFile << "        type zeroGradient;\n";
    outFile << "    }\n";
    // Add other boundary conditions as needed
    outFile << "}\n";
    outFile << "\n// ************************************************************************* //\n";

    outFile.close();
}*/
// ************************************************************************* //
