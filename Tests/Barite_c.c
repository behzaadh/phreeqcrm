#if defined(USE_MPI)
#include <mpi.h>
#endif
//#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RM_interface_C.h"
void Barite_c()
{
    int nCells = 10;
    int nthreads = 1;
    int id = RM_Create(nCells, nthreads);
    int status  = RM_UseSolutionDensityVolume(id, 1);
    status = RM_LoadDatabase(id, "phreeqc.dat");
    status = RM_RunFile(id, 1, 1, 1, "barite.pqi");
    status = RM_SetSelectedOutputOn(id, 1);

        // Set array of initial conditions
    int* ic1 = (int*)malloc((size_t)(7 * nCells * sizeof(int)));
    int* ic2 = (int*)malloc((size_t)(7 * nCells * sizeof(int)));
    double* f1 = (double*)malloc((size_t)(7 * nCells * sizeof(double)));
    if (ic1 == NULL || ic2 == NULL || f1 == NULL) exit(4);
    for (int i = 0; i < nCells; i++)
    {
        ic1[i] = 1;       // Solution 1
        // if (i > 5) {
            // ic1[nCells + i] = 2;
        // } else {
            ic1[nCells + i] = i+1;      // Equilibrium phases 1
        // }
        ic1[2 * nCells + i] = -1;       // Exchange none
        ic1[3 * nCells + i] = -1;      // Surface none
        ic1[4 * nCells + i] = -1;      // Gas phase none
        ic1[5 * nCells + i] = -1;      // Solid solutions none
        ic1[6 * nCells + i] = -1;      // Kinetics none
        ic2[i] = -1;      // Solution none
        ic2[nCells + i] = -1;      // Equilibrium phases none
        ic2[2 * nCells + i] = -1;      // Exchange none
        ic2[3 * nCells + i] = -1;      // Surface none
        ic2[4 * nCells + i] = -1;      // Gas phase none
        ic2[5 * nCells + i] = -1;      // Solid solutions none
        ic2[6 * nCells + i] = -1;      // Kinetics none
        f1[i] = 1.0;      // Mixing fraction ic1 Solution
        f1[nCells + i] = 1.0;      // Mixing fraction ic1 Equilibrium phases
        f1[2 * nCells + i] = 1.0;      // Mixing fraction ic1 Exchange 1
        f1[3 * nCells + i] = 1.0;      // Mixing fraction ic1 Surface
        f1[4 * nCells + i] = 1.0;      // Mixing fraction ic1 Gas phase
        f1[5 * nCells + i] = 1.0;      // Mixing fraction ic1 Solid solutions
        f1[6 * nCells + i] = 1.0;      // Mixing fraction ic1 Kinetics
    }

    int nPPAssemblage = RM_GetPPAssemblageCount(id);
    fprintf(stderr, "Assemblage components:\n");
    for (int i = 0; i < nPPAssemblage; ++i) {
        char str[100];
        status = RM_GetPPAssemblageComp(id, i, str, 100);
        fprintf(stderr, "\t%d: %s\n", i, str);
    }

    double assemblageMole[] = {1,2,3,4,5,6,7,8,9,10,
                               1,2,3,4,5,6,7,8,9,10,
                               10,2,3,4,5,6,7,8,9,10,
                               1,2,3,4,5,6,7,8,9,10,
                               1,2,3,4,5,6,7,8,9,10};
    status = RM_SetPPAssemblageMoles(id, assemblageMole);
    double assemblageSI[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
                             0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
                             1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
                             0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
                             0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
    status = RM_SetPPAssemblageSI(id, assemblageSI);

    status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1);

    // int nPPAssemblage = RM_GetPPAssemblageCount(id);
    // fprintf(stderr, "Assemblage components:\n");
    // for (int i = 0; i < nPPAssemblage; ++i) {
    //     char str[100];
    //     status = RM_GetPPAssemblageComp(id, i, str, 100);
    //     fprintf(stderr, "\t%d: %s\n", i, str);
    // }

    int ncomps = RM_FindComponents(id);
    char str[100];
    char str1[200];

    // Get component information
    char** components = (char**)malloc((size_t)(ncomps * sizeof(char*)));
    if (components == NULL) exit(4);
    for (int i = 0; i < ncomps; i++)
    {
        components[i] = (char*)malloc((size_t)(100 * sizeof(char*)));
        if (components[i] == NULL) exit(4);
        status = RM_GetComponent(id, i, components[i], 100);
        snprintf(str, sizeof(str), "%10s\n", components[i]);
        status = RM_OutputMessage(id, str);
    }

    int ngas = RM_GetGasComponentsCount(id);

    // Get gas component names
    char** gas_comps = (char**)malloc((size_t)(ngas * sizeof(char*)));
    for (int i = 0; i < ngas; i++)
    {
        gas_comps[i] = (char*)malloc((size_t)(100 * sizeof(char*)));
        status = RM_GetGasComponentsName(id, i, gas_comps[i], 100);
    }

    status = RM_SetSpeciesSaveOn(id, 1);

    // Set initial porosity
    double* por = (double*)malloc((size_t)(nCells * sizeof(double)));
    for (int i = 0; i < nCells; i++) por[i] = 0.3;
    status = RM_SetPorosity(id, por);

    // Set initial saturation
    double* sat = (double*)malloc((size_t)(nCells * sizeof(double)));
    for (int i = 0; i < nCells; i++) sat[i] = 1.0;
    status = RM_SetSaturationUser(id, sat);

    // Set representative volume
    double* rv = (double*)malloc((size_t)(nCells * sizeof(double)));
    for (int i = 0; i < nCells; i++) rv[i] = 1.0;
    status = RM_SetRepresentativeVolume(id, rv);

    status = RM_SetUnitsSolution(id, 2); // 1, mg/L; 2, mol/L; 3, kg/kgs

    // double* C_components;
    // C_components = (double*)malloc((size_t)(ncomps * nCells * sizeof(double)));

    double C_components[] = {
        55.5253251007785, 55.5253251007785, 55.5253251007785, 55.5253251007785,
        55.5253251007785, 55.5253251007785, 55.5253251007785, 55.5253251007785,
        55.5253251007785, 55.5253251007785, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.999999999999998, 0.999999999999962, 0.999999999999188, 0.999999999981096,
        0.999999999522609, 0.999999986852127, 0.999999602677431, 0.999986727917285,
        0.999505492364312, 0.979214498573520, 0.0100000000000002, 0.0100000000000034,
        0.0100000000000731, 0.0100000000017014, 0.0100000000429656, 0.0100000011833204,
        0.0100000357593888, 0.0100011944993894, 0.0100445061322733, 0.0118707138355216,
        0.0999999999999998, 0.0999999999999966, 0.0999999999999269, 0.0999999999982986,
        0.0999999999570344, 0.0999999988166796, 0.0999999642406112, 0.0999988055006106,
        0.0999554938677268, 0.0981292861644784, 0.0999999999999998, 0.0999999999999966,
        0.0999999999999269, 0.0999999999982986, 0.0999999999570344, 0.0999999988166796,
        0.0999999642406112, 0.0999988055006106, 0.0999554938677268, 0.0981292861644784,
        0.0999999999999998, 0.0999999999999966, 0.0999999999999269, 0.0999999999982986,
        0.0999999999570344, 0.0999999988166796, 0.0999999642406112, 0.0999988055006106,
        0.0999554938677268, 0.0981292861644784, 0.00100000000000019, 0.00100000000000373,
        0.00100000000008041, 0.00100000000187150, 0.00100000004726217, 0.00100000130165246,
        0.00100003933532771, 0.00100131394932830, 0.00104895674550058, 0.00305778521907373
    };

    status = RM_SetConcentrations(id, C_components);

    double T[] = {76.6666666666667,
                  76.6666666666667,
                  76.6666666666667,
                  76.6666666666667,
                  76.6666666666667,
                  76.6666666666667,
                  76.6666666666667,
                  76.6665988498264,
                  76.6608513726128,
                  76.1097971598307};
    status = RM_SetTemperature(id, T);
    status = RM_RunCells(id);

    int n_user = RM_GetNthSelectedOutputUserNumber(id, 0);
    int col = RM_GetSelectedOutputColumnCount(id);

    // Allocate(selected_out(nxyz,col))
    double* selected_out = (double*)malloc((size_t)(col * nCells * sizeof(double)));
    status = RM_GetSelectedOutput(id, selected_out);
    // Print results

    char heading[100];
    for (int i = 0; i < 1; i++)
    {
        fprintf(stderr, "Cell number %d\n", i);
        fprintf(stderr, "     Components: \n");
        fprintf(stderr, "     Selected output: \n");
        for (int j = 0; j < col; j++)
        {
            status = RM_GetSelectedOutputHeading(id, j, heading, 100);
            fprintf(stderr, "\t%2d %10s: %10.4f\n", j, heading, selected_out[j * nCells + i]);
        }
    }
    double* c_out_comps = (double*)malloc((size_t)(ncomps * nCells * sizeof(double)));
    status = RM_GetConcentrations(id, c_out_comps);    

    fprintf(stderr, "---------- Components ----------\n");

    int k=0;
    for (int j = 0; j < ncomps; j++)
    {
        status = RM_GetComponent(id, j, str, 100);
        fprintf(stderr, "%s\n", str);
        for (int i = 0; i < nCells; i++)
        {
            fprintf(stderr, "     Cell %d: %10.4f\n", i, c_out_comps[k]);
            k++;
        }
    }

    int nspecies = RM_GetSpeciesCount(id);

    double* species_c= (double*)malloc((size_t)(nspecies * nCells * sizeof(double)));
    status = RM_GetSpeciesConcentrations(id, species_c);

    fprintf(stderr, "---------- Species ----------\n");

    k=0;
    for (int j = 0; j < nspecies; j++)
    {
        status = RM_GetSpeciesName(id, j, str, 100);
        fprintf(stderr, "%s\n", str);
        for (int i = 0; i < nCells; i++)
        {
            fprintf(stderr, "     Cell %d: %10.4f\n", i, species_c[k]);
            k++;
        }
    }

    status = RM_Destroy(id);
}
