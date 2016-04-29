#include <iostream>
#include "EigenSetup.h"
#include <fstream>
#include <cmath>

#include "include_deal.h"
dealii::Threads::Mutex mutex_deal;
#include "Basics.h"
#include "exact_solution.h"
#include "equation_gen.h"
#include "mesh_gen.h"
#include "PostProc.h"
#include "Solver.h"
#include <time.h>

using namespace dealii;
using namespace std;
using namespace SolverDG;

int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,num_threads);

	const int dim = 2;
	const int refine_cycles = 5;			// total number of refinement cycles
	const int p = 2;
	const string mesh_file_name = "mesh/mesh_file";
	const unsigned int mapping_order = 2;
	
	enum System_Types
	{systemA,systemB};

	System_Types system_type = systemA;

	switch(system_type)
	{
		case systemA:
		{
			const unsigned int nEqnA = 6;
			exact_solutionA<dim> exact_solutionA(nEqnA);

			Solver_DG<dim> solver(p,mapping_order,Solver_DG<dim>::global,&exact_solutionA);
			solver.run(mesh_file_name,refine_cycles);

			break;
		}

		case systemB:
		{
			const unsigned int nEqnB = 10;
			exact_solutionB<dim> exact_solutionB(nEqnB);

			Solver_DG<dim> solver(p,mapping_order,
								 Solver_DG<dim>::global,&exact_solutionB);

			solver.run(mesh_file_name,refine_cycles);

			break;
		}
	}


}
	
