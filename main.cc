#include <iostream>
#include "EigenSetup.h"
#include <fstream>
#include <cmath>

#include "include_deal.h"
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
	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, ::numbers::invalid_unsigned_int);
	clock_t start_t,end_t;

	const int dim = 2;
	const int refine_cycles = 4;			// total number of refinement cycles
	const int p = 1;
	const string mesh_file_name = "mesh/mesh_file";
	const unsigned int mapping_order = 2;

	start_t = clock();

	Solver_DG<dim> solver(p,mapping_order,Solver_DG<dim>::global);

	solver.run(mesh_file_name,refine_cycles);

	end_t = clock();

	cout << ".......Total time taken: " << (double)(end_t-start_t)/(CLOCKS_PER_SEC) << endl;

}
	
