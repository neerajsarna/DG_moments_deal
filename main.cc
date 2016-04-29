#include <iostream>
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include <fstream>
#include <cmath>

#include "include_deal.h"
dealii::Threads::Mutex mutex_deal;
#include "Basics.h"
#include "equation_gen.h"
#include "exact_solution.h"
#include "mesh_gen.h"
#include "PostProc.h"
#include "Solver.h"
#include <time.h>

using namespace dealii;
using namespace std;


int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,num_threads);
	//MultithreadInfo::set_thread_limit(num_threads);

	const int dim = 2;
	const int refine_cycles = 7;			// total number of refinement cycles
	const int p = 2;
	const string mesh_file_name = "mesh/mesh_file";
	const unsigned int mapping_order = 2;
	
	nEqn_data num_equations;
	num_equations.no_of_total_systems = 1;
	num_equations.total_nEqn.resize(num_equations.no_of_total_systems);
	num_equations.total_nEqn[0] = 6;
	num_equations.system_id_nEqn[0] = 'A';

	assert(num_equations.no_of_total_systems == 1);
	const unsigned int solve_system = 0;				// id of the system we wish to solve

	const System_Type system_type = un_symmetric;
	const Num_Flux num_flux = Upwind;

	EquationGenerator::Base_EquationGenerator<system_type,num_flux,dim> system_of_equations(num_equations);

	ExactSolution::Base_ExactSolution<system_type,num_flux,dim> exact_solution(0,6);
	SolverDG::Solver_DG<system_type,num_flux,dim> solver(p,mapping_order,
															SolverDG::Solver_DG<system_type,num_flux,dim>::global,
															&system_of_equations,
															&exact_solution,
															solve_system);

	solver.run(mesh_file_name,refine_cycles);

}
	
