#include <iostream>
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include <fstream>
#include <cmath>

#include "include_deal.h"
#include "mkl.h"
#include "parse_command_line.h"
#include "input_parameters.h"
#include "Basics.h"
#include "equation_gen.h"
#include "exact_solution.h"
#include "mesh_gen.h"
#include "PostProc.h"
#include "Solver.h"

using namespace dealii;
using namespace std;


int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,num_threads);

	parse_command_line commmand_line_parser(argc,argv);
	string input_file_name;
	commmand_line_parser.read_input_filename(input_file_name);

	Input_parameters::parameter_handling container_parameter;
	container_parameter.read_input_file(input_file_name);

	const int dim = 2;
	string mesh_file_name = "mesh/mesh_file";
	string output_dir_name;

	numerical_data numerical_constants;
	nEqn_data num_equations;
	physical_data physical_constants;

	container_parameter.read_numerical_constants(numerical_constants);
	container_parameter.read_system_data(num_equations);
	container_parameter.read_physical_constants(physical_constants);
	container_parameter.read_output_dir_name(output_dir_name);

	const Num_Flux num_flux = Upwind;

	EquationGenerator::Base_EquationGenerator<num_flux,dim> system_of_equations(num_equations,
																				physical_constants,
																				output_dir_name);

	ExactSolution::Base_ExactSolution<dim> exact_solution(num_equations.system_to_solve,
															num_equations.total_nEqn[num_equations.system_to_solve],
															system_of_equations.system_data[num_equations.system_to_solve].S_half,
															num_equations.system_type,
															physical_constants,
															output_dir_name);

	SolverDG::Solver_DG<num_flux,dim> solver(numerical_constants.p,
									         numerical_constants.mapping_order,
											 numerical_constants.refinement,
											  &system_of_equations,
											  &exact_solution,
											 num_equations.system_to_solve,
											 physical_constants,
											 output_dir_name);

	solver.run(mesh_file_name,numerical_constants.refine_cycles);
}
	
