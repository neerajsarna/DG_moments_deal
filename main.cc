// Running the file
//1. cmake .
//2. make debug or make release
//3. cd "directory for mkl"
//4. source ./mklvars.sh intel64 mod ilp64
//5. in the current directory ./mycode.out #threads -p "location of the input file"
#include <iostream>
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include <fstream>
#include <cmath>

/*The following file contains all
 the header files needed from deal*/
#include "include_deal.h"
#include "mkl.h"
/*Used to order the eigen system obtained from EigenSolver.
 Ordering is done in ascending order*/
#include "ordering.h"
/*Get command line arguments,
 the number of processors and the name of the input file*/
#include "parse_command_line.h"

/*file for handling
 the input parameters*/
#include "input_parameters.h"

// develops all the info related to tensors
#include "Tensors_Info.h"


/*contains the most basic of the code i.e.
 the variables like the wall velocity and the names of the folder and 
subdirectoreis*/
#include "Basics.h"
/*generates all the data 
required for the equations*/
#include "equation_gen.h"
/*develops the exact solution for
 different systems generated in the above equation*/
#include "exact_solution.h"

/*functions for generating
 and printing the mesh*/
#include "mesh_gen.h"
/*abstract class for
 post processing*/
#include "PostProc.h"
/*combines all the above 
class to solve the resulting system*/
#include "Solver.h"

using namespace dealii;
using namespace std;


int main(int argc,char **argv)
{
	/*maximum number of threads to be used
	 by the program. Does not control the number of threads for Pardiso,
	 enter the number of threads for pardiso from export OMP_NUM_THREADS*/
	const unsigned int num_threads = atoi(argv[1]);

	/*Though the whole code does not exploit distributed memory but we still need
	 to initialize MPI for trilinos*/
	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,num_threads);

	/*read the command line arguments*/
	parse_command_line commmand_line_parser(argc,argv);
	string input_file_name;
	/*store the filename read from the command line*/
	commmand_line_parser.read_input_filename(input_file_name);

	/*Class which reads all the parameters*/
	Input_parameters::parameter_handling container_parameter;
	
	/*initialize the reading from the input file*/
	container_parameter.read_input_file(input_file_name);

	


	/*number of dimensions in the system*/
	const int dim = 2;
	

	/*main output directory name, to be read from the input file*/
	string output_dir_name;

	/*all the data related to mapping order, degree of the polynomial etc.*/
	numerical_data numerical_constants;
	nEqn_data num_equations;
	physical_data physical_constants;
	tensor_data tensor_info;
	mesh_data mesh_info;

	/*reading of all the parameters,
	  will be used by the future classes*/
	container_parameter.read_numerical_constants(numerical_constants);
	container_parameter.read_system_data(num_equations);
	container_parameter.read_physical_constants(physical_constants);
	container_parameter.read_output_dir_name(output_dir_name);
	container_parameter.read_mesh_info(mesh_info);

	// given the tensorial degree of the tensor, the following routine allocates the total number of 
	// independent variables.
	Tensor_info::tensor_development<dim> develop_tensors;
	develop_tensors.map_free_indices_to_variable(tensor_info);	

	// Type of numerical flux, Upwind and LLF available
	const Num_Flux num_flux = Upwind;

	// generate the data corresponding to different to the equations to be solved
	EquationGenerator::Base_EquationGenerator<num_flux,dim> system_of_equations(num_equations,
																				physical_constants,
																				tensor_info,
																				mesh_info,
	 																			output_dir_name);

	// due to the availability of exact solution the implementation differs
	switch (mesh_info.mesh_type)
	{
		// the exact solution is only known for a ring geometry
		case ring:
		{
				ExactSolution::Base_ExactSolution<dim> exact_solution(num_equations.system_to_solve,
															num_equations.total_nEqn[num_equations.system_to_solve],
															system_of_equations.system_data[num_equations.system_to_solve].S_half,
															num_equations.system_type,
															physical_constants,
															output_dir_name);


				SolverDG::Solver_DG<num_flux,dim> solver(numerical_constants,
											  			 &system_of_equations,
											  			 &exact_solution,
											  			 num_equations,
											  			 physical_constants,
											  			 output_dir_name,
											  			 mesh_info);


				solver.run(numerical_constants.refine_cycles);
				
			break;
		}

		// no known exact solution for a periodic square
		case periodic_square:
		{
			switch(num_equations.system_id[num_equations.system_to_solve])
			{
				//G26A system
				case 6:
				{
					ExactSolution::systemA_period_sqr<dim> exact_solution(num_equations.system_to_solve,
																	num_equations.total_nEqn[num_equations.system_to_solve],
																	system_of_equations.system_data[num_equations.system_to_solve].S_half,
																	num_equations.system_type,
																	physical_constants,
																	output_dir_name);

			    	SolverDG::Solver_DG<num_flux,dim> solver(numerical_constants,
											  			 &system_of_equations,
											  			 &exact_solution,
											  			 num_equations,
											  			 physical_constants,
											  			 output_dir_name,
											  			 mesh_info);

					solver.run(numerical_constants.refine_cycles);

					break;
				}

				// the G26 system
				case 17:
				{
					ExactSolution::G26_period_sqr<dim> exact_solution(num_equations.system_to_solve,
																	num_equations.total_nEqn[num_equations.system_to_solve],
																	system_of_equations.system_data[num_equations.system_to_solve].S_half,
																	num_equations.system_type,
																	physical_constants,
																	output_dir_name);

					SolverDG::Solver_DG<num_flux,dim> solver(numerical_constants,
											  			 &system_of_equations,
											  			 &exact_solution,
											  			 num_equations,
											  			 physical_constants,
											  			 output_dir_name,
											  			 mesh_info);

					solver.run(numerical_constants.refine_cycles);
					break;
				}

				default:
				{
				
					Assert(1 == 0 , ExcNotImplemented());
					break;
				}
			}

	
			break;
		}

		default:
		{
			Assert(1 == 0, ExcNotImplemented());
			break;
		}
	}


	
}
	
