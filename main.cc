#include "include_deal.h"
#include "gtest/gtest.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include "parse_command_line.h"
#include "mkl.h"

// The input file containing the computational parameters
std::string input_file;

// core routines
#include "Base_Constants.h"
#include "Base_TensorInfo.h"
#include "Base_ForceType.h"
#include "Base_MatrixOpt.h"
#include "Base_BoundaryHandler.h"
#include "Base_BCrhs.h"
#include "BCrhs_systemA.h"
#include "Base_ExactSolution.h"
#include "ExactSolution_PoissonHeatConduction_G26.h"
#include "ExactSolution_SystemA.h"
#include "Base_EquationGenerator.h"
#include "LinearSolver.h"
#include "SystemA.h"
#include "G26.h"
#include "Base_MeshGenerator.h"
#include "Base_PeriodicityHandler.h"
#include "Base_PostProc.h"
#include "Base_Solver.h"

using namespace dealii;

int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);
	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv, num_threads);

	parse_command_line commmand_line_parser(argc,argv);

	/*store the filename read from the command line*/
	commmand_line_parser.read_input_filename(input_file);

	std::string folder_name = "../system_matrices/";
	const unsigned int dim = 2;
	Constants::Base_Constants constants(input_file);

	// division based upon the number of equations we are solving
	switch(constants.constants.nEqn)
	{
		case 6:
		{
			// first develop the class which loads the systems
			SystemA::SystemA<dim> systemA(constants.constants,folder_name);			

			// routine to develop the exact solution
			ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

			// the solver which combines the above two classes
			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &systemA,
													 &exact_solution_systemA);

			// now run the simulations 
			switch(systemA.constants.mesh_type)
			{
				case ring:
				{
					base_solver.run_ring();
					break;
				}

				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
				
			}
			break;
		}

		case 17:
		{
			G26::G26<dim> G26(constants.constants,folder_name);

			ExactSolution::G26_PoissonHeat<dim>  G26_PoissonHeat(constants.constants,G26.base_tensorinfo.S_half);

			FEM_Solver::Base_Solver<dim> base_solver("mesh_file_name",
													 "grid",
													 constants.constants,
													 &G26,
													 &G26_PoissonHeat);

			switch(G26.constants.mesh_type)
			{
				case periodic_square:
				{
					base_solver.run_periodic();
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			break;
		}

		default:
		{
			Assert(1 == 0,ExcMessage("Should not have reached here"));
			break;
		}
	}



	return 0;
}
