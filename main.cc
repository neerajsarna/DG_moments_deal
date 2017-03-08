#include "include_deal.h"
#include "EigenSetup.h"
#include "basic_data_structures.h"
#include "parse_command_line.h"
#include "mkl.h"

// The input file containing the computational parameters
std::string input_file;

#include "include_files.h"

using namespace dealii;

int main(int argc,char **argv)
{
	const unsigned int num_threads = atoi(argv[1]);
	dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,argv, num_threads);

	parse_command_line commmand_line_parser(argc,argv);

	/*store the filename read from the command line*/
	commmand_line_parser.read_input_filename(input_file);

	//Run the simulations 	
	
	Run();

	return 0;
}
