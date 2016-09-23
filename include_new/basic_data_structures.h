
// meshworker works on distributed memory type
enum Assembly_Type
{meshworker, manuel};

enum System_Type
{symmetric,un_symmetric};

enum Num_Flux
{Upwind, LLF};

enum Force_Type
{type1,type2,type3,type4};

enum Refinement
{global,adaptive,adaptive_kelly};

enum BC_Type
{characteristic,odd};

// ring domain for systemA and systemB
// periodic_square for head conduction and Couette flow problem
enum Mesh_type
{ring,periodic_square};

// how to generate using dealii or generate using gmsh
enum Meshing_Options 
{read_msh, generate_internal};


/*storing data related to tensors*/
struct tensor_data
{
	std::vector<unsigned int> free_index_info;
};


struct constant_data
{
	int p;
	int mapping_order;
	int refine_cycles;
	Refinement refinement;
	Assembly_Type assembly_type;

		// data set for system info
	int nEqn;
	int nBC;
	int system_id;

	System_Type system_type;
	Force_Type force_type;
	BC_Type bc_type;

		// physical constants
	double tau;
	double zeta;
	double chi;
	double theta0;
	double theta1;
	double uW;
	double A0;
	double A1;
	double A2;
	double kappa;
	
		// coefficient for normal relaxational velocity
	double epsilon;

		// variable in which error has to be found
	std::string error_variable;

		// forcing term for poisson heat conduction
	double alpha;

		// the type of mesh to be used
	Mesh_type mesh_type;

		// different options for meshing
	Meshing_Options mesh_options;
	std::string mesh_filename;

	// to be only used for square domain
	// xcord of left boundary
	double xl;

	// xcord of the right boundary
	double xr;

	// ycord of the bottom boundary
	double yb;

	// ycord of the top boundary
	double yt;

	// number of partisions per dimension for the square domain
	unsigned int part_x;

	unsigned int part_y;
};

