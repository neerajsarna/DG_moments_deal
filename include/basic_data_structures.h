
// meshworker works on distributed memory type
enum Assembly_Type
{meshworker, manuel};

enum System_Type
{symmetric,un_symmetric};

enum Num_Flux
{Upwind, LLF};

enum Force_Type
{type1,type2,type3};

enum Refinement
    {global,adaptive,adaptive_kelly};

enum BC_type
	{characteristic,odd};

// ring domain for systemA and systemB
// periodic_square for head conduction and Couette flow problem
enum Mesh_type
	{ring,periodic_square};

// how to generate using dealii or generate using gmsh
enum Meshing_Options 
	{read_msh, generate_internal};

struct numerical_data
{
	int p;
	int mapping_order;
	int refine_cycles;
	Refinement refinement;
	Assembly_Type assembly_type;
};

struct nEqn_data
{
	unsigned int no_of_total_systems;
	unsigned int system_to_solve;

	std::vector<unsigned int> total_nEqn;
	std::vector<unsigned int> nBC;
	std::map<unsigned int,char> system_id_nEqn;			// maps the system_id with the number of equations
	
	Force_Type force_type;
	System_Type system_type;			
	BC_type bc_type;
};

struct physical_data
{
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
};

/*storing data related to tensors*/
struct tensor_data
{
	std::vector<unsigned int> free_index_info;
};

struct mesh_data
{
	Mesh_type mesh_type;
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

};