// meshworker works on distributed memory type
enum Assembly_Type
{meshworker, manuel};

enum Num_Flux
{Upwind, LLF};

enum Force_Type
{type1,type2,type3,type4};

enum Refinement
{global,adaptive,adaptive_kelly};

enum BC_Type
{characteristic,odd};

// types of mesh which will be generated 
enum Mesh_type
{ring,square_domain,square_circular_cavity,NACA5012,line,box_cylinder};

// types of test cases
enum Problem_type
{heat_conduction,inflow_outflow,periodic,lid_driven_cavity};

// how to generate using dealii or generate using gmsh
enum Meshing_Options 
{read_msh, generate_internal};

enum Collision_Operator
{BGK,Boltzmann_MM};


/*storing data related to tensors*/
struct tensor_data
{
	std::vector<unsigned int> free_index_info;
};

	// details regarding the set of equations
	struct constant_system
	{
		// data set for system info
		int total_systems;
		
		// we consider vectors so that we can load more than one system
		std::vector<int> nEqn;
		std::vector<int> nBC;
		std::vector<int> system_id;
		std::vector<int> Ntensors;

	};


	// details regarding the numerical constants
	struct constant_numerics
	{
		int p;
		int mapping_order;
		int refine_cycles;
		int refine_cycles_c;
		double fraction_refine;
		Refinement refinement;
		Assembly_Type assembly_type;
	
		Force_Type force_type;
		BC_Type bc_type;		

				// physical constants
		double tau;
		double zeta;
		double chi;
		double uW;
		double A0;
		double A1;
		double A2;
		double kappa;

		// coefficient for normal relaxational velocity
		double epsilon;

		// variable in which error has to be found
		std::string error_variable;

	// the equation in which force has to be applied
		std::string force_variable;

		// forcing term for poisson heat conduction
		double alpha;

			// data for the boundary
	// first wall
		double theta0;
		double vx0;
		double vy0;
		double vz0;


	// data for the second wall
		double theta1;
		double vx1;
		double vy1;
		double vz1;

	// data for the third wall
		double theta2;
		double vx2;
		double vy2;
		double vz2;

	// data for the fourth wall
		double theta3;
		double vx3;
		double vy3;
		double vz3;

	// data for the fourth wall
		double theta4;
		double vx4;
		double vy4;
		double vz4;

	// first inflow boundary
		double rho101;
		double theta101;
		double vx101;
		double vy101;
		double vz101;

	// second inflow boundary
		double rho102;
		double theta102;
		double vx102;
		double vy102;
		double vz102;
	
		// the type of mesh to be used
		Mesh_type mesh_type;

	// type of test case 
		Problem_type problem_type;

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

	// inner and outer radius for the ring geometry
		double outer_radius;
		double inner_radius;

	// initial refinement level
		unsigned int initial_refinement;



	// names of the sub directories
		std::vector<std::string> sub_directory_names;

	// output directory name
		std::string main_output_dir;

	// decide whether to print for all the refine steps or not
		bool print_all;

	// decide whether to print the solution or not
		bool print_solution;

	// decide whether to print error or not
		bool print_error;

	// decide whether to print exact solution or not
		bool print_exactsolution;

	// decide whether to print the convergence table or not
		bool print_convergence_table;

		bool print_velocity_space_error;

		bool print_fe_index;


	// we create a map between the id of the variable and it's name
	// for every dimension we have a different variable map
		std::vector<std::map<std::string,unsigned int>> variable_map;
	//std::map<std::string,unsigned int> variable_map_1D;

		Collision_Operator coll_op;

	};

