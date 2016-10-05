// the following namespace just reads all the parameters
namespace Constants
{
	using namespace dealii;
	using namespace std;
	

	class Base_Constants
	{
	public:
		Base_Constants(string &input_file);

		// object which handles the parameters
		ParameterHandler prm;

			// the following routine prepares the input file to be read
		void read_input_file(const string &input_file_name);

			// the following routine reads the numerical parameters 
		void read_numerical_constants();

		// read the system data
		void read_system_data();

		// read physical constants
		void read_physical_constants();

		// read mesh info
		void read_mesh_info();

		// read output director
		void read_output_directory();

		// read printing options
		void read_printing_options();

		// set of numerical constants
		constant_data constants;
		private:
			void declare_parameters();
			void allocate_variable_map();
			void allocate_subdirector_names();
	};

	Base_Constants::
	Base_Constants(string &input_file)
	{
		// first we declare all the parameters which could be found in the input file
		declare_parameters();

		// then we prepare the reading of the input file
		read_input_file(input_file);

		// now we read the numerical constants 
		read_numerical_constants();

		// now we read the data corresponding to the system
		// like the number of equations, the number of boundary conditions etc. 
		read_system_data();

		// read the physical constants
		read_physical_constants();

		// read the mesh info
		read_mesh_info();

		// read output directory name
		read_output_directory();

		// printing options
		read_printing_options();

		// we create a map between the id of the variable and it's name,
		// proves to be helpful during error computation
		allocate_variable_map();

		allocate_subdirector_names();
	}

	void Base_Constants
	::read_input_file(const string &input_file)
	{
		prm.read_input(input_file.c_str());
	}

	void Base_Constants
	::read_numerical_constants()
	{
		bool entered = false;
		prm.enter_subsection("Numerical Constants");
		{
			entered = true;
			constants.p = prm.get_integer("polynomial degree");
			constants.mapping_order = prm.get_integer("mapping order");
			constants.refine_cycles = prm.get_integer("total h refinement cycles");
			constants.refinement = (Refinement)prm.get_integer("type of refinement");
			constants.assembly_type = (Assembly_Type)prm.get_integer("assembly type");
		}
		prm.leave_subsection();

		Assert(entered, ExcMessage("did not read numerical constants"));
	}


	void Base_Constants
	::declare_parameters()
	{
		prm.enter_subsection("Numerical Constants");
		{
			prm.declare_entry("polynomial degree",
				"1",
				Patterns::Integer(0,10),
				"polynomial order per dimension");

			prm.declare_entry("mapping order",
				"1",
				Patterns::Integer(1,10),
				"mapping order for isoparametrics");

			prm.declare_entry("total h refinement cycles",
				"1",
				Patterns::Integer(1,10),
				"total number of h refinement cycles");

			prm.declare_entry("type of refinement",
				"0",
				Patterns::Integer(0,5),
				"Type of refinement to be used");

				// default is meshworker
			prm.declare_entry("assembly type",
				"0",
				Patterns::Integer(0,1),
				"type of assembly");
		}
		prm.leave_subsection();

		prm.enter_subsection("Physical Constants");
		{
			prm.declare_entry("tau",
				"0.1",
				Patterns::Double(0.001,10),
				"The knudsen number");

			prm.declare_entry("zeta",
				"1.0",
				Patterns::Double(0.0,20.0),
				"Boundary condition handling");


			prm.declare_entry("chi",
				"1.0",
				Patterns::Double(0.0,20.0),
				"Boundary condition handling");
			
				// could be negative if one is using different variables
			prm.declare_entry("theta0",
				"2.0",
				Patterns::Double(-20.0,20.0),
				"Temperature of inner wall");


			prm.declare_entry("theta1",
				"1.0",
				Patterns::Double(-20.0,20.0),
				"Temperature of outer wall");


			prm.declare_entry("uW",
				"0.1",
				Patterns::Double(0.0,10.0),
				"Wall velocity");



			prm.declare_entry("A0",
				"0.0",
				Patterns::Double(0,10.0),
				"Coefficient in the force");

			prm.declare_entry("A1",
				"0.0",
				Patterns::Double(0.0,20.0),
				"Coefficient in the force");

			prm.declare_entry("A2",
				"0.1",
				Patterns::Double(0.0,10.0),
				"Coefficient in the force");

			prm.declare_entry("kappa",
				"0.0",
				Patterns::Double(0,100.0),
				"Coefficient for the boundary");

			prm.declare_entry("epsilon",
				"0.00001",
				Patterns::Double(0,0.0001),
				"Coefficient for normal relaxational velocity");

			prm.declare_entry("error_variable",
				"theta",
				Patterns::DirectoryName(),
				"name of the variable for error evaluation");		

			prm.declare_entry("force_variable",
							  "rho",
							  Patterns::DirectoryName(),
							"name of the variable to which external force has to be applied");		

			prm.declare_entry("alpha",
				"0",
				Patterns::Double(0,10.0),
				"forcing term poisson heat conduction");										  		
		}
		prm.leave_subsection();

		prm.enter_subsection("Output Directory Name");
		{
			prm.declare_entry("Main Output Directory",
				"outputs",
				Patterns::DirectoryName(),
				"Directory name for writting results");
		}
		prm.leave_subsection();

		prm.enter_subsection("System Properties");
		{


			prm.declare_entry("system_id",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm.declare_entry("equations in the system",
				"6",
				Patterns::Integer(6,20),
				"total number of equations in system");

			prm.declare_entry("nBC",
				"2",
				Patterns::Integer(2,10),
				"total number of boundary conditions needed for the system");

				/*0 for force of system A and 1 for force of system B*/
			prm.declare_entry("force type",
				"0",
				Patterns::Integer(0,4),
				"Type of force on the system");

				/*0 for characteristic and 1 for odd*/
			prm.declare_entry("BC type",
				"0",
				Patterns::Integer(0,1),
				"Type of BC");
		}
		prm.leave_subsection();

		prm.enter_subsection("Mesh Info");
		{
			prm.declare_entry("filename",
				"mesh/mesh_file",
				Patterns::DirectoryName(),
				"name of mesh file");

				/*0 for ring and 1 for periodic square*/
			prm.declare_entry("mesh type",
				"0",
				Patterns::Integer(0,1),
				"type of mesh");

				/*0 to read from gmsh
				  1 to generate using deal*/
			prm.declare_entry("mesh options",
				"0",
				Patterns::Integer(0,1),
				"how to generate the mesh");

				// the left and right boundary
			prm.declare_entry("left boundary",
				"-1",
				Patterns::Double(-10.0,-0.1),
				"left boundary");

			prm.declare_entry("right boundary",
				"1",
				Patterns::Double(0.1,10.0),
				"right boundary");

			prm.declare_entry("bottom boundary",
				"-1",
				Patterns::Double(-10.0,-0.1),
				"bottom boundary of square domain");

			prm.declare_entry("top boundary",
				"1",
				Patterns::Double(0.1,10.0),
				"top boundary of square domain");	

			prm.declare_entry("part_x",
				"2",
				Patterns::Integer(1,50),
				"parts in x");

			prm.declare_entry("part_y",
				"2",
				Patterns::Integer(1,1000),
				"parts in y");

		}
		prm.leave_subsection();

		prm.enter_subsection("Printing Options");
		{
			prm.declare_entry("print_all",
							"false",
							Patterns::Anything(),
							"decide to print when");

			prm.declare_entry("print_solution",
							  "false",
							  Patterns::Anything(),
							  "decide to print solution");

			prm.declare_entry("print_error",
							  "false",
							  Patterns::Anything(),
							  "decide to print error");

			prm.declare_entry("print_exactsolution",
							  "false",
							  Patterns::Anything(),
							  "decide to print exact solution");

			prm.declare_entry("print_convergence_table",
							  "true",
							  Patterns::Anything(),
							  "decide to print the convergence table");
		}
		prm.leave_subsection();

	}


	void Base_Constants
	::read_system_data()
	{
		bool entered = false;
		prm.enter_subsection("System Properties");
		{
			entered = true;

			constants.nEqn = prm.get_integer("equations in the system");
			constants.nBC = prm.get_integer("nBC");
			constants.system_id = prm.get_integer("system_id");

			constants.force_type = (Force_Type)prm.get_integer("force type");
			constants.bc_type = (BC_Type)prm.get_integer("BC type");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read system data"));
	
	}


	void Base_Constants
	::read_physical_constants()
	{
		bool entered = false;
		prm.enter_subsection("Physical Constants");
		{
			entered = true;
			constants.tau = prm.get_double("tau");
			constants.zeta = prm.get_double("zeta");
			constants.chi = prm.get_double("chi");
			

			constants.theta0 =  prm.get_double("theta0");
		 	constants.theta1  = prm.get_double("theta1");
			constants.uW = prm.get_double("uW");
			constants.A0 = prm.get_double("A0");
		    constants.A1 = prm.get_double("A1");
			constants.A2 = prm.get_double("A2");
			constants.kappa  = prm.get_double("kappa");
			constants.epsilon = prm.get_double("epsilon");
			constants.error_variable = prm.get("error_variable");
			constants.force_variable = prm.get("force_variable");
			constants.alpha = prm.get_double("alpha");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read physical constants"));
	}


	void Base_Constants
	::read_mesh_info()
	{
		bool entered = false;

		prm.enter_subsection("Mesh Info");
		{
			entered = true;
			constants.mesh_filename = prm.get("filename");
			constants.mesh_type = (Mesh_type)prm.get_integer("mesh type");
			constants.mesh_options = (Meshing_Options)prm.get_integer("mesh options");
			constants.xl = prm.get_double("left boundary");
			constants.xr = prm.get_double("right boundary");
			constants.yt = prm.get_double("top boundary");
			constants.yb = prm.get_double("bottom boundary");
			constants.part_x = prm.get_integer("part_x");
			constants.part_y = prm.get_integer("part_y");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read mesh info"));
	}

	void Base_Constants
	::read_output_directory()
	{
		bool entered = false;

		// enter the subsection which reads the output directory name
		prm.enter_subsection("Output Directory Name");
		{
			entered = true;

			// input the output directory name
			constants.main_output_dir = prm.get("Main Output Directory");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("Did not read the output directory name"));
	}

	void Base_Constants
	::allocate_variable_map()
	{

	  	Assert(constants.nEqn <= 17,ExcNotImplemented());
	  	std::vector<std::string> var_names = {"rho","vx","vy","theta","sigmaxx","sigmaxy","sigmayy","qx","qy","mxxx","mxxy","mxyy","myyy",
	  								"Delta","Rxx","Rxy","Ryy"};

	  	for (unsigned int i = 0 ; i < 17 ; i++)
	  		constants.variable_map[var_names[i]] = i;
	  
	}

	void Base_Constants
	::allocate_subdirector_names()
	{
		const unsigned int num_outputs = 5;
		constants.sub_directory_names.resize(num_outputs);

		constants.sub_directory_names[0] = constants.main_output_dir + "/grids";

		
		constants.sub_directory_names[1] = constants.main_output_dir + "/solution";

		
		constants.sub_directory_names[2] = constants.main_output_dir + "/convergence_tables";

		
		constants.sub_directory_names[3] = constants.main_output_dir + "/sparsity_patterns";

		
		constants.sub_directory_names[4] = constants.main_output_dir + "/matrix_read";			
		
	}

	void Base_Constants
	::read_printing_options()
	{
		bool entered = false;

		prm.enter_subsection("Printing Options");
		{
			constants.print_all = prm.get_bool("print_all");
			constants.print_solution = prm.get_bool("print_solution");
			constants.print_error = prm.get_bool("print_error");
			constants.print_exactsolution = prm.get_bool("print_exactsolution");
			constants.print_convergence_table = prm.get_bool("print_convergence_table");
		}
		prm.leave_subsection();
	}
}