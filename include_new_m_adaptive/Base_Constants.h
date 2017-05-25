// the following namespace just reads all the parameters
namespace Constants
{
	using namespace dealii;
	using namespace std;
	

	class Base_Constants
	{
	public:
		Base_Constants(std::vector<std::string> &input_file);

		// object which handles all the parameters apart from the system info
		ParameterHandler prm;

		// parameter handler which handles the data corresponding to the system
		ParameterHandler prm_system_info;

			// the following routine prepares the input file to be read
		void read_input_file(std::vector<string> &input_file_name);

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
		constant_numerics constants_num;
		constant_system constants_sys;
		private:
			void declare_parameters();
			void allocate_variable_map();
			void allocate_variable_map_1D();
			void allocate_subdirector_names();
	};

	Base_Constants::
	Base_Constants(std::vector<string> &input_file)
	{
		// we currently expect two input file names
		AssertDimension(input_file.size(),2);

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

		allocate_variable_map_1D();

		allocate_subdirector_names();
	}

	void Base_Constants
	::read_input_file(std::vector<std::string> &input_file)
	{
		AssertDimension(input_file.size(),2);
		//prm.parse_input(input_file.c_str());

		// first file is for the system info
		prm_system_info.read_input(input_file[0].c_str());
		prm.read_input(input_file[1].c_str());
	}

	void Base_Constants
	::read_numerical_constants()
	{
		bool entered = false;
		prm.enter_subsection("Numerical Constants");
		{
			entered = true;
			constants_num.p = prm.get_integer("polynomial degree");
			constants_num.mapping_order = prm.get_integer("mapping order");
			constants_num.refine_cycles = prm.get_integer("total h refinement cycles");
			constants_num.refine_cycles_c = prm.get_integer("refine_cycles_c");
			constants_num.refinement = (Refinement)prm.get_integer("type of refinement");
			constants_num.assembly_type = (Assembly_Type)prm.get_integer("assembly type");
		}
		prm.leave_subsection();

		Assert(entered, ExcMessage("did not read numerical constants"));
	}


	void Base_Constants
	::declare_parameters()
	{
		// first we declare entries for the system
		prm_system_info.enter_subsection("System Properties");
		{
			prm_system_info.declare_entry("total_systems",
										  "1",
										  Patterns::Integer(1,10),
										  "total systems to be loaded");

			prm_system_info.declare_entry("Ntensors0",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id0",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations0",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC0",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			prm_system_info.declare_entry("Ntensors1",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id1",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations1",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC1",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			prm_system_info.declare_entry("Ntensors2",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id2",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations2",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC2",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			//*****************
			prm_system_info.declare_entry("Ntensors3",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id3",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations3",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC3",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			//////********************
			prm_system_info.declare_entry("Ntensors4",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id4",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations4",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC4",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			//*********************
			prm_system_info.declare_entry("Ntensors5",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id5",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations5",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC5",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

			//**************************
			prm_system_info.declare_entry("Ntensors6",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id6",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations6",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC6",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");


			//*************************
				prm_system_info.declare_entry("Ntensors7",
				"2",
				Patterns::Integer(2,100),
				"number of tensors");

			prm_system_info.declare_entry("system_id7",
				"6",
				Patterns::Integer(-100,100),
				"system ID(distinguish between regularized and normal theories)");

			prm_system_info.declare_entry("num_equations7",
				"3",
				Patterns::Integer(3,100),
				"total number of equations in system");

			prm_system_info.declare_entry("nBC7",
				"2",
				Patterns::Integer(2,100),
				"total number of boundary conditions needed for the system");

		}
		prm_system_info.leave_subsection();

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
				Patterns::Integer(1,100),
				"total number of h refinement cycles");

			// total number of refinement cycles in the velocity space
			prm.declare_entry("refine_cycles_c",
							"1",Patterns::Integer(0,100),
							"total number of refinement cycles in the velocity space");

			prm.declare_entry("type of refinement",
				"0",
				Patterns::Integer(0,5),
				"Type of refinement to be used");

				// default is meshworker
			prm.declare_entry("assembly type",
				"0",
				Patterns::Integer(0,1),
				"type of assembly");

			prm.declare_entry("matrix_type",
				"0",
				Patterns::Integer(0,1),
				"type of global matrix");
		}
		prm.leave_subsection();

		prm.enter_subsection("Physical Constants");
		{

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

			prm.declare_entry("theta0",
				"0",
				Patterns::Double(-100.0,100.0),
				"temperature");								  		

			prm.declare_entry("vx0",
				"0",
				Patterns::Double(0,100.0),
				"normal velocity");

			prm.declare_entry("vy0",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");

			prm.declare_entry("theta1",
				"0",
				Patterns::Double(-100.0,100.0),
				"temperature of the second wall");								  		

			prm.declare_entry("vx1",
				"0",
				Patterns::Double(0,100.0),
				"normal velocity");

			prm.declare_entry("vy1",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");


			prm.declare_entry("theta2",
				"0",
				Patterns::Double(0,100.0),
				"temperature");								  		

			prm.declare_entry("vx2",
				"0",
				Patterns::Double(0,100.0),
				"normal velocity");

			prm.declare_entry("vy2",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");

			prm.declare_entry("theta3",
				"0",
				Patterns::Double(0,100.0),
				"temperature");								  		

			prm.declare_entry("vx3",
				"0",
				Patterns::Double(0,100.0),
				"normal velocity");

			prm.declare_entry("vy3",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");


			prm.declare_entry("theta4",
				"0",
				Patterns::Double(0,100.0),
				"temperature");								  		

			prm.declare_entry("vx4",
				"0",
				Patterns::Double(0,100.0),
				"normal velocity");

			prm.declare_entry("vy4",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");

			prm.declare_entry("rho101",
				"0",
				Patterns::Double(0,10.0),
				"temperature");								  		

			prm.declare_entry("theta101",
				"0",
				Patterns::Double(0,100.0),
				"temperature");								  		

			prm.declare_entry("vx101",
				"0",
				Patterns::Double(-100,100.0),
				"normal velocity");

			prm.declare_entry("vy101",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");


			prm.declare_entry("rho102",
				"0",
				Patterns::Double(0,10.0),
				"temperature");								  		

			prm.declare_entry("theta102",
				"0",
				Patterns::Double(0,100.0),
				"temperature");								  		

			prm.declare_entry("vx102",
				"0",
				Patterns::Double(-100,100.0),
				"normal velocity");

			prm.declare_entry("vy102",
				"0",
				Patterns::Double(-100.0,100.0),
				"tangential velocity");

			prm.declare_entry("Collision_Operator",
											"0",
										Patterns::Integer(0,1),
										"Type of Collision Operator");
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


		prm.enter_subsection("Mesh Info");
		{
			// filename for the mesh
			prm.declare_entry("mesh_filename",
							"mesh/mesh_file",
				Patterns::DirectoryName(),
				"name of mesh file");

				/*0 for ring, 
				  1 square,
				  2 square with a circle
				  3 for NACA5012 in a channel 
				  4 for a line in 1D
				*/

			prm.declare_entry("mesh type",
				"0",
				Patterns::Integer(0,4),
				"type of mesh");

			/*0 heat_conduction,
			  1 inflow_outflow,
			  2 periodic,
			  3 lid_driven_cavity*/

			prm.declare_entry("problem type",
				"0",
				Patterns::Integer(0,4),
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
				Patterns::Integer(1,1000),
				"parts in x");

			prm.declare_entry("part_y",
				"2",
				Patterns::Integer(1,1000),
				"parts in y");

			prm.declare_entry("initial_refinement",
				"1",
				Patterns::Integer(1,50),
				"initial level of refinement");

			prm.declare_entry("inner_radius",
							 "0.5",
							 Patterns::Double(0.5,5.0),
							 "inner radius for the ring");

			prm.declare_entry("outer_radius",
							  "2.0",
							  Patterns::Double(0.1,10.0),
							  "outer radius of the ring");

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

			prm.declare_entry("print_velocity_space_error",
							  "false",
							  Patterns::Anything(),
							  "decide whether we want to print the error in the velocity space or not");

			prm.declare_entry("print_fe_index",
							  "false",
							  Patterns::Anything(),
							  "decide to print the fe index ");

		}
		prm.leave_subsection();

	}


	void Base_Constants
	::read_system_data()
	{
		bool entered = false;
		prm_system_info.enter_subsection("System Properties");
		{
			entered = true;

			constants_sys.total_systems = prm_system_info.get_integer("total_systems");
			

		}
		prm_system_info.leave_subsection();


			// we check the total number of systems we are loading into our program
			AssertThrow(constants_sys.total_systems == 3,
						ExcMessage("Incorrect number of systems being loaded"));


			constants_sys.nEqn.resize(constants_sys.total_systems);
			constants_sys.nBC.resize(constants_sys.total_systems);
			constants_sys.system_id.resize(constants_sys.total_systems);
			constants_sys.Ntensors.resize(constants_sys.total_systems);


		prm_system_info.enter_subsection("System Properties");
		{
			// details of the very first system
			if (constants_sys.total_systems == 1)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");				
			}

			if (constants_sys.total_systems == 2)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");				
			}

			if (constants_sys.total_systems == 3)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		
			}


			if (constants_sys.total_systems == 4)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		

				// details of the fourth system
				constants_sys.nEqn[3] = prm_system_info.get_integer("num_equations3");
				constants_sys.Ntensors[3] = prm_system_info.get_integer("Ntensors3");
				constants_sys.nBC[3] = prm_system_info.get_integer("nBC3");
				constants_sys.system_id[3] = prm_system_info.get_integer("system_id3");
			}

			if (constants_sys.total_systems == 5)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		

				// details of the fourth system
				constants_sys.nEqn[3] = prm_system_info.get_integer("num_equations3");
				constants_sys.Ntensors[3] = prm_system_info.get_integer("Ntensors3");
				constants_sys.nBC[3] = prm_system_info.get_integer("nBC3");
				constants_sys.system_id[3] = prm_system_info.get_integer("system_id3");

				constants_sys.nEqn[4] = prm_system_info.get_integer("num_equations4");
				constants_sys.Ntensors[4] = prm_system_info.get_integer("Ntensors4");
				constants_sys.nBC[4] = prm_system_info.get_integer("nBC4");
				constants_sys.system_id[4] = prm_system_info.get_integer("system_id4");

			}


			if (constants_sys.total_systems == 6)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		

				// details of the fourth system
				constants_sys.nEqn[3] = prm_system_info.get_integer("num_equations3");
				constants_sys.Ntensors[3] = prm_system_info.get_integer("Ntensors3");
				constants_sys.nBC[3] = prm_system_info.get_integer("nBC3");
				constants_sys.system_id[3] = prm_system_info.get_integer("system_id3");

				constants_sys.nEqn[4] = prm_system_info.get_integer("num_equations4");
				constants_sys.Ntensors[4] = prm_system_info.get_integer("Ntensors4");
				constants_sys.nBC[4] = prm_system_info.get_integer("nBC4");
				constants_sys.system_id[4] = prm_system_info.get_integer("system_id4");


				constants_sys.nEqn[5] = prm_system_info.get_integer("num_equations5");
				constants_sys.Ntensors[5] = prm_system_info.get_integer("Ntensors5");
				constants_sys.nBC[5] = prm_system_info.get_integer("nBC5");
				constants_sys.system_id[5] = prm_system_info.get_integer("system_id5");
			}


			if (constants_sys.total_systems == 7)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		

				// details of the fourth system
				constants_sys.nEqn[3] = prm_system_info.get_integer("num_equations3");
				constants_sys.Ntensors[3] = prm_system_info.get_integer("Ntensors3");
				constants_sys.nBC[3] = prm_system_info.get_integer("nBC3");
				constants_sys.system_id[3] = prm_system_info.get_integer("system_id3");

				constants_sys.nEqn[4] = prm_system_info.get_integer("num_equations4");
				constants_sys.Ntensors[4] = prm_system_info.get_integer("Ntensors4");
				constants_sys.nBC[4] = prm_system_info.get_integer("nBC4");
				constants_sys.system_id[4] = prm_system_info.get_integer("system_id4");


				constants_sys.nEqn[5] = prm_system_info.get_integer("num_equations5");
				constants_sys.Ntensors[5] = prm_system_info.get_integer("Ntensors5");
				constants_sys.nBC[5] = prm_system_info.get_integer("nBC5");
				constants_sys.system_id[5] = prm_system_info.get_integer("system_id5");

				constants_sys.nEqn[6] = prm_system_info.get_integer("num_equations6");
				constants_sys.Ntensors[6] = prm_system_info.get_integer("Ntensors6");
				constants_sys.nBC[6] = prm_system_info.get_integer("nBC6");
				constants_sys.system_id[6] = prm_system_info.get_integer("system_id6");
			}

			if (constants_sys.total_systems == 8)
			{
				constants_sys.nEqn[0] = prm_system_info.get_integer("num_equations0");
				constants_sys.Ntensors[0] = prm_system_info.get_integer("Ntensors0");
				constants_sys.nBC[0] = prm_system_info.get_integer("nBC0");
				constants_sys.system_id[0] = prm_system_info.get_integer("system_id0");			

				// details of the second system
				constants_sys.nEqn[1] = prm_system_info.get_integer("num_equations1");
				constants_sys.Ntensors[1] = prm_system_info.get_integer("Ntensors1");
				constants_sys.nBC[1] = prm_system_info.get_integer("nBC1");
				constants_sys.system_id[1] = prm_system_info.get_integer("system_id1");						


				// details of the second system
				constants_sys.nEqn[2] = prm_system_info.get_integer("num_equations2");
				constants_sys.Ntensors[2] = prm_system_info.get_integer("Ntensors2");
				constants_sys.nBC[2] = prm_system_info.get_integer("nBC2");
				constants_sys.system_id[2] = prm_system_info.get_integer("system_id2");		

				// details of the fourth system
				constants_sys.nEqn[3] = prm_system_info.get_integer("num_equations3");
				constants_sys.Ntensors[3] = prm_system_info.get_integer("Ntensors3");
				constants_sys.nBC[3] = prm_system_info.get_integer("nBC3");
				constants_sys.system_id[3] = prm_system_info.get_integer("system_id3");

				constants_sys.nEqn[4] = prm_system_info.get_integer("num_equations4");
				constants_sys.Ntensors[4] = prm_system_info.get_integer("Ntensors4");
				constants_sys.nBC[4] = prm_system_info.get_integer("nBC4");
				constants_sys.system_id[4] = prm_system_info.get_integer("system_id4");


				constants_sys.nEqn[5] = prm_system_info.get_integer("num_equations5");
				constants_sys.Ntensors[5] = prm_system_info.get_integer("Ntensors5");
				constants_sys.nBC[5] = prm_system_info.get_integer("nBC5");
				constants_sys.system_id[5] = prm_system_info.get_integer("system_id5");

				constants_sys.nEqn[6] = prm_system_info.get_integer("num_equations6");
				constants_sys.Ntensors[6] = prm_system_info.get_integer("Ntensors6");
				constants_sys.nBC[6] = prm_system_info.get_integer("nBC6");
				constants_sys.system_id[6] = prm_system_info.get_integer("system_id6");

				constants_sys.nEqn[7] = prm_system_info.get_integer("num_equations7");
				constants_sys.Ntensors[7] = prm_system_info.get_integer("Ntensors7");
				constants_sys.nBC[7] = prm_system_info.get_integer("nBC7");
				constants_sys.system_id[7] = prm_system_info.get_integer("system_id7");
			}
		}
		prm_system_info.leave_subsection();




		Assert(entered,ExcMessage("did not read system data"));
	
	}


	void Base_Constants
	::read_physical_constants()
	{
		bool entered = false;
		prm.enter_subsection("Physical Constants");
		{
			entered = true;
			constants_num.tau = prm.get_double("tau");
			constants_num.zeta = prm.get_double("zeta");
			constants_num.chi = prm.get_double("chi");
		
			constants_num.uW = prm.get_double("uW");
			constants_num.A0 = prm.get_double("A0");
		    constants_num.A1 = prm.get_double("A1");
			constants_num.A2 = prm.get_double("A2");
			constants_num.kappa  = prm.get_double("kappa");
			constants_num.epsilon = prm.get_double("epsilon");
			constants_num.error_variable = prm.get("error_variable");
			constants_num.force_variable = prm.get("force_variable");
			constants_num.alpha = prm.get_double("alpha");

			// densities of the incoming distributions
			constants_num.rho101 = prm.get_double("rho101");
			constants_num.rho102 = prm.get_double("rho102");

			constants_num.theta0 =  prm.get_double("theta0");
		 	constants_num.theta1  = prm.get_double("theta1");
		 	constants_num.theta2  = prm.get_double("theta2");
		 	constants_num.theta3  = prm.get_double("theta3");
		 	constants_num.theta4 = prm.get_double("theta4");

		 	// temperature of the incoming distributions
		 	constants_num.theta101 = prm.get_double("theta101");
		 	constants_num.theta102 = prm.get_double("theta102");

		 	constants_num.vx0 =  prm.get_double("vx0");
		 	constants_num.vx1  = prm.get_double("vx1");
		 	constants_num.vx2  = prm.get_double("vx2");
		 	constants_num.vx3  = prm.get_double("vx3");
		 	constants_num.vx4 = prm.get_double("vx4");

		 	// normal velocity of the incoming distributions
		 	constants_num.vx101 = prm.get_double("vx101");
		 	constants_num.vx102 = prm.get_double("vx102");


		 	constants_num.vy0 =  prm.get_double("vy0");
		 	constants_num.vy1  = prm.get_double("vy1");
		 	constants_num.vy2  = prm.get_double("vy2");
		 	constants_num.vy3  = prm.get_double("vy3");
		 	constants_num.vy4  = prm.get_double("vy4");

		 	// tangential velocity of the incoming distribution
		 	constants_num.vy101 = prm.get_double("vy101");
		 	constants_num.vy102 = prm.get_double("vy102");

		 	constants_num.coll_op = (Collision_Operator)prm.get_integer("Collision_Operator");

		 	constants_num.force_type = (Force_Type)prm.get_integer("force type");
		 	constants_num.bc_type = (BC_Type)prm.get_integer("BC type");
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
			constants_num.mesh_filename = prm.get("mesh_filename");
			constants_num.mesh_type = (Mesh_type)prm.get_integer("mesh type");
			constants_num.problem_type = (Problem_type)prm.get_integer("problem type");
			constants_num.initial_refinement = prm.get_integer("initial_refinement");
			constants_num.mesh_options = (Meshing_Options)prm.get_integer("mesh options");
			constants_num.xl = prm.get_double("left boundary");
			constants_num.xr = prm.get_double("right boundary");
			constants_num.yt = prm.get_double("top boundary");
			constants_num.yb = prm.get_double("bottom boundary");
			constants_num.part_x = prm.get_integer("part_x");
			constants_num.part_y = prm.get_integer("part_y");
			constants_num.inner_radius = prm.get_double("inner_radius");
			constants_num.outer_radius = prm.get_double("outer_radius");
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

			// input the output directory name, which comes from the 
			constants_num.main_output_dir = prm.get("Main Output Directory");
		}
		prm.leave_subsection();


		std::string temp = "output";

		for (int sys = 0 ; sys < constants_sys.total_systems ; sys++)
			temp +=  "_N" +  std::to_string(constants_sys.Ntensors[sys]);

		constants_num.main_output_dir = temp + "_" + constants_num.main_output_dir;

		Assert(entered,ExcMessage("Did not read the output directory name"));
	}

	void Base_Constants
	::allocate_variable_map()
	{

	  	
	  	std::vector<std::string> var_names = {"rho","vx","vy","theta","sigmaxx","sigmaxy","sigmayy","qx","qy","mxxx","mxxy","mxyy","myyy",
	  								"Delta","Rxx","Rxy","Ryy"};

	  	for (unsigned int i = 0 ; i < var_names.size() ; i++)
	  		constants_num.variable_map[var_names[i]] = i;
	  
	}

	void Base_Constants
	::allocate_variable_map_1D()
	{

	  	
	  	std::vector<std::string> var_names = {"rho","vx","theta","qx","mxxx","Delta","Rxx"};

	  	for (unsigned int i = 0 ; i < var_names.size() ; i++)
	  		constants_num.variable_map_1D[var_names[i]] = i;
	  
	}

	void Base_Constants
	::allocate_subdirector_names()
	{
		const unsigned int num_outputs = 5;
		constants_num.sub_directory_names.resize(num_outputs);

		constants_num.sub_directory_names[0] = constants_num.main_output_dir + "/grids";

		
		constants_num.sub_directory_names[1] = constants_num.main_output_dir + "/solution";

		
		constants_num.sub_directory_names[2] = constants_num.main_output_dir + "/convergence_tables";

		
		constants_num.sub_directory_names[3] = constants_num.main_output_dir + "/sparsity_patterns";

		
		constants_num.sub_directory_names[4] = constants_num.main_output_dir + "/matrix_read";			
		
	}

	void Base_Constants
	::read_printing_options()
	{
		bool entered = false;

		prm.enter_subsection("Printing Options");
		{
			constants_num.print_all = prm.get_bool("print_all");
			constants_num.print_solution = prm.get_bool("print_solution");
			constants_num.print_error = prm.get_bool("print_error");
			constants_num.print_exactsolution = prm.get_bool("print_exactsolution");
			constants_num.print_convergence_table = prm.get_bool("print_convergence_table");
			constants_num.print_velocity_space_error = prm.get_bool("print_velocity_space_error");
			constants_num.print_fe_index = prm.get_bool("print_fe_index");
            entered = true;
		}
		prm.leave_subsection();
        
        Assert(entered,ExcMessage("Did not read printing options."));
            
        
	}
}
