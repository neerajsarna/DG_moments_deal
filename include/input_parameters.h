namespace Input_parameters
{
	using namespace dealii;
	using namespace std;

	class parameter_handling
	{
		public:

			ParameterHandler prm;
			parameter_handling();
			void read_input_file(const string &input_file);
			void read_numerical_constants(numerical_data &numerical_constants);
			void read_system_data(nEqn_data &num_equations);
			void read_physical_constants(physical_data &physical_constants);
			void read_output_dir_name(string &output_dir);
			void read_mesh_info(mesh_data &mesh_info);

		private:
			void declare_parameters();
						
	};

	parameter_handling
	::parameter_handling()
	{
		declare_parameters();
	}

	void parameter_handling
	::declare_parameters()
	{
			prm.enter_subsection("Numerical Constants");
			{
				prm.declare_entry("polynomial degree",
								  "1",
								  Patterns::Integer(1,10),
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
			
				prm.declare_entry("theta0",
								  "2.0",
								  Patterns::Double(0.001,20.0),
								  "Temperature of inner wall");

		
				prm.declare_entry("theta1",
								  "1.0",
								  Patterns::Double(0.001,20.0),
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
				prm.declare_entry("total systems to load",
								   "1",
								   Patterns::Integer(1,2),
								   "total number of systems to load");

				prm.declare_entry("equations in the system",
					 			  "6",
					 			  Patterns::Integer(6,10),
					 			  "total number of equations in system");

				prm.declare_entry("nBC",
								  "2",
								  Patterns::Integer(2,10),
								  "total number of boundary conditions needed for the system");

				/*by default we will solve the symmetric system*/
				prm.declare_entry("symmetric or un symmetric",
								  "0",
								  Patterns::Integer(0,1),
								  "to solve the symmetric system or the un symmetric system");

				/*0 for force of system A and 1 for force of system B*/
				prm.declare_entry("force type",
								  "0",
								  Patterns::Integer(0,1),
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
								   Patterns::Double(-10.0,-0.5),
								   "left boundary");

				prm.declare_entry("right boundary",
								  "1",
								  Patterns::Double(0.5,10.0),
								   "right boundary");

				prm.declare_entry("bottom boundary",
								  "-1",
								  Patterns::Double(-10.0,-0.5),
								   "bottom boundary of square domain");

				prm.declare_entry("top boundary",
									"1",
									Patterns::Double(0.5,10.0),
									"top boundary of square domain");				
			}
			prm.leave_subsection();

	}

	void parameter_handling
	::read_input_file(const string &input_file)
	{
		prm.read_input(input_file.c_str());
	}

	void parameter_handling
	::read_numerical_constants(numerical_data &numerical_constants)
	{
		bool entered = false;
		prm.enter_subsection("Numerical Constants");
		{
			entered = true;
			numerical_constants.p = prm.get_integer("polynomial degree");
			numerical_constants.mapping_order = prm.get_integer("mapping order");
			numerical_constants.refine_cycles = prm.get_integer("total h refinement cycles");
			numerical_constants.refinement = (Refinement)prm.get_integer("type of refinement");
			numerical_constants.assembly_type = (Assembly_Type)prm.get_integer("assembly type");
		}
		prm.leave_subsection();

		Assert(entered, ExcMessage("did not read numerical constants"));
	}

	void parameter_handling
	::read_system_data(nEqn_data &num_equations)
	{
		bool entered = false;
		prm.enter_subsection("System Properties");
		{
			entered = true;

			num_equations.no_of_total_systems = prm.get_integer("total systems to load");
			Assert(num_equations.no_of_total_systems == 1,ExcNotImplemented());

			num_equations.total_nEqn.resize(num_equations.no_of_total_systems);
			num_equations.nBC.resize(num_equations.no_of_total_systems);

			num_equations.total_nEqn[0] = prm.get_integer("equations in the system");
			num_equations.nBC[0] = prm.get_integer("nBC");
			num_equations.system_to_solve = 0;

			num_equations.system_type = (System_Type)prm.get_integer("symmetric or un symmetric");
			num_equations.force_type = (Force_Type)prm.get_integer("force type");
			num_equations.bc_type = (BC_type)prm.get_integer("BC type");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read system data"));
	}

	void parameter_handling
	::read_physical_constants(physical_data &physical_constants)
	{
		bool entered = false;
		prm.enter_subsection("Physical Constants");
		{
			entered = true;
			physical_constants.tau = prm.get_double("tau");

		 	physical_constants.zeta = prm.get_double("zeta");
			physical_constants.chi  = prm.get_double("chi");
			physical_constants.theta0 =  prm.get_double("theta0");
		 	physical_constants.theta1  = prm.get_double("theta1");
			physical_constants.uW = prm.get_double("uW");
			physical_constants.A0 = prm.get_double("A0");
		    physical_constants.A1 = prm.get_double("A1");
			physical_constants.A2 = prm.get_double("A2");
			physical_constants.kappa  = prm.get_double("kappa");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read physical constants"));
	}

	void parameter_handling
	::read_output_dir_name(string &output_dir)
	{
		bool entered = false;

		prm.enter_subsection("Output Directory Name");
		{
			entered = true;
			output_dir = prm.get("Main Output Directory");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read output dir name"));
	}

	void parameter_handling
	::read_mesh_info(mesh_data &mesh_info)
	{
		bool entered = false;

		prm.enter_subsection("Mesh Info");
		{
			entered = true;
			mesh_info.mesh_filename = prm.get("filename");
			mesh_info.mesh_type = (Mesh_type)prm.get_integer("mesh type");
			mesh_info.mesh_options = (Meshing_Options)prm.get_integer("mesh options");
			mesh_info.xl = prm.get_double("left boundary");
			mesh_info.xr = prm.get_double("right boundary");
			mesh_info.yt = prm.get_double("top boundary");
			mesh_info.yb = prm.get_double("bottom boundary");
		}
		prm.leave_subsection();

		Assert(entered,ExcMessage("did not read mesh info"));
	}
}