namespace Test_BaseConstants
{
	TEST(ReadingInputs,HandlingInputs)
	{
		std::string input_file = "../test_input_files/input1.in";
		Constants::Base_Constants base_constants(input_file);
		
		// first we allocate the values for the variables and then 
		// we check whether that is what we expect or not 
				// set of numerical constants
		int p = 1;
		int mapping_order = 1;
		int refine_cycles = 10;
		Refinement refinement = global;
		Assembly_Type assembly_type = meshworker;

		// data set for system info
		int nEqn = 6;
		int nBC = 2;
		int system_id = 6;

		System_Type system_type = un_symmetric;
		Force_Type force_type = type1;
		BC_Type bc_type = odd;

		// physical constants
		double tau = 0.1;
		double zeta = 1.0;
		double chi = 1.0;
		double theta0 = 2.0;
		double theta1 = 1.0;
		double uW = tau;
		double A0 = 0.0;
		double A1 = 0.2;
		double A2 = 0.1;
		double kappa = 0.0;
		
		// coefficient for normal relaxational velocity
		double epsilon = 1e-5;

		// variable in which error has to be found
		std::string error_variable = "theta";

		// forcing term for poisson heat conduction
		double alpha = 1.0;

		// the type of mesh to be used
		Mesh_type mesh_type = periodic;

		// different options for meshing
		Meshing_Options mesh_options = generate_internal;

		// name of the mesh file 
		std::string mesh_filename = "mesh";
		

	// to be only used for square domain
	// xcord of left boundary
		double xl = -1.0;

	// xcord of the right boundary
		double xr = 1.0;

	// ycord of the bottom boundary
		double yb = -1.0;

	// ycord of the top boundary
		double yt = 1.0;

	// number of partisions per dimension for the square domain
		unsigned int part_x = 10;

		unsigned int part_y = 10;

/***************NOW WE TEST THE VALUES WHICH HAVE BEEN READ ************/

		EXPECT_EQ(p,base_constants.constants.p);
		EXPECT_EQ(mapping_order,base_constants.constants.mapping_order);
		EXPECT_EQ(refine_cycles,base_constants.constants.refine_cycles);
		EXPECT_EQ(refinement,base_constants.constants.refinement);
		EXPECT_EQ(assembly_type,base_constants.constants.assembly_type);

		// data set for system info
		EXPECT_EQ(nEqn,base_constants.constants.nEqn);
		EXPECT_EQ(nBC,base_constants.constants.nBC);
		EXPECT_EQ(system_id,base_constants.constants.system_id);

		EXPECT_EQ(system_type,base_constants.constants.system_type);
		EXPECT_EQ(force_type,base_constants.constants.force_type);
		EXPECT_EQ(bc_type,base_constants.constants.bc_type);

		
		// physical constants
		EXPECT_EQ(tau,base_constants.constants.tau);
		EXPECT_EQ(zeta,base_constants.constants.zeta);
		EXPECT_EQ(chi,base_constants.constants.chi);
		EXPECT_EQ(theta0,base_constants.constants.theta0);
		EXPECT_EQ(theta1,base_constants.constants.theta1);
		EXPECT_EQ(uW,base_constants.constants.tau);
		EXPECT_EQ(A0,base_constants.constants.A0);
		EXPECT_EQ(A1,base_constants.constants.A1);
		EXPECT_EQ(A2,base_constants.constants.A2);
		EXPECT_EQ(kappa,base_constants.constants.kappa);
		
		// coefficient for normal relaxational velocity
		EXPECT_EQ(epsilon,base_constants.constants.epsilon);
		// variable in which error has to be found
		EXPECT_STREQ(error_variable.c_str(),base_constants.constants.error_variable.c_str());

		// forcing term for poisson heat conduction
		EXPECT_EQ(alpha,base_constants.constants.alpha);

		// the type of mesh to be used
		EXPECT_EQ(mesh_type,base_constants.constants.mesh_type);

		// different options for meshing
		EXPECT_EQ( mesh_options ,base_constants.constants.mesh_options);
		

	// to be only used for square domain
	// xcord of left boundary
		EXPECT_EQ(xl,base_constants.constants.xl);

	// xcord of the right boundary
		EXPECT_EQ(xr,base_constants.constants.xr);

	// ycord of the bottom boundary
		EXPECT_EQ(yb,base_constants.constants.yb);

	// ycord of the top boundary
		EXPECT_EQ(yt,base_constants.constants.yt);

	// number of partisions per dimension for the square domain
		EXPECT_EQ(part_x,base_constants.constants.part_x);

		EXPECT_EQ(part_y,base_constants.constants.part_y);
	}

}