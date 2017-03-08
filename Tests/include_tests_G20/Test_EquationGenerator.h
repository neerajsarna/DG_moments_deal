// in the following namespace we test the Equation Generator with the help of A matrix
namespace Test_EquationGenerator
{
	using namespace dealii;

	void develop_A1(Full_matrix &A1, const unsigned int nEqn)
	{
		A1.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
	
		A1 << 0.,1.,0.,0.,0.,0.,1.,0.,1.4142135623730951,0.,0.,0.,0.,1.4142135623730951,0.,
   			  1.7320508075688772,0.,0.,0.,0.,1.7320508075688772,0.,2.,0.,0.,0.,0.,2.,0.,
   			 2.23606797749979,0.,0.,0.,0.,2.23606797749979,0.;

	}

	void develop_B(Full_matrix &B,const unsigned int nEqn, const unsigned int nBC)
	{
		B.resize(nBC,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		
		B << -1e-5,1.,-1.4142135623730951 * 1e-5,0.,0.,0.,0.,0.,-0.9213177319235614,1.,-0.5319230405352439,
   			0.,0.,0.,0.41202581549140227,0.,-0.7136496464611086,1.;
	}

	void develop_S_half(Full_matrix &S_half,const unsigned int nEqn)
	{
		S_half.resize(nEqn,nEqn);
		S_half.setZero();

		for (unsigned int i = 0 ; i < nEqn ; i ++)
			S_half.coeffRef(i,i) = 1;
	}

	void develop_S_half_inv(Full_matrix &S_half_inv,const unsigned int nEqn)
	{
		S_half_inv.resize(nEqn,nEqn);
		S_half_inv.setZero();

		for (unsigned int i = 0 ; i < nEqn ; i ++)
			S_half_inv.coeffRef(i,i) = 1;

	}


	void develop_P(Full_matrix &P,const unsigned int nEqn,const double tau)
	{
		P.resize(nEqn,nEqn);
		P.setZero();
		

		for (unsigned int i = 3 ; i < nEqn ; i++)
			P.coeffRef(i,i) =1.0 / tau;
	}

	void develop_Binflow(Full_matrix &Binflow,const unsigned int nEqn, const unsigned int nBC)
	{
		Binflow.resize(nBC,nEqn);

		Binflow << -0.7978845608028654,1.,-0.5641895835477563,0.,0.16286750396763977,0.,
   					0.32573500793527993,0.,-0.6909882989426709,1.,-0.5984134206021489,0.,-0.2185096861184158,
   					0.,0.2575161346821262,0.,-0.6690465435572894,1.;

	}

	TEST(DevelopingSystem,HandlingDevelopingSystems)
	{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);

		G20::G20<dim> G20(constants.constants,folder_name);

		Assert(constants.constants.tau<10.0,ExcNotImplemented());
		Assert(dim == 1,ExcNotImplemented());


		fflush(stdout);

		// no we will manually enter the matrices so as to compare with the read values
		// we will now input all the matrices manually so as to compare with the results
		Full_matrix A1;
		Full_matrix B;
		Full_matrix S_half;
		Full_matrix S_half_inv;
		Full_matrix Binflow;
		Full_matrix P;
		MatrixUI odd_ID;
		
		S_half.setZero();
		S_half_inv.setZero();
		P.setZero();

		// develop A1
		develop_A1(A1,G20.constants.nEqn); 
		develop_B(B,G20.constants.nEqn,G20.constants.nBC);
		develop_S_half(S_half,G20.constants.nEqn);
		develop_S_half_inv(S_half_inv,G20.constants.nEqn);
		develop_P(P,G20.constants.nEqn,G20.constants.tau);
		develop_Binflow(Binflow,constants.constants.nEqn,constants.constants.nBC);

		
		Compare_Float_Mat(G20.system_data.A[0].matrix,A1);
		
		
		Compare_Float_Mat(G20.base_tensorinfo.S_half,S_half);
		Compare_Float_Mat(G20.base_tensorinfo.S_half_inv,S_half_inv);
		Compare_Float_Mat(G20.system_data.P.matrix,P);
		Compare_Float_Mat(G20.system_data.B.matrix,B);
		Compare_Float_Mat(G20.system_data.Binflow.matrix,Binflow);

		for (unsigned int i = 0 ; i < constants.constants.nBC ; i++)
			std::cout << "Odd ID " << G20.system_data.odd_ID(i,0) << std::endl;
	
	}
}
