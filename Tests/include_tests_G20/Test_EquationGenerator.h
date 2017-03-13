// in the following namespace we test the Equation Generator with the help of A matrix
namespace Test_EquationGenerator
{
	using namespace dealii;

	void develop_A1(Full_matrix &A1, const unsigned int nEqn)
	{
		A1.resize(nEqn,nEqn);
		Assert(nEqn == 6 || nEqn == 8 || nEqn == 9 || nEqn == 11, ExcNotImplemented());
	
		if (nEqn == 6)
			A1 << 0.,1.,0.,0.,0.,0.,1.,0.,1.4142135623730951,0.,0.,0.,0.,1.4142135623730951,0.,
   			  1.7320508075688772,0.,0.,0.,0.,1.7320508075688772,0.,2.,0.,0.,0.,0.,2.,0.,
   			 2.23606797749979,0.,0.,0.,0.,2.23606797749979,0.;

   		if (nEqn == 8)
   			A1 << 0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,1.4142135623730951,0.,0.,0.,0.,0.,0.,
   1.4142135623730951,0.,1.7320508075688772,0.,0.,0.,0.,0.,0.,1.7320508075688772,0.,2.,
   0.,0.,0.,0.,0.,0.,2.,0.,2.23606797749979,0.,0.,0.,0.,0.,0.,2.23606797749979,0.,
   2.449489742783178,0.,0.,0.,0.,0.,0.,2.449489742783178,0.,2.6457513110645907,0.,0.,0.,
   0.,0.,0.,2.6457513110645907,0.;

   		if (nEqn == 9)
   			A1 << 0.,1.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.4142135623730951,0.,0.,0.,0.,0.,0.,0.,
   1.4142135623730951,0.,1.7320508075688772,0.,0.,0.,0.,0.,0.,0.,1.7320508075688772,0.,
   2.,0.,0.,0.,0.,0.,0.,0.,2.,0.,2.23606797749979,0.,0.,0.,0.,0.,0.,0.,2.23606797749979,
   0.,2.449489742783178,0.,0.,0.,0.,0.,0.,0.,2.449489742783178,0.,2.6457513110645907,0.,
   0.,0.,0.,0.,0.,0.,2.6457513110645907,0.,2.8284271247461903,0.,0.,0.,0.,0.,0.,0.,
   2.8284271247461903,0.;


   		if (nEqn == 11)
   			A1 << 0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.4142135623730951,0.,0.,0.,0.,0.,0.,0.,0.,
   0.,1.4142135623730951,0.,1.7320508075688772,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   1.7320508075688772,0.,2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.,0.,2.23606797749979,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,2.23606797749979,0.,2.449489742783178,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   2.449489742783178,0.,2.6457513110645907,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   2.6457513110645907,0.,2.8284271247461903,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   2.8284271247461903,0.,3.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3.,0.,3.1622776601683795,0.,0.,
   0.,0.,0.,0.,0.,0.,0.,3.1622776601683795,0.;


	}

	void develop_B(Full_matrix &B,const unsigned int nEqn, const unsigned int nBC)
	{
		B.resize(nBC,nEqn);
		Assert(nEqn == 6 || nEqn == 8 || nEqn == 9 || nEqn == 11, ExcNotImplemented());

		if (nEqn == 6)
			B << -1e-5,1.,-1.4142135623730951 * 1e-5,0.,0.,0.,0.,0.,-0.9213177319235614,1.,-0.5319230405352439,
   				0.,0.,0.,0.41202581549140227,0.,-0.7136496464611086,1.;

   		if(nEqn == 8)
   			B << -1e-5,1.,-1.4142135623730951 * 1e-5,0.,0.,0.,0.,0.,0.,0.,-0.9213177319235614,1.,
   -0.5319230405352436,0.,0.14567312407894378,0.,0.,0.,0.41202581549140227,0.,
   -0.7136496464611084,1.,-0.5863230142835033,0.,0.,0.,-0.2860963362053377,0.,
   0.27529632787052893,0.,-0.6785370394786211,1.;

   		if(nEqn == 9)
   			B << -1e-5,1.,-1.4142135623730951 * 1e-5,0.,0.,0.,0.,0.,0.,0.,0.,-0.9213177319235614,1.,
   -0.5319230405352436,0.,0.14567312407894378,0.,-0.4671933606552149,0.,0.,
   0.41202581549140227,0.,-0.7136496464611084,1.,-0.5863230142835033,0.,
   0.7660958162452001,0.,0.,-0.2860963362053377,0.,0.27529632787052893,0.,
   -0.6785370394786211,1.,-1.434653512450009;


   		if (nEqn == 11)
   			B << -1e-5,1.,-1.4142135623730951 * 1e-5,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.9213177319235614,1.,
   -0.5319230405352436,0.,0.14567312407894378,0.,-0.07786556010920254,0.,
   0.4103875353830485,0.,0.,0.41202581549140227,0.,-0.7136496464611084,1.,
   -0.5863230142835033,0.,0.17411268551027304,0.,-0.6240050098398657,0.,0.,
   -0.2860963362053377,0.,0.27529632787052893,0.,-0.6785370394786211,1.,
   -0.6044888395154532,0.,0.8750703998273118,0.,0.,0.22477851044824515,0.,
   -0.18168630692147242,0.,0.24878422172674358,0.,-0.6649038006690545,1.,
   -1.5258944299836208;

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

		if (nEqn == 6)
			Binflow << -0.7978845608028654,1.,-0.5641895835477563,0.,0.16286750396763977,0.,
   						0.32573500793527993,0.,-0.6909882989426709,1.,-0.5984134206021489,0.,-0.2185096861184158,
   						0.,0.2575161346821262,0.,-0.6690465435572894,1.;

   		if (nEqn == 8)
   			Binflow << -0.7978845608028654,1.,-0.5641895835477563,0.,0.16286750396763977,0.,
   						-0.08920620580763855,0.,0.32573500793527993,0.,-0.6909882989426709,1.,
   						-0.5984134206021489,0.,0.18209140509867972,0.,-0.2185096861184158,0.,
   						0.2575161346821262,0.,-0.6690465435572894,1.,-0.6107531398786498,0.,
   						0.16858388283618386,0.,-0.16688952945311364,0.,0.24088428688671293,0.,
   						-0.6596887883819926,1.;

   		if(nEqn == 9)
   			Binflow << -0.7978845608028654,1.,-0.5641895835477563,0.,0.16286750396763977,0.,
   						-0.08920620580763855,0.,0.47682722700889624,0.32573500793527993,0.,
   						-0.6909882989426709,1.,-0.5984134206021489,0.,0.18209140509867972,0.,
   						-0.6618572609282211,-0.2185096861184158,0.,0.2575161346821262,0.,-0.6690465435572894,
   						1.,-0.6107531398786498,0.,0.8966803303779044,0.16858388283618386,0.,
   						-0.16688952945311364,0.,0.24088428688671293,0.,-0.6596887883819926,1.,
   						-1.5354016523692513;


   		if (nEqn == 11)
   			Binflow << -0.7978845608028654,1.,-0.5641895835477563,0.,0.16286750396763977,0.,
   -0.08920620580763855,0.,0.05960340337611214,0.,-0.43979252558799503,
   0.32573500793527993,0.,-0.6909882989426709,1.,-0.5984134206021489,0.,
   0.18209140509867972,0.,-0.10219854764332836,0.,0.5899320821131322,
   -0.2185096861184158,0.,0.2575161346821262,0.,-0.6690465435572894,1.,
   -0.6107531398786498,0.,0.1904357497768605,0.,-0.7444471532832225,0.16858388283618386,
   0.,-0.16688952945311364,0.,0.24088428688671293,0.,-0.6596887883819926,1.,
   -0.6170823570053586,0.,0.967993530876935,-0.13907460787759468,0.,0.1264379121271379,
   0.,-0.15329782146499238,0.,0.23323520786882213,0.,-0.6545146787836008,1.,
   -1.6025521020936084;

	}

	TEST(DevelopingSystem,HandlingDevelopingSystems)
	{
		const unsigned int dim = 1;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);

	// 	// we construct a vector of all the system data being considered 
		std::vector<Develop_System::System<dim>> System;

		// initialize the vector containing all the systems
		for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
			 System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
			 											constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));

		Assert(constants.constants_num.tau<10.0,ExcNotImplemented());
		Assert(dim == 1,ExcNotImplemented());


		fflush(stdout);

		std::vector<Full_matrix> A1(constants.constants_sys.total_systems);
		std::vector<Full_matrix> B(constants.constants_sys.total_systems);
		std::vector<Full_matrix> S_half(constants.constants_sys.total_systems);
		std::vector<Full_matrix> S_half_inv(constants.constants_sys.total_systems);
		std::vector<Full_matrix> Binflow(constants.constants_sys.total_systems);
		std::vector<Full_matrix> P(constants.constants_sys.total_systems);
		std::vector<MatrixUI> odd_ID(constants.constants_sys.total_systems);

		// develop A1
		for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		{
			develop_A1(A1[i],constants.constants_sys.nEqn[i]); 
			develop_B(B[i],constants.constants_sys.nEqn[i],constants.constants_sys.nBC[i]);
			develop_S_half(S_half[i],constants.constants_sys.nEqn[i]);
			develop_S_half_inv(S_half_inv[i],constants.constants_sys.nEqn[i]);
			develop_P(P[i],constants.constants_sys.nEqn[i],constants.constants_num.tau);
			develop_Binflow(Binflow[i],constants.constants_sys.nEqn[i],constants.constants_sys.nBC[i]);	

			Compare_Float_Mat(System[i].system_data.A[0].matrix,A1[i]);
			Compare_Float_Mat(System[i].base_tensorinfo.S_half,S_half[i]);
			Compare_Float_Mat(System[i].base_tensorinfo.S_half_inv,S_half_inv[i]);
			Compare_Float_Mat(System[i].system_data.P.matrix,P[i]);
			Compare_Float_Mat(System[i].system_data.B.matrix,B[i]);
			Compare_Float_Mat(System[i].system_data.Binflow.matrix,Binflow[i]);

			for (unsigned int j = 0 ; j < constants.constants_sys.nBC[i] ; j++)
				std::cout << "Odd ID " << System[i].system_data.odd_ID(j,0) << std::endl;				
		}

	
	 }
}
