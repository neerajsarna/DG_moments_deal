// in the following namespace we test the Equation Generator with the help of A matrix
namespace Test_EquationGenerator
{
	using namespace dealii;

	void develop_A1(Full_matrix &A1, const unsigned int nEqn)
	{
		A1.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
	
		A1 << 0.,1.,0.,0.,0.,0.,1.,0.,0.,0.7886751345948129,0.,-0.21132486540518713,0.,0.,0.,0.,
   			  0.7071067811865476,0.,0.,0.7886751345948128,0.,0.,0.,0.,0.,0.,0.7071067811865476,0.,
   			  0.,0.,0.,-0.2113248654051871,0.,0.,0.,0.;

	}

	void develop_A2(Full_matrix &A2,const unsigned int nEqn)
	{
		A2.resize(nEqn,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
	
		A2 << 0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.7071067811865476,0.,1.,0.,0.,-0.21132486540518713,
   			  0.,0.7886751345948129,0.,0.,-0.2113248654051871,0.,0.,0.,0.,0.7071067811865476,0.,0.,
   			  0.,0.,0.,0.,0.7886751345948128,0.,0.,0.;


	}

	void develop_B(Sparse_matrix &B,const unsigned int nEqn, const unsigned int nBC)
	{
		B.resize(nBC,nEqn);
		Assert(nEqn == 6, ExcNotImplemented());
		switch(nEqn)
		{
			case 6:
			{
				Assert(B.rows() !=0,ExcNotInitialized());
				B.coeffRef(0, 0) =  -1.0000000000000;
				B.coeffRef(0, 1) =  1.0000000000000;
				B.coeffRef(0, 3) =  -1.0000000000000;
				B.coeffRef(1, 2) =  -1.0000000000000;
				B.coeffRef(1, 4) =   1.0000000000000;

				B.makeCompressed();
				break;
			}
			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
			}
		}
	}

	void develop_S_half(Full_matrix &S_half,const unsigned int nEqn)
	{
		S_half.resize(nEqn,nEqn);

		S_half  << 1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.3660254037844386,
   					0.,0.3660254037844386,0.,0.,0.,0.,1.4142135623730951,0.,0.,0.,0.,0.3660254037844386,
   					0.,1.3660254037844386;
	}

	void develop_S_half_inv(Full_matrix &S_half_inv,const unsigned int nEqn)
	{
		S_half_inv.resize(nEqn,nEqn);
		
		S_half_inv << 1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.7886751345948129,
   0.,-0.21132486540518713,0.,0.,0.,0.,0.7071067811865476,0.,0.,0.,0.,
   -0.21132486540518713,0.,0.7886751345948129;
	}


	void develop_P(Full_matrix &P,const unsigned int nEqn,const double tau)
	{
		P.resize(nEqn,nEqn);
		
		P << 0.,0.,0.,0.,0.,0.,0.,1.0/tau,0.,0.,0.,0.,0.,0.,1./tau,0.,0.,0.,0.,0.,0.,1./tau,0.,
   			-1.1102230246251565e-15,0.,0.,0.,0.,1.000000000000002/tau,0.,0.,0.,0.,
   			-1.1102230246251565e-15,0.,1./tau;
	}

	// develops the ID of odd variables
	void develop_odd_ID(MatrixUI &vector,const unsigned int nEqn)
	{
		switch(nEqn)
		{
			case 6:
			{
				vector.resize(2,1);
				vector(0,0) = 1;
				vector(1,0) = 4;
				break;
			}
			default:
			{
				ASSERT_EQ(1,0)<<"Should not have reached here";
				break;
			}
		}
	}

	TEST(DevelopingSystem,HandlingDevelopingSystems)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		Assert(constants.constants.tau<10.0,ExcNotImplemented());

		// no we will manually enter the matrices so as to compare with the read values
		// we will now input all the matrices manually so as to compare with the results
		Full_matrix A1;
		Full_matrix A2;
		Sparse_matrix B;
		Full_matrix S_half;
		Full_matrix S_half_inv;
		Full_matrix P;
		MatrixUI odd_ID;
		

		// develop A1
		develop_A1(A1,systemA.constants.nEqn);
		develop_A2(A2,systemA.constants.nEqn);
		develop_B(B,systemA.constants.nEqn,systemA.constants.nBC);
		develop_S_half(S_half,systemA.constants.nEqn);
		develop_S_half_inv(S_half_inv,systemA.constants.nEqn);
		develop_P(P,systemA.constants.nEqn,systemA.constants.tau);
		develop_odd_ID(odd_ID,systemA.constants.nEqn);

		
		Compare_Float_Mat(systemA.system_data.A[0].matrix,A1);
		Compare_Float_Mat(systemA.system_data.A[1].matrix,A2);
		Compare_Float_Mat(systemA.base_tensorinfo.S_half,S_half);
		Compare_Float_Mat(systemA.base_tensorinfo.S_half_inv,S_half_inv);
		Compare_Float_Mat(systemA.system_data.P.matrix,P);
		Compare_Float_Mat(systemA.system_data.B.matrix,B);
	
		for (unsigned int i = 0 ; i < 2 ; i ++)
			EXPECT_EQ(odd_ID(i,0),systemA.system_data.odd_ID(i,0));
	}
}
