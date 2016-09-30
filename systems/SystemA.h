// In the following file we declare all the attributes of the A system
namespace SystemA
{
	using namespace dealii;


	// In the following class we will construct all the matrices associated with the systemA
	template<int dim>
	class
	Base_SystemA:public EquationGenerator::Base_EquationGenerator<dim>
	{
		public:
			Base_SystemA(const constant_data &constants,
						std::string &folder_name);
			
			virtual void build_P(Sparse_matrix &P);

			// special function to create varIdx for this particular system
			MatrixUI build_varIdx();

			// prepare varIdx for systemA, since we do not have a one to one correspondence with 
			const unsigned int Ntensors = 3;
	};

	// folder name is the name of the folder in which the system_matrices have been kept
	template<int dim>
	Base_SystemA<dim>
	::Base_SystemA(const constant_data &constants,
				   std::string &folder_name)
	:
	// we need to first initialize the constant structure
	EquationGenerator::Base_EquationGenerator<dim>(constants)
	{
		// we generate all the matrices for this system
		this->generate_matrices(this->system_data,this->constants.nEqn,this->constants.nBC,folder_name);

		// we now develop the P matrix
		build_P(this->system_data.P.matrix);

		this->varIdx = build_varIdx();

	}

	template<int dim>
	void 
	Base_SystemA<dim>
	::build_P(Sparse_matrix &P)
	{
		// the first one is the conservation law
		for (unsigned int i = 1 ; i < P.rows() ; i ++)
			P.coeffRef(i,i) = 1/this->constants.tau;
	}

	template<int dim>
	MatrixUI
	Base_SystemA<dim>
	::build_varIdx()
	{
		MatrixUI result;
		result.resize(Ntensors,2);

		Assert(result.rows() !=0 ,ExcNotInitialized());
		Assert(result.cols() !=0 ,ExcNotInitialized());
		
		// theta
		result(0,0) = 1;
		result(0,1) = 0;

		// q
		result(1,0) = 0;
		result(1,1) = 1;

		// R
		result(2,0) = 1;
		result(2,1) = 2;

		return(result);
	}

	// accomodation of different test case.
	// the test cases vary in the sense of the force type
	// and the other variation is the boundary conditions
	template<int dim>
	class
	SystemA:public Base_SystemA<dim>
	{
		public:
			// we need the input parameters to initialize the base class
			SystemA(const constant_data &constants,
					std::string &folder_name);

			// routines for tensor development
			TensorInfo::Base_TensorInfo<dim> base_tensorinfo;

			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value);

			virtual Sparse_matrix build_Projector(const Tensor<1,dim> &normal_vector);
			virtual Sparse_matrix build_InvProjector(const Tensor<1,dim> &normal_vector);

			virtual void build_BCrhs(const Tensor<1,dim,double> p,
									const Tensor<1,dim,double> normal_vector,
									Vector<double> &bc_rhs);

			virtual Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector);
	};

	template<int dim>
	SystemA<dim>
	::SystemA(const constant_data &constants,
			  std::string &folder_name)
	:
	Base_SystemA<dim>(constants,folder_name),
	base_tensorinfo(this->varIdx,this->Ntensors)
	{


		// now we switch depending upon the type of force
		switch(this->constants.force_type)
		{
			case type1:
			{
				this->force = &this->force1;
				break;
			}

			case type2:
			{
				this->force = &this->force2;
				break;
			}

			case type3:
			{
				this->force = &this->force3;
				break;
			}

			default:
			{
				AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
				break;

			}
		}

		// develop the matrices which are independent of the test case
		switch(this->constants.bc_type)
		{
			// characteristic boundary conditions
			case characteristic:
			{
				BoundaryHandler::Base_BoundaryHandler_Char<dim> boundary_handler_char(this->system_data.Ax.matrix,
																				this->system_data.B.matrix,
																				this->constants.nBC);

				// we do not need to fix B for the present problem since we do not have a velocity in the system

				this->B_tilde_inv.resize(this->constants.nBC,this->constants.nBC);
				this->B_hat.resize(this->constants.nEqn,this->constants.nEqn);

				this->B_tilde_inv = boundary_handler_char.build_B_tilde_inv();
				this->B_hat = boundary_handler_char.build_B_hat(this->B_tilde_inv);

				break;
			}

			// picking up the odd variables
			case odd:
			{
				BoundaryHandler::Base_BoundaryHandler_Odd<dim> boundary_handler_odd(this->system_data.B.matrix,
																				   this->system_data.odd_ID);

				 this->BC = boundary_handler_odd.develop_BC();
				break;
			}

			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
				break;
			}
		}

		//the following implementation is mesh dependent
		switch(this->constants.mesh_type)
		{
			case ring:
			{
				switch (this->constants.bc_type)
				{
					case characteristic:
					{
						this->base_bcrhs = &this->bcrhs_ring_char;
						break;
					}
					case odd:
					{
						this->base_bcrhs = &this->bcrhs_ring_odd;
						break;
					}
				}
				break;
			}

			case periodic_square:
			{
				switch(this->constants.bc_type)
				{
					case characteristic:
					{
						this->base_bcrhs = &this->bcrhs_periodic_char;
						break;
					}

					case odd:
					{
						this->base_bcrhs = &this->bcrhs_periodic_odd;
						break;
					}

					default:
					{
						Assert(1 == 0, ExcMessage("Should not have reached here"));
						break;
					}
				}

				break;
			}
			default:
			{
				Assert(1 ==0,ExcNotImplemented());
				break;
			}
		}

		// first we create the class which handles matrix orperations
		MatrixOpt::Base_MatrixOpt matrixopt;

		// we will always send Ax independent of symmetric or unsymmetric system
		this->Aminus_1D = matrixopt.compute_Aminus(this->system_data.Ax.matrix);

		this->X_minus = matrixopt.compute_Xminus(this->system_data.Ax.matrix,this->constants.nBC);

	}
	

	template<int dim>
	void 
	SystemA<dim>
	::source_term(const std::vector<Point<dim>> &p,
				 std::vector<Vector<double>> &value)
	{
		// now we simply pass on the value to the base class
		this->force->source_term(p,value);
	}

	// develops the normal vector corresponding to a normal vector
	template<int dim>
	Sparse_matrix
	SystemA<dim>
	::build_Projector(const Tensor<1,dim> &normal_vector)
	{

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		this->base_tensorinfo.reinit_global(nx,ny);

		return(base_tensorinfo.global_Projector);
	}


	template<int dim>
	Sparse_matrix
	SystemA<dim>
	::build_InvProjector(const Tensor<1,dim> &normal_vector)
	{

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		this->base_tensorinfo.reinit_Invglobal(nx,ny);

		return(this->base_tensorinfo.global_Projector);
	}


	template<int dim>
	void
	SystemA<dim>
	::build_BCrhs(const Tensor<1,dim,double> p,
				const Tensor<1,dim,double> normal_vector,
				Vector<double> &bc_rhs)
	{

		this->base_bcrhs->BCrhs(p,normal_vector,bc_rhs);
	}

	// computation of Aminus for An
	template<int dim>
	Full_matrix
	SystemA<dim>
	::build_Aminus(const Tensor<1,dim,double> normal_vector)
	{
		Full_matrix Aminus;
		Aminus.resize(this->constants.nEqn,this->constants.nEqn);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		Assert(this->Aminus_1D.rows() != 0 || this->Aminus_1D.cols() != 0,
			  ExcNotInitialized());

		// first multiply from the left
		base_tensorinfo.reinit_Invglobal(nx,ny);
		Aminus = base_tensorinfo.global_Projector * this->Aminus_1D;

		// now multiply from the right
		base_tensorinfo.reinit_global(nx,ny);
		Aminus = Aminus * base_tensorinfo.global_Projector;

		return(Aminus);
	}
	

}


