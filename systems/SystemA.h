// In the following file we declare all the attributes of the A system
namespace SystemA
{
	using namespace dealii;


	// In the following class we will construct the system and all of it's attributes.
	template<int dim>
	class
	Base_SystemA:public EquationGenerator::Base_EquationGenerator<dim>
	{
		public:
			Base_SystemA(const constant_data &constants,
						 std::string &folder_name);

			const constant_data constants;
			typename EquationGenerator::Base_EquationGenerator<dim>::equation_data systemA_data;	
			std::vector<std::string> basefile_systemA;
			virtual void build_P(Sparse_matrix &P);
	};

	// folder name is the name of the folder in which the system_matrices have been kept
	template<int dim>
	Base_SystemA<dim>
	::Base_SystemA(const constant_data &constants,
				   std::string &folder_name)
	:
	// we need to first initialize the constant structure
	constants(constants)
	{
		// we generate all the matrices for this system
		this->generate_matrices(systemA_data,constants.nEqn,constants.nBC,folder_name);

		// we now develop the P matrix
		build_P(systemA_data.P.matrix);

	}

	template<int dim>
	void 
	Base_SystemA<dim>
	::build_P(Sparse_matrix &P)
	{
		// the first one is the conservation law
		for (unsigned int i = 1 ; i < P.rows() ; i ++)
			P.coeffRef(i,i) = 1/constants.tau;
	}

	// accomodation of different test case.
	// the test cases vary in the sense of the force type
	// and the other variation is the boundary conditions
	template<int dim>
	class
	SystemA:public Base_SystemA<dim>
	{
		public:
			SystemA(const constant_data &constants,
					std::string &folder_name);

			
			// we have a pointer to the kind of force we want in the system
			ForceType::Base_ForceType<dim> *force;

			// we will use them to initialize force 
			ForceType::ForceType1<dim> force1;
			ForceType::ForceType2<dim> force2;
			ForceType::ForceType3<dim> force3;

			void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value);

			// matrices for the boundary conditions
			Full_matrix B_tilde_inv;
			Full_matrix B_hat;

	};

	template<int dim>
	SystemA<dim>
	::SystemA(const constant_data &constants,
			  std::string &folder_name)
	:
	Base_SystemA<dim>(constants,folder_name),
	force1(constants),
	force2(constants),
	force3(constants)
	{
		// now we switch depending upon the type of force
		switch(this->constants.force_type)
		{
			case type1:
			{
				force = &force1;
				break;
			}

			case type2:
			{
				force = &force2;
				break;
			}

			case type3:
			{
				force = &force3;
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
			case characteristic:
			{
				BoundaryHandler::Base_BoundaryHandler_Char<dim> boundary_handler_char(this->systemA_data.Ax.matrix,
																				this->systemA_data.B.matrix,
																				this->constants.nBC);

				// decision whether to change the boundary conditions for vx or not
				bool action = true;
				boundary_handler_char.fix_B_vx(this->constants.epsilon,action);

				B_tilde_inv.resize(this->constants.nBC,this->constants.nBC);
				B_hat.resize(this->constants.nEqn,this->constants.nBC);

				B_tilde_inv = boundary_handler_char.build_B_tilde_inv();
				B_hat = boundary_handler_char.build_B_hat(B_tilde_inv);

				// we need to check whether B has been fixed or not
				if (action)
					Assert(fabs(this->systemA_data.B.matrix.coeffRef(0,0)) < 1e-3,ExcMessage("B not changed for vx"));

				// we compress B again
				this->systemA_data.B.matrix.makeCompressed();
				break;
			}

			case odd:
			{
				Assert(1 == 0, ExcNotImplemented());
				break;
			}

			default:
			{
				Assert(1 == 0,ExcNotImplemented());
				break;
			}
		}
	}

	template<int dim>
	void 
	SystemA<dim>
	::source_term(const std::vector<Point<dim>> &p,
				 std::vector<Vector<double>> &value)
	{
		// now we simply pass on the value to the base class
		force->source_term(p,value);
	}

}


