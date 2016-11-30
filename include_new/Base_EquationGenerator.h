// the following is the base of the equation generator, it contains all the routines which
// which are equation dependent. This includes the data for the equations and the Projector matrices.
// It is an abstract base class with a specifiied routine for the reading of the matrices. 
namespace EquationGenerator
{
	using namespace dealii;

	template<int dim>
	class Base_EquationGenerator
	{
		public:
			Base_EquationGenerator(const constant_data &constants);
			
			const constant_data constants;
			struct equation_data
			{

			// matrix to store the jacobian of the flux
				 system_matrix A[dim];

		  	// The Production term
				  system_matrix P;

			// The Jacobian in the x-direction
				system_matrix Ax;

			// The matrix for the boundary
				system_matrix B;

			// ID of the odd variables 
				MatrixUI odd_ID;

			};

			// a data structure to store all the system data
			equation_data system_data;

			// the base filenames for the system matrices
			std::vector<std::string> basefile;

			// Given the filename, read the tripet for the sparse matrix 
			void build_triplet(std::vector<triplet> &Row_Col_Value,
								const std::string filename);

			// Develop the matrix from the triplet
	  		void build_matrix_from_triplet(Sparse_matrix &matrix,
	  									   std::vector<triplet> &Row_Col_Value);

	  		// reads a vector from a file
	  		void build_Vector(MatrixUI &vector,const std::string filename);

	  		// Print the given matrix to the desired file
	  		void print_matrix(Sparse_matrix &matrix,const std::string filename);


	  		// generate the matrix for a particular system
	  		void generate_matrices(equation_data &system_data,
	  							   const unsigned int nEqn,
	  							   const unsigned int nBC,
	  							   std::string folder_name);

	  		// in the following function we print the location from which we will be 
	  		// reading the system matrices.
	  		void fancy_print_filename(std::vector<std::string> &filenames);

			// we have a pointer to the kind of force we want in the system
			ForceType::Base_ForceType<dim> *force;

			// we will use them to initialize force 
			// the for
			//ce which corresponds to the exact solution of systemA
			ForceType::ForceType1<dim> force1;

			// the force which corresponds to the exact solution of systemB
			ForceType::ForceType2<dim> force2;

			// the force which corresponds to the poisson heat conduction problem
			ForceType::ForceType3<dim> force3;

			void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value);

			
			// matrices for the boundary conditions
			Full_matrix B_tilde_inv;
			Full_matrix B_hat;
			Full_matrix X_minus;
			Sparse_matrix BC;

			// first we create the class which handles matrix orperations
			MatrixOpt::Base_MatrixOpt matrixopt;

			// routines for tensor development
			TensorInfo::Base_TensorInfo<dim> base_tensorinfo;

			// development of projectors
			// build the projector corresponding to a particular normal vector
			Sparse_matrix build_Projector(const Tensor<1,dim> &normal_vector);
			Sparse_matrix build_InvProjector(const Tensor<1,dim> &normal_vector);

			// base class for the boundary routines. The decision regarding which boundary system to load
			// has to be done inside specific systems
			BCrhs::Base_BCrhs<dim> *base_bcrhs;


			// The following function already knows the Aminus_1D game
			Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector);

			// Matrices for the flux
			Full_matrix Aminus_1D;

			void build_BCrhs(const Tensor<1,dim,double> p,
									const Tensor<1,dim,double> normal_vector,
									Vector<double> &bc_rhs,
									const unsigned int b_id);


			//convert a matrix to a symmetric matrix
			void create_symmetric_matrix(Sparse_matrix &matrix,Sparse_matrix &S_half,
										 Sparse_matrix &S_half_inv);

			void reinit_system(const std::string &folder_name);

			// reinit the Aminus_1D matrix
			void reinit_Aminus1D();

			void reinit_Xminus();

			void symmetrize_system();

			double forcing_factor();

			void reinit_force();

			void reinit_BoundaryMatrices();

			virtual void reinit_BCrhs() = 0;

			double force_factor;
			void reinit_P(const std::string &folder_name);

	};

	template<int dim>
	Base_EquationGenerator<dim>
	::Base_EquationGenerator(const constant_data &constants)
	:
	constants(constants),
	force1(constants),
	force2(constants),
	force3(constants)
	{
		basefile.resize(dim + 2);

		unsigned int entry;

		for (unsigned int i = 0 ; i < dim ; i ++)
			basefile[i] = "A" + std::to_string(i + 1) + "_";

		entry = dim;
		Assert(entry < dim + 3,ExcNotInitialized());
		basefile[entry] = "B_";

		entry = dim + 1;
		Assert(entry < dim + 3,ExcNotInitialized());
		basefile[entry] = "odd_ID_";
	}

	template<int dim>
	void Base_EquationGenerator<dim>
	::build_triplet(std::vector<triplet> &Row_Col_Value,const std::string filename)
	{

		// as a precautionary measure, remove all the old value
		Row_Col_Value.clear();		

		// number of nonzeros in the system
		unsigned int nz;			

		// create an fstrem to read the file
		std::fstream in(filename.c_str());

		// string which will read the file line by line
		std::string line;

		// check whether we can open the file or not
		AssertThrow(in.is_open(),ExcMessage("Triplet cant be read, file could not be opened"));

		// get the first line in "in"
		std::getline(in,line);
		std::stringstream ss(line);
		ss >> nz;

		// shout out if the matrix is empty
		Assert(nz !=0 ,ExcMessage("Total number of non zeros should be > 0"));
		Row_Col_Value.reserve(nz);

		unsigned int counter = 0;
		while (getline(in, line))
		{
			std::stringstream ss2(line);

			unsigned int row;
			unsigned int col;
			double value;

			ss2 >> row;
			ss2 >> col;
			ss2 >> value;
			Row_Col_Value.push_back(triplet(row,col,value));

			counter ++;
		}

		// we check whether the correct number of values have been read or not
		AssertDimension(counter,nz);

	}

	template<int dim>
	void Base_EquationGenerator<dim>
	::fancy_print_filename(std::vector<std::string> &filenames)
	{
		const unsigned int num_matrices = filenames.size();
		const unsigned int Width = 15;

		std::cout << std::left << std::setw(Width) << "Matrix Name" << std::setw(Width) << "Read From" << std::endl;

		for(unsigned int i = 0 ; i < num_matrices ; i++)
			std::cout << std::left  << std::setw(Width)	<< basefile[i]
									<< std::setw(Width) << "From :"
								   << std::setw(Width) << filenames[i]
								   << std::endl;
	}

	// the following routines reads a vector from a file
	template<int dim>
	void
	Base_EquationGenerator<dim>
	::build_Vector(MatrixUI &vector,const std::string filename)
	{

		// number of entries in the vector
		unsigned int entries;			

		// create an fstrem to read the file
		std::fstream in(filename.c_str());

		// string which will read the file line by line
		std::string line;

		// check whether we can open the file or not
		AssertThrow(in.is_open(),ExcMessage("Vector cant be read, file could not be opened"));

		// get the first line in "in"
		std::getline(in,line);
		std::stringstream ss(line);
		ss >> entries;

		// shout out if the matrix is empty
		Assert(entries !=0 ,ExcMessage("Total number of entries should be > 0"));
		vector.resize(entries,1);

		unsigned int counter = 0;
		while (getline(in, line))
		{
			std::stringstream ss2(line);

			int value;
			ss2 >> value;

			Assert(value > 0,ExcMessage("Odd variable ID should not be less than zero"));
			vector(counter,0) = value;
			counter ++;
		}

		// we check whether the correct number of values have been read or not
		AssertDimension(counter,entries);

	}

	template<int dim>
	void Base_EquationGenerator<dim>
	::build_matrix_from_triplet(Sparse_matrix &matrix,std::vector<triplet> &Row_Col_Value)
	{
		// first we check whether Row_Col_Value has some size or not
		Assert(Row_Col_Value.size() != 0,ExcMessage("Triplet not constructed"));

		// Now we check whether the sparse matrix has been resized or not
		Assert(matrix.rows() != 0 && matrix.cols() !=0,ExcNotInitialized());

		matrix.setFromTriplets(Row_Col_Value.begin(),Row_Col_Value.end());

		// compress the system matrix
		matrix.makeCompressed();

		// clear the memory consumed by the triplet
		Row_Col_Value.clear();
	}

	template<int dim>
	void Base_EquationGenerator<dim>
	::print_matrix(Sparse_matrix &matrix,const std::string filename)
	{
		FILE *fp;
    	fp = fopen(filename.c_str(),"w+");
    	AssertThrow(fp != NULL,ExcMessage("Could not open file for writting. "));


		for (unsigned int i = 0 ; i < matrix.rows() ; i++)
		{
			for (unsigned int j = 0 ; j < matrix.cols() ; j++)
				 fprintf(fp, " %f ",matrix.coeffRef(i,j));

			fprintf(fp, "\n");
			
		}

		// coeffRef distorts the compressed form of matrices
		matrix.makeCompressed();
		fclose(fp);
	}

	// In the following function, we devlop the Jacobians of the flux, the square root of the 
	// Symmetrizer, the inverse of the square root of the symmetrizer and the boundary matrix.
	template<int dim>
	void Base_EquationGenerator<dim>
	::generate_matrices(equation_data &system_data,
	  					const unsigned int nEqn,
	  					const unsigned int nBC,
	  					const std::string folder_name)
	{
		// we need to now define the correct filenames corresponding to the 
		// system we are dealing with
		std::vector<std::string> basefile_system;

		// Shout out no basefile names have been allocated
		Assert(basefile.size() !=0,ExcMessage("basefile names not allocated"));

		// allocate basefile_system in case 
		basefile_system.resize(this->basefile.size());


		// allocate memory to all the matrices
		for (unsigned int i = 0 ; i < dim ; i ++)
			system_data.A[i].matrix.resize(nEqn,nEqn);

		system_data.Ax.matrix.resize(nEqn,nEqn);
		system_data.P.matrix.resize(nEqn,nEqn);
		system_data.B.matrix.resize(nBC,nEqn);
		system_data.Ax.matrix.resize(nEqn,nEqn);

		// now we loop over all the names in the basefile vector
		// by default at the system information is expected to be stored in a folder
		// labeled system_matrices
		for (unsigned int i = 0 ; i < this->basefile.size(); i++)
			basefile_system[i] = folder_name+this->basefile[i] + std::to_string(nEqn) + ".txt";

		bool print_fancy = false;		

		if (print_fancy)
		{
			std::cout << "Printing folder information " << std::endl;
			fancy_print_filename(basefile_system);
		}

		// since the number of jacobians = number of space dimensions, so that's what we do here
		for (unsigned int i = 0 ; i < dim ; i ++)
		{
			// develop the triplet
			this->build_triplet(system_data.A[i].Row_Col_Value,basefile_system[i]);

			// develop the matrix 
			this->build_matrix_from_triplet(system_data.A[i].matrix,system_data.A[i].Row_Col_Value);
		}

			system_data.Ax = system_data.A[0];

			// develop the B matrix
			this->build_triplet(system_data.B.Row_Col_Value,basefile_system[dim]);
			this->build_matrix_from_triplet(system_data.B.matrix,system_data.B.Row_Col_Value);

			// develop the ID of the odd variables
			this->build_Vector(system_data.odd_ID,basefile_system[dim + 1]);


	}

	// builds the Projector matrix to be used during computation
	template<int dim>
	Sparse_matrix
	Base_EquationGenerator<dim>
	::build_Projector(const Tensor<1,dim> &normal_vector)
	{

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		return(base_tensorinfo.reinit_global(nx,ny));
	}


	template<int dim>
	Sparse_matrix
	Base_EquationGenerator<dim>
	::build_InvProjector(const Tensor<1,dim> &normal_vector)
	{

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		return(base_tensorinfo.reinit_Invglobal(nx,ny));
	}

	template<int dim>
	Full_matrix 
	Base_EquationGenerator<dim>::
	build_Aminus(const Tensor<1,dim,double> normal_vector)
	{
		Full_matrix Aminus;
		Aminus.resize(this->constants.nEqn,this->constants.nEqn);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		Assert(this->Aminus_1D.rows() != 0 || this->Aminus_1D.cols() != 0,
			  ExcMessage("Aminus_1D has not been built yet"));

		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		return(base_tensorinfo.reinit_Invglobal(nx,ny) 
			   * this->Aminus_1D 
			   * base_tensorinfo.reinit_global(nx,ny));
	}

	// convert a matrix to a symmetric matrix
	template<int dim>
	void 
	Base_EquationGenerator<dim>::
	create_symmetric_matrix(Sparse_matrix &matrix,
							Sparse_matrix &S_half,
							Sparse_matrix &S_half_inv)
	{
		Assert(matrix.rows() != 0|| matrix.cols() != 0,ExcNotInitialized());
		Assert(S_half.rows() != 0|| S_half.cols() != 0,ExcNotInitialized());
		Assert(S_half_inv.rows() != 0|| S_half_inv.cols() != 0,ExcNotInitialized());
		
		// every symmetrix matrix = S_half * matrix * S_half_inv
		matrix = S_half * matrix * S_half_inv;
	}

	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::reinit_Aminus1D()
	{
		Aminus_1D = matrixopt.compute_Aminus(system_data.Ax.matrix);
	}

	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::reinit_Xminus()
	{
		X_minus = matrixopt.compute_Xminus(system_data.Ax.matrix,constants.nBC);
	}

	template<int dim>
	void
	Base_EquationGenerator<dim>
	::symmetrize_system()
	{
		// now we symmetrize the system
		// First symmetrize the jacobians
		for (unsigned int i = 0 ; i < dim ; i ++)
			create_symmetric_matrix(system_data.A[i].matrix,base_tensorinfo.S_half,base_tensorinfo.S_half_inv);

		// now symmetrize the production term
		create_symmetric_matrix(system_data.P.matrix,base_tensorinfo.S_half,base_tensorinfo.S_half_inv);
	}

	// The following routine computes S_half.f where f is the external forcing into the system
	template<int dim>
	double
	Base_EquationGenerator<dim>
	::forcing_factor()
	{
		const unsigned int var_force = constants.variable_map.find(constants.force_variable)->second;

		Assert(base_tensorinfo.S_half.rows() != 0 || base_tensorinfo.S_half.cols() != 0,ExcNotInitialized());
		const double force_factor = base_tensorinfo.S_half.coeffRef(var_force,var_force);

		base_tensorinfo.S_half.makeCompressed();

		return(force_factor);
	}

	template<int dim>
	void
	Base_EquationGenerator<dim>
	::reinit_force()
	{
		switch(constants.force_type)
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
	}

	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::reinit_BoundaryMatrices()
	{
		
		// develop the matrices which are independent of the test case
		switch(constants.bc_type)
		{
			// characteristic boundary conditions
			case characteristic:
			{

				BoundaryHandler::Base_BoundaryHandler_Char<dim> boundary_handler_char(system_data.Ax.matrix,
																				system_data.B.matrix,
																				constants.nBC);



				B_tilde_inv.resize(constants.nBC,constants.nBC);
				B_hat.resize(constants.nEqn,this->constants.nEqn);

				B_tilde_inv = boundary_handler_char.build_B_tilde_inv();
				B_hat = boundary_handler_char.build_B_hat(this->B_tilde_inv);

				break;
			}

			// picking up the odd variables
			case odd:
			{
				BoundaryHandler::Base_BoundaryHandler_Odd<dim> boundary_handler_odd(system_data.B.matrix,
																				   system_data.odd_ID);

				 BC = boundary_handler_odd.develop_BC();
				break;
			}

			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here"));
				break;
			}
		}

	}



	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::source_term(const std::vector<Point<dim>> &p,
				 std::vector<Vector<double>> &value)
	{

		// now we simply pass on the value to the base class
		force->source_term(p,value,force_factor);
	}

	template<int dim>
	void
	Base_EquationGenerator<dim>
	::build_BCrhs(const Tensor<1,dim,double> p,
				const Tensor<1,dim,double> normal_vector,
				Vector<double> &bc_rhs,
				const unsigned int b_id)
	{

		base_bcrhs->BCrhs(p,normal_vector,bc_rhs,b_id);
	}

	template<>
	void
	Base_EquationGenerator<2>
	::reinit_P(const std::string &folder_name)
	{
		AssertDimension(system_data.P.matrix.rows(),constants.nEqn);
		AssertDimension(system_data.P.matrix.cols(),constants.nEqn);

		if (constants.nEqn == 6)
		{
			for (int i = 1 ; i < constants.nEqn ; i ++)
				system_data.P.matrix.coeffRef(i,i) = 1/constants.tau;

			
		}
		else
		{
			const unsigned int num_conserved = 4;

			switch(constants.coll_op)
			{
				case BGK:
				{
					for (int i =  num_conserved; i < constants.nEqn ; i++)
						system_data.P.matrix.coeffRef(i,i) = 1/constants.tau;					

					break;
				}
				case Boltzmann_MM:
				{
					std::string file_for_P = folder_name + "MM_P_"+std::to_string(constants.nEqn)+".txt";
					this->build_triplet(system_data.P.Row_Col_Value,file_for_P);

					// develop the matrix 
					this->build_matrix_from_triplet(system_data.P.matrix,system_data.P.Row_Col_Value);
					system_data.P.matrix /= constants.tau;
					break;
				}
				default:
				{
					Assert(1 == 0,ExcMessage("Should not have reached here"));
					break;
				}
			}
			
		}

		system_data.P.matrix.makeCompressed();
	}

	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::reinit_system(const std::string &folder_name)
	{
		this->generate_matrices(this->system_data,this->constants.nEqn,this->constants.nBC,folder_name);

		// we now develop the P matrix
		reinit_P(folder_name);	


		// if the number of equations are not equal to 6 then we fix the normal relaxation velocity
		if (constants.nEqn != 6)
		{
				const bool fix_B = true;
				const bool check_B = true;

				BoundaryHandler::Base_BoundaryHandler_Char<dim>::fix_B_vx(constants.epsilon,fix_B,system_data.B.matrix);

				// check whether B has been fixed or not
				if(check_B)
					for (unsigned int i = 0 ; i < system_data.B.matrix.cols() ; i++)
					 if(i != constants.variable_map.find("vx")->second)			
						Assert(fabs(system_data.B.matrix.coeffRef(0,i)) < 1e-3,ExcMessage("relaxational normal velocity not accommodated in B"));
		}

		// we need to compres B again
		system_data.B.matrix.makeCompressed();
	}
}
