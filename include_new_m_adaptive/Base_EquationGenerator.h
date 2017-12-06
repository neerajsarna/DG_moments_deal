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
			Base_EquationGenerator(const constant_numerics &constants,
								   const int nEqn,
								   const int nBC,
								   const int Ntensors);
			
			DeclException1 (ExcFileNotOpen, std::string,
                        << "Could not open " << arg1 );

			const constant_numerics constants;
			const int nEqn;
			const int nBC;
			const int Ntensors;

			struct equation_data
			{

			// matrix to store the jacobian of the flux
				 system_matrix A[dim];

		  	// The Production term
				  system_matrix P;

			// The Jacobian in the x-direction
				system_matrix Ax;

			// The matrix for the boundary, full accomodation
				system_matrix B;


			// The matrix for the prescribed inflow
				system_matrix B_prescribedInflow;

			// The matrix for the boundary conditions, specular reflection
				system_matrix B_specular;

			// penalty matrix for the odd boundary conditions 
				system_matrix Sigma; 

			// ID of the odd variables 
				MatrixUI odd_ID;

			// matrix which models the inflow
				system_matrix Binflow;


			// matrix containing the rhoW and rhoInflow
				system_matrix rhoW;
				system_matrix rhoInflow;

			// the A mod coming from the kinetic scheme
				system_matrix Amod_kinetic;

			};

			bool initialized_system = false;

			// a data structure to store all the system data
			equation_data system_data;

			BCrhs::Base_BCrhs<dim> *base_bcrhs;


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
			
			// penalty matrices for the characteristic implementation of boundary conditions
			Full_matrix penalty_char_wall;
			Full_matrix penalty_char_inflow;
			Full_matrix X_minus;


			// first we create the class which handles matrix orperations
			MatrixOpt::Base_MatrixOpt matrixopt;

			// routines for tensor development
			TensorInfo::Base_TensorInfo<dim> base_tensorinfo;

			// development of projectors
			// build the projector corresponding to a particular normal vector
			Sparse_matrix build_Projector(const Tensor<1,dim> &normal_vector);
			Sparse_matrix build_InvProjector(const Tensor<1,dim> &normal_vector);


			// The following function already knows the Aminus_1D game
			Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector);

			// build the Aminus matrix using a kinetic flux
			Full_matrix build_Aminus_kinetic(const Tensor<1,dim,double> normal_vector);

			Full_matrix build_An(const Tensor<1,dim,double> normal_vector);

			// Matrices for the flux
			Full_matrix Aminus_1D;

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

			void reinit_penalty_matrices(const Sparse_matrix &Ax, Sparse_matrix &B,Full_matrix &penalty_matrix);
			void reinit_BoundaryMatrices();

			void reinit_Bspecular();

			double force_factor;
			void reinit_P(const std::string &folder_name);
			void reinit_rhoInflow();

	};

	template<int dim>
	Base_EquationGenerator<dim>
	::Base_EquationGenerator(const constant_numerics &constants,
				 			 const int nEqn,
							 const int nBC,
							 const int Ntensors)
	:
	constants(constants),
	nEqn(nEqn),
	nBC(nBC),
	Ntensors(Ntensors),
	force1(constants,nEqn),
	force2(constants,nEqn),
	force3(constants,nEqn),
	base_tensorinfo(nEqn,Ntensors)
	{
		// the maximum number of matrices in the system
		// 3 flux matrices (Ax,Ay and Az)
		// 2 boundary matrix (B and B_inflow)
		// 1 sigma (Penalty matrix Sigma)
		// 1 odd variables (Vector containing all the odd variables)
		const int max_matrices = dim + 7;

		// ending for the name of different matrices
		std::string name_ending = + "_" + std::to_string(dim)+"D_";

		// the base file names, to be updated with the equation number later on
		basefile.resize(max_matrices);

		Assert(basefile.size() == max_matrices,ExcMessage("incorrect vector size"));

		unsigned int entry;

		for (unsigned int i = 0 ; i < dim ; i ++)
			basefile[i] = "A" + std::to_string(i + 1) + name_ending;

		entry = dim;
		Assert(entry < max_matrices,ExcNotInitialized());
		basefile[entry] = "B"+name_ending;

		entry = dim + 1;
		Assert(entry < max_matrices,ExcNotInitialized());
		basefile[entry] = "odd_ID"+name_ending;

		entry = dim + 2;
		Assert(entry < max_matrices,ExcNotInitialized());
		basefile[entry] = "Sigma" + name_ending;

		entry = dim + 3;
		Assert(entry < max_matrices,ExcNotInitialized());
		basefile[entry] = "Binflow" + name_ending;

		entry = dim + 4;
		basefile[entry] = "rhoW" + name_ending;

		entry = dim + 5;
		basefile[entry] = "rhoInflow" + name_ending;

		entry = dim + 6;
		basefile[entry] = "kineticAmod" + name_ending;


		// we can check the number of equations and the number of boundary conditions for the 1D case
		if (dim == 1)
		{
					// now we check the number of boundary conditions 
		// if we have even number of variables in the system
		if (nEqn %2 == 0)
			Assert(nBC == (nEqn / 2),ExcMessage("Incorrect number of boundary conditions"));

		// if we have odd number of variables in the system
		if (nEqn %2 != 0)
			Assert(nBC == (nEqn - 1)  / 2,ExcMessage("Incorrect number of boundary conditions"));	
		}
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
		AssertThrow(in.is_open(),ExcFileNotOpen(filename));

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

	// In the following function we build a sparse matrix from a given triplet which has been read from a file.
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

	// In the following function, we develop the Jacobians of the flux, the square root of the 
	// Symmetrizer, the inverse of the square root of the symmetrizer and the boundary matrix.
	// No specialization needed for different dimensions.
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

		// Shout out if no basefile names have been allocated
		Assert(basefile.size() !=0,ExcMessage("basefile names not allocated"));

		// allocate basefile_system in case 
		basefile_system.resize(this->basefile.size());


		// allocate memory to all the matrices
		for (unsigned int i = 0 ; i < dim ; i ++)
			system_data.A[i].matrix.resize(nEqn,nEqn);

		// allocate the memory for all the other matrices
		system_data.Ax.matrix.resize(nEqn,nEqn);
		system_data.P.matrix.resize(nEqn,nEqn);
		system_data.B.matrix.resize(nBC,nEqn);
		system_data.B_prescribedInflow.matrix.resize(nBC,nEqn);
		system_data.Sigma.matrix.resize(nEqn,nBC);
		system_data.Binflow.matrix.resize(nBC,nEqn);

		system_data.rhoW.matrix.resize(nEqn,nEqn);
		system_data.rhoInflow.matrix.resize(nEqn,nEqn);
		system_data.Amod_kinetic.matrix.resize(nEqn,nEqn);

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

			// we store the unsymmetrized flux matrix for numerical flux generation
			system_data.Ax = system_data.A[0];

			// develop the B matrix
			this->build_triplet(system_data.B.Row_Col_Value,basefile_system[dim]);
			this->build_matrix_from_triplet(system_data.B.matrix,system_data.B.Row_Col_Value);
			this->build_triplet(system_data.B.Row_Col_Value,basefile_system[dim]);
			this->build_matrix_from_triplet(system_data.B_prescribedInflow.matrix,system_data.B.Row_Col_Value);


			// develop the ID of the odd variables
			this->build_Vector(system_data.odd_ID,basefile_system[dim + 1]);


			// develop the Sigma matrix
			this->build_triplet(system_data.Sigma.Row_Col_Value,basefile_system[dim+2]);
			this->build_matrix_from_triplet(system_data.Sigma.matrix,system_data.Sigma.Row_Col_Value);

			this->build_triplet(system_data.Binflow.Row_Col_Value,basefile_system[dim+3]);
			this->build_matrix_from_triplet(system_data.Binflow.matrix,system_data.Binflow.Row_Col_Value);

			// this->build_triplet(system_data.rhoW.Row_Col_Value,basefile_system[dim+4]);
			// this->build_matrix_from_triplet(system_data.rhoW.matrix,system_data.rhoW.Row_Col_Value);
		
		// we only need this if the dimension is not equal to 1. When dim ==1, the upwind flux itself corresponds to
	    // a kinetic flux.
		if(dim != 1)
		{
			this->build_triplet(system_data.Amod_kinetic.Row_Col_Value,basefile_system[dim+6]);
			this->build_matrix_from_triplet(system_data.Amod_kinetic.matrix,system_data.Amod_kinetic.Row_Col_Value);			
		}


		system_data.B_prescribedInflow.matrix = system_data.B.matrix;
			
	}

	// builds the Projector matrix to be used during computation. Specialization for the 2D case
	template<int dim>
	Sparse_matrix
	Base_EquationGenerator<dim>
	::build_Projector(const Tensor<1,dim> &normal_vector)
	{


		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		Assert(initialized_system,ExcMessage("You are trying to access Projector development without initializing the system. Please initialize the system first."));

		return(base_tensorinfo.reinit_global(normal_vector));
	}


	// specialization for the 2D case
	template<int dim>
	Sparse_matrix
	Base_EquationGenerator<dim>
	::build_InvProjector(const Tensor<1,dim> &normal_vector)
	{


		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		return(base_tensorinfo.reinit_Invglobal(normal_vector));
	}

	
	template<int dim>
	Full_matrix 
	Base_EquationGenerator<dim>::
	build_Aminus(const Tensor<1,dim,double> normal_vector)
	{

		Assert(this->Aminus_1D.rows() != 0 || this->Aminus_1D.cols() != 0,
			  ExcMessage("Aminus_1D has not been built yet"));

		Assert(base_tensorinfo.varIdx.rows() != 0 || base_tensorinfo.varIdx.cols() !=0 ,
				ExcMessage("Base tensor info not initialized"));

		return(base_tensorinfo.reinit_Invglobal(normal_vector) 
			   	* this->Aminus_1D 
			   	* base_tensorinfo.reinit_global(normal_vector));

	}

	// Aminus using kinetic flux, for the 1D case it is the same as the normal Aminus
	template<>
	Full_matrix
	Base_EquationGenerator<1>::
	build_Aminus_kinetic(const Tensor<1,1,double> normal_vector)
	{
		// the normal upwind flux is the kinetic flux for this case due to the underlying discrete veloctiy 
		// structure.
		return(build_Aminus(normal_vector));
	}


	template<>
	Full_matrix
	Base_EquationGenerator<2>::
	build_Aminus_kinetic(const Tensor<1,2,double> normal_vector)
	{
		// the factor of half is taken care of separately. The symmetrizer is taken care of by the projector itself
		return(base_tensorinfo.reinit_Invglobal(normal_vector) 
			   	* (system_data.Amod_kinetic.matrix - system_data.Ax.matrix)
			   	* base_tensorinfo.reinit_global(normal_vector));
	}

	// the following matrix returns the flux in a particular direction
	template<int dim>
	Full_matrix
	Base_EquationGenerator<dim>::
	build_An(const Tensor<1,dim,double> normal_vector)
	{
		return(base_tensorinfo.reinit_Invglobal(normal_vector) 
			   	* system_data.Ax.matrix
			   	* base_tensorinfo.reinit_global(normal_vector));
	}

	// Convert a matrix to a symmetric matrix. The following function returns 
	// S_half * matrix * S_half_inv. So it symmetrizes the system under a similarity transform
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
		X_minus = matrixopt.compute_Xminus(system_data.Ax.matrix,nBC);
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

	
	template<int dim>
	double
	Base_EquationGenerator<dim>
	::forcing_factor()
	{
		const unsigned int var_force = constants.variable_map[dim-1].find(constants.force_variable)->second;

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
		bool initialized_force = false;
		switch(constants.force_type)
		{
			case type1:
			{
				initialized_force = true;
				force = &force1;
				break;
			}

			case type2:
			{
				initialized_force = true;
				force = &force2;
				break;
			}

			case type3:
			{
				initialized_force = true;
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
	::
	reinit_penalty_matrices(const Sparse_matrix &Ax,Sparse_matrix &B,Full_matrix &penalty_matrix)
	{
		BoundaryHandler::Base_BoundaryHandler_Char boundary_handler_char(Ax,
																		 B,
																		 nBC);



		penalty_matrix = boundary_handler_char.build_penalty_char();

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

				//boundary matrices for the wall boundary
				reinit_penalty_matrices(system_data.A[0].matrix,system_data.B.matrix,penalty_char_wall);

				// boundary matrices for the inflow boundary
				reinit_penalty_matrices(system_data.A[0].matrix,system_data.Binflow.matrix,penalty_char_inflow);

				break;
			}

			// the penalty matrix for the odd implementation is read from a file
			case odd:
			{
				break;
			}

			// we don't need to do anything for the kinetic boundary implementation
			case kinetic:
			{
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


	// specialization for the 2D case
	template<int dim>
	void
	Base_EquationGenerator<dim>
	::reinit_P(const std::string &folder_name)
	{
		AssertDimension(system_data.P.matrix.rows(),nEqn);
		AssertDimension(system_data.P.matrix.cols(),nEqn);


		// the total number of conserved variables
		const unsigned int num_conserved = dim + 2;

		// system A
		if (nEqn == 6 && dim == 2)
		{
			for (int i = 1 ; i < nEqn ; i ++)
				system_data.P.matrix.coeffRef(i,i) = 1/constants.tau;


		}
		else
		{
			switch(constants.coll_op)
			{
				case BGK:
				{
					for (int i =  num_conserved; i < nEqn ; i++)
						system_data.P.matrix.coeffRef(i,i) = 1/constants.tau;					

					break;
				}
				case Boltzmann_MM:
				{
					Assert(dim != 1,ExcMessage("Not implemented for 1D systems"));
					std::string file_for_P = folder_name + "MM_P"+"_"+std::to_string(dim) + "D_"
											+std::to_string(nEqn)+".txt";
					this->build_triplet(system_data.P.Row_Col_Value,file_for_P);

					// develop the matrix 
					this->build_matrix_from_triplet(system_data.P.matrix,system_data.P.Row_Col_Value);

					// need to scale the production terms as per the Knudsen number
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

	// we do not precsribe any particular moment at the inflow or the outflow therefore rho inflow is simply zero
	template<int dim>
	void
	Base_EquationGenerator<dim>
	::reinit_rhoInflow()
	{

		for (int i = 0 ; i < nEqn ; i++)
			system_data.rhoInflow.matrix.coeffRef(i,i) = 0;

		system_data.rhoInflow.matrix.makeCompressed();
	}


	// initialize the system_data corresponding to a particular moment system
	// Specialization for the 2D case
	template<int dim>
	void 
	Base_EquationGenerator<dim>
	::reinit_system(const std::string &folder_name)
	{
		this->generate_matrices(this->system_data,this->nEqn,this->nBC,folder_name);

		// we now develop the P matrix
		reinit_P(folder_name);	

		// we now develop the rhoInflow matrix, which has nothin but zeros on the diagonal
		reinit_rhoInflow();

		//don't fix the boundary conditions for the A system
		if (nEqn == 6 && dim ==2)
		{
				// do nothing
		}
		else
		{

			const bool fix_B = true;
			const bool check_B = true;

			std::cout << "value of epsilon " << constants.epsilon << std::endl;

			// we fix the wall boundary and the prescribed inflow boundary with different epsilons
			BoundaryHandler::Base_BoundaryHandler_Char::fix_B_vx(constants.epsilon,
																fix_B,system_data.B.matrix);

			// the epsilon for the inflow
			BoundaryHandler::Base_BoundaryHandler_Char::fix_B_vx(constants.epsilon_inflow,
																fix_B,system_data.B_prescribedInflow.matrix);
				// check whether B has been fixed or not
			if(check_B)
				for (unsigned int i = 0 ; i < system_data.B.matrix.cols() ; i++)
					// in both the cases, 1D, 2D or 3D the index for vx does not changed and remains at one
					if(i != constants.variable_map[dim-1].find("vx")->second)	
						Assert(fabs(system_data.B.matrix.coeffRef(0,i)) < 1e-3,
								ExcMessage("relaxational normal velocity not accommodated in B"));						
				
		}


		// 

		// now we can initialize the boundary matrix for specular reflection
		reinit_Bspecular();



		// we need to compress B again so that it consumes less memory
		system_data.B.matrix.makeCompressed();

	}



	template<int dim>
	void
	Base_EquationGenerator<dim>
	::reinit_Bspecular()
	{
		system_data.B_specular.matrix.resize(nBC,nEqn);

		// we initialize the specular matrix with the same matrix as full accommodation
		system_data.B_specular.matrix = system_data.B.matrix;

		bool odd_variable = false;

		// we do not change the coefficients of the first row since we need them for stability reasons
		for (int i = 1 ; i < nBC ; i++)
		{

			for (unsigned int j = 0 ; j < (unsigned int)nEqn ; j++)
			{
				odd_variable = false;

				// if the coefficient we have captured corresponds to an odd variable then do not do anything, else we 
				// put the coefficient in the boundary matrix to zero.
				for (int k = 0 ; k < system_data.odd_ID.rows() ; k++)
					if (j == system_data.odd_ID(k,0))
					{
						odd_variable = true;
						break;
					}

				// put all the even variables to zero
				if (!odd_variable)
					system_data.B_specular.matrix.coeffRef(i,j) = 0;
						
			}

		}

		
		system_data.B_specular.matrix.makeCompressed();

	}
	

}
