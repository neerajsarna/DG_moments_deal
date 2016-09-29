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

			// Square root of the symmetrizer
				system_matrix S_half;			// sqrt(S)

			// Inverse of the sqaure root of the symmetrizer
				system_matrix S_half_inv;

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

			// The matrix from the production term
			virtual void build_P(Sparse_matrix &P) = 0;

			MatrixUI varIdx;

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

			virtual void source_term(const std::vector<Point<dim>> &p,
									 std::vector<Vector<double>> &value) = 0;

			
			// matrices for the boundary conditions
			Full_matrix B_tilde_inv;
			Full_matrix B_hat;
			Sparse_matrix BC;

			// development of projectors
			// build the projector corresponding to a particular normal vector
			virtual Sparse_matrix build_Projector(const Tensor<1,dim> &normal_vector) = 0;
			virtual Sparse_matrix build_InvProjector(const Tensor<1,dim> &normal_vector) = 0;

			// this system is stupid and special. So special treatment of the boundary
			BCrhs_systemA::Base_BCrhs_systemA<dim> *base_bcrhs;
			BCrhs_systemA::BCrhs_ring_char_systemA<dim> bcrhs_ring_char;
			BCrhs_systemA::BCrhs_ring_odd_systemA<dim> bcrhs_ring_odd;

			virtual Sparse_matrix build_Aminus(const Tensor<1,dim,double> normal_vector) = 0;

			// Matrices for the flux
			Sparse_matrix Aminus_1D;

			virtual void build_BCrhs(const Tensor<1,dim,double> p,
									const Tensor<1,dim,double> normal_vector,
									Vector<double> &bc_rhs) = 0;
	};

	template<int dim>
	Base_EquationGenerator<dim>
	::Base_EquationGenerator(const constant_data &constants)
	:
	constants(constants),
	force1(constants),
	force2(constants),
	force3(constants),
	bcrhs_ring_char(constants),
	bcrhs_ring_odd(constants)
	{
		basefile.resize(dim + 4);
		unsigned int entry;

		for (unsigned int i = 0 ; i < dim ; i ++)
			basefile[i] = "A" + std::to_string(i + 1) + "_";

		entry = dim;
		Assert(entry < dim + 3,ExcNotInitialized());
		basefile[entry] = "B_";

		entry = dim + 1;
		Assert(entry < dim + 3,ExcNotInitialized());
		basefile[entry] = "S_half_";

		entry = dim + 2;
		Assert(entry < dim + 3,ExcNotInitialized());
		basefile[entry] = "S_half_inv_";

		entry = dim + 3;
		Assert(entry <= dim + 3,ExcNotInitialized());
		basefile[entry] = "odd_ID_";
	};

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
		system_data.S_half.matrix.resize(nEqn,nEqn);
		system_data.S_half_inv.matrix.resize(nEqn,nEqn);
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

			// develop S_half
			this->build_triplet(system_data.S_half.Row_Col_Value,basefile_system[dim+1]);
			this->build_matrix_from_triplet(system_data.S_half.matrix,system_data.S_half.Row_Col_Value);

			// develop S_half_inv
			this->build_triplet(system_data.S_half_inv.Row_Col_Value,basefile_system[dim+2]);
			this->build_matrix_from_triplet(system_data.S_half_inv.matrix,system_data.S_half_inv.Row_Col_Value);

			// develop the ID of the odd variables
			this->build_Vector(system_data.odd_ID,basefile_system[dim + 3]);
	}

}
