// the following is the base of the equation generator, it contains all the routines which
// which are equation dependent. This includes the data for the equations and the Projector matrices.
// It is an abstract base class with a specifiied routine for the reading of the matrices. 
namespace EquationGenerator
{
	using namespace dealii;

	template<int dim>
	class Base_EquationGenerator
	{
		protected:
			Base_EquationGenerator();
			
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

			};

			// the base filenames for the system matrices
			std::vector<std::string> basefile;

			// Given the filename, read the tripet for the sparse matrix 
			void build_triplet(std::vector<triplet> &Row_Col_Value,
								const std::string filename);

			// Develop the matrix from the triplet
	  		void build_matrix_from_triplet(Sparse_matrix &matrix,
	  									   std::vector<triplet> &Row_Col_Value);

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


	};

	template<int dim>
	Base_EquationGenerator<dim>
	::Base_EquationGenerator()
	{
		basefile.resize(dim + 3);

		for (unsigned int i = 0 ; i < dim ; i ++)
			basefile[i] = "A" + std::to_string(i + 1) + "_";

		basefile[dim] = "B_";
		basefile[dim+1] = "S_half_";
		basefile[dim + 2] = "S_half_inv_";
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
		

		fancy_print_filename(basefile_system);

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
	}

}
