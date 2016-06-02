namespace EquationGenerator
{
	using namespace dealii;
	using namespace std;
	using namespace Basics;

	

	template<int num_flux,int dim>
	class Base_EquationGenerator:public Base_Basics
	{
	   public:
		const nEqn_data num_equations;
		const System_Type system_type;
		const Force_Type force_type;

	  struct equation_data
	  {
	  		bool is_symmetric;
	  		unsigned int nEqn;
			system_matrix A[dim];
			system_matrix P;
			system_matrix BC;
			system_matrix S_half;			// sqrt(S)
			system_matrix S_half_inv;
			system_matrix Ax;

	  		Full_matrix Aminus_1D_Int;
	  		Full_matrix Aminus_1D_Bound;
	  };

		
	  
	  		Base_EquationGenerator(nEqn_data const&num_equations,
	  							   physical_data &physical_constants,
	  							   string &output_dir);

	  		vector<equation_data> system_data;

	  		void build_BCrhs(const Tensor<1,dim,double> p,
							const Tensor<1,dim,double> normal_vector,Vector<double> &bc_rhs,
							const unsigned int system_id) const; 

			Sparse_matrix build_Projector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id) const;

			void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value,const unsigned int system_id);
			Sparse_matrix build_InvProjector(const Tensor<1,dim,double> normal_vector,
											const unsigned int system_id) const;		
			Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector,
										const unsigned int system_id) ;		
	
	  protected:
	  	void build_triplet(system_matrix &matrix_info,const string filename);
	  	void build_matrix_from_triplet(system_matrix &matrix_info);
	  	void print_matrix(const system_matrix matrix_info,const string filename);

		void generate_matrices(equation_data &system_data,const unsigned int system_id);
		Tensor<1,dim,double> mirror(const Tensor<1,dim,double> normal_vector) const;			
		void build_Aminus1D(Full_matrix &Aminus_1D_Int,
							Full_matrix &Aminus_1D_Bound,
							const unsigned int system_id);										
		void build_P(system_matrix &P,const unsigned int system_id);												
		void build_BC(system_matrix &BC,const unsigned int system_id);											
	};

	template<int num_flux,int dim> 
	void Base_EquationGenerator<num_flux,dim>::
	build_triplet(system_matrix &matrix_info,const string filename)
	{

		cout << "......reading triplet from: " << filename << endl;
		matrix_info.Row_Col_Value.clear();						// remove all the old values from the input vector
		unsigned int nz;							// number of non-zeros in the system

		fstream in(filename.c_str());
		string line;
		assert(in.is_open());

		getline(in,line);
		stringstream ss(line);
		ss >> nz;

		matrix_info.Row_Col_Value.reserve(nz);

		unsigned int counter = 0;
		while (getline(in, line))
		{
			stringstream ss2(line);

			unsigned int row;
			unsigned int col;
			double value;

			ss2 >> row;
			ss2 >> col;
			ss2 >> value;
			matrix_info.Row_Col_Value.push_back(triplet(row,col,value));

			counter ++;
		}

		assert(counter == nz);

	}

	template<int num_flux,int dim>
	void Base_EquationGenerator< num_flux, dim>::
	build_matrix_from_triplet(system_matrix &matrix_info)
	{
		cout << "developing matrix from triplet......" << endl;
		assert(matrix_info.Row_Col_Value.size() != 0);
		assert(matrix_info.matrix.cols() != 0 || matrix_info.matrix.rows() != 0);
		matrix_info.matrix.setFromTriplets(matrix_info.Row_Col_Value.begin(), matrix_info.Row_Col_Value.end());
		
	}

	template< int num_flux,int dim> 
	void Base_EquationGenerator<num_flux,dim>::
	print_matrix(system_matrix matrix_info,const string filename_to_write)
	{
    	FILE *fp;
    	fp = fopen(filename_to_write.c_str(),"w+");
    	assert(fp != NULL);


		for (unsigned int i = 0 ; i < matrix_info.matrix.rows() ; i++)
		{
			for (unsigned int j = 0 ; j < matrix_info.matrix.cols() ; j++)
				 fprintf(fp, " %f ",matrix_info.matrix.coeffRef(i,j));

			fprintf(fp, "\n");
			
		}

		fclose(fp);
		cout << "writting the read matrices to: " << filename_to_write<< "\n" << endl;
	}


	/*The constructor of the class*/
	template<int num_flux,int dim>
	 Base_EquationGenerator<num_flux,dim>
	 ::Base_EquationGenerator(nEqn_data const&num_equations,
	 					      physical_data &physical_constants,
	 					      string &output_dir)
	 :
	 Base_Basics(physical_constants,output_dir),
	 num_equations(num_equations),
	 system_type(num_equations.system_type),
	 force_type(num_equations.force_type)
	{
		cout << "loading " << num_equations.no_of_total_systems << " systems" << endl;
		system_data.resize(num_equations.no_of_total_systems);

		for (unsigned int i = 0 ; i < num_equations.no_of_total_systems ; i++)
			generate_matrices(system_data[i],i);
					
	}

	#include "develop_systems.h"
//	#include "generate_systemB.h"
	
}
