namespace EquationGenerator
{
	using namespace dealii;
	using namespace std;
	using namespace Basics;

	

	template<int dim> class Base_EquationGenerator
	{
		// 
	  public:	
		virtual Sparse_matrix build_Projector(const Tensor<1,dim,double> normal_vector) = 0;
		virtual Sparse_matrix build_InvProjector(const Tensor<1,dim,double> normal_vector) = 0;
		virtual MatrixXd build_Aminus(const Tensor<1,dim,double> normal_vector)  =0;
		virtual void build_BCrhs(const Tensor<1,dim,double> p,const Tensor<1,dim,double> normal_vector,Vector<double> &bc_rhs) = 0;
		virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) = 0;
		bool exists(vector<triplet> Row_Col_Value,const int row_index,const int col_index) ;

	  protected:
	  	virtual void build_Aminus1D() = 0;
	  	void build_triplet(system_matrix &matrix_info,const string filename);
	  	void build_matrix_from_triplet(system_matrix &matrix_info);
	  	void print_matrix(const system_matrix matrix_info,const string filename);
		virtual void build_P(system_matrix &P) = 0;
		virtual void build_BC(system_matrix &BC) = 0;
		virtual Tensor<1,dim,double> mirror(const Tensor<1,dim,double> normal_vector)  = 0;
		/*virtual void build_Aminus1D(const Sparse_matrix Ax) const = 0; */
	};

	template<int dim> void Base_EquationGenerator<dim>::build_triplet(system_matrix &matrix_info,const string filename)
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
		cout << "done reading triplet..." << endl;

	}

	template<int dim> void Base_EquationGenerator<dim>::build_matrix_from_triplet(system_matrix &matrix_info)
	{
		cout << "developing matrix from triplet......" << endl;
		assert(matrix_info.Row_Col_Value.size() != 0);
		assert(matrix_info.matrix.cols() != 0 || matrix_info.matrix.rows() != 0);

		matrix_info.matrix.setFromTriplets(matrix_info.Row_Col_Value.begin(), matrix_info.Row_Col_Value.end());
		cout << "done developing matrix......" << endl;
		
	}

	template<int dim> bool Base_EquationGenerator<dim>::exists(vector<triplet> Row_Col_Value,const int row_index,const int col_index) 
	{
		assert(Row_Col_Value.size() != 0);
		vector<triplet>::iterator it = Row_Col_Value.begin(),it_end = Row_Col_Value.end();

		for (; it != it_end ;it++)
			if (row_index == it->row() && col_index == it->col())
				return true;
			
		return false;

	}


	template<int dim> void Base_EquationGenerator<dim>::print_matrix(system_matrix matrix_info,const string filename_to_write)
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


	#include "generate_systemA.h"
	#include "generate_systemB.h"
	
}
