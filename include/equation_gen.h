namespace EquationGenerator
{
	using namespace dealii;
	using namespace std;
	using namespace Basics;

	template<int dim> class Base_EquationGenerator
	{
		// 
	  public:	
		virtual Sparse_matrix build_Projector(const Point<dim> normal_vector) const = 0;
		virtual Sparse_matrix build_InvProjector(const Point<dim> normal_vector) const = 0;
		virtual MatrixXd build_Aminus(const Point<dim> normal_vector) const =0;
		virtual void build_BCrhs(const Point<dim> p,const Point<dim> normal_vector,Vector<double> &bc_rhs)const = 0;
		virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) const = 0;
		bool exists(vector<triplet> Row_Col_Value,const int row_index,const int col_index) const;

	  protected:
	  	void build_triplet(system_matrix &matrix_info,const string filename);
	  	void build_matrix_from_triplet(system_matrix &matrix_info);
	  	void print_matrix(const system_matrix matrix_info,const string filename);
		virtual void build_P(system_matrix &P) = 0;
		virtual void build_BC(system_matrix &BC) = 0;
		virtual Point<dim> mirror(const Point<dim> normal_vector) const = 0;
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

	template<int dim> bool Base_EquationGenerator<dim>::exists(vector<triplet> Row_Col_Value,const int row_index,const int col_index) const
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

	// need to now read the triplets for the sparse matrices from a file
	template<int dim> class generate_systemA:public Base_EquationGenerator<dim>,protected constants_SystemA
	{
		public:
			generate_systemA();
			system_matrix A[dim];
			system_matrix P;
			system_matrix BC;
			Full_matrix Aminus1D;

			virtual void build_BCrhs(const Point<dim> p,const Point<dim> normal_vector,Vector<double> &bc_rhs) const;
			virtual Sparse_matrix build_Projector(const Point<dim> normal_vector) const;
			virtual Sparse_matrix build_InvProjector(const Point<dim> normal_vector) const;
			virtual Full_matrix build_Aminus(const Point<dim> normal_vector) const;
			virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) const;		

	  protected:
	  		string generate_filename_to_write(const string folder_name,const string filename);
			virtual void build_P(system_matrix &P);
			virtual void build_BC(system_matrix &BC);
			virtual Point<dim> mirror(const Point<dim> normal_vector) const ;
			/*virtual void build_Aminus1D(const Sparse_matrix Ax) const; */
	};

	template<int dim> void generate_systemA<dim>::source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) const
	{
		assert(value.size() == p.size());

		for (unsigned int i = 0 ; i < value.size() ; i++)
		{
			double norm = sqrt(p[i].square());

			value[i] = 0;
			value[i][0] = A0 + A2*norm*norm + A1*p[i][0]/norm;	
		}

	}

	template<int dim> string generate_systemA<dim>::generate_filename_to_write(const string folder_name,const string filename)
	{
			assert(sub_directory_names[4].compare("outputs/matrix_read") == 0);

			string filename_out;
			filename_out = sub_directory_names[4] + "/" + to_string(nEqn) + filename.substr(folder_name.size());

			return filename_out;
	}

	template<int dim> void generate_systemA<dim>::build_P(system_matrix &matrix_info)
	{
			assert(matrix_info.matrix.cols() !=0 || matrix_info.matrix.rows() !=0);

			matrix_info.matrix.coeffRef(1,1) = zeta/tau;
			for (unsigned int i = 2 ; i < matrix_info.matrix.rows() ; i++)
				matrix_info.matrix.coeffRef(i,i) = 1/tau;
	}

	template<int dim> void generate_systemA<dim>::build_BC(system_matrix &matrix_info)
	{
			assert(matrix_info.matrix.cols() !=0 || matrix_info.matrix.rows() !=0);
  		
  			matrix_info.matrix.coeffRef(0,0) = 1.0;
  			matrix_info.matrix.coeffRef(1,0) = chi;
  			matrix_info.matrix.coeffRef(1,3) = 1.0*chi;
  			matrix_info.matrix.coeffRef(2,2) = 1.0;
  			matrix_info.matrix.coeffRef(3,3) = 1.0;
  			matrix_info.matrix.coeffRef(4,2) = chi;
		    matrix_info.matrix.coeffRef(5,5) = 1.0;			
	}
	
	template<int dim> void generate_systemA<dim>::build_BCrhs(const Point<dim> p,const Point<dim> normal_vector,Vector<double> &bc_rhs) const
	{
		assert(bc_rhs.size() == nEqn);
		double norm = sqrt(p.square());

		  if( norm > 0.7 ) {
 		   bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta
  		  } else {
    		bc_rhs(1) = chi*theta0;
    		bc_rhs(4) = uW*normal_vector(1);
  		  }; 
 
	}

	template<int dim> Point<dim> generate_systemA<dim>::mirror(const Point<dim> normal_vector)const
	{
		Point<dim> mirrored_vector(normal_vector(0),-normal_vector(1));
		return mirrored_vector;
	}

	template<int dim> Sparse_matrix generate_systemA<dim>::build_Projector(const Point<dim> normal_vector) const
	{
		Sparse_matrix Projector;
		Projector.resize(nEqn,nEqn);

		double nx = normal_vector(0);
		double ny = normal_vector(1);
		assert(Projector.rows() == nEqn || Projector.cols() == nEqn);

		  Projector.coeffRef(0,0) = 1.0;
 		  Projector.coeffRef(1,1) = nx;
  		  Projector.coeffRef(1,2) = ny;
  	      Projector.coeffRef(2,1) = -ny;
  		  Projector.coeffRef(2,2) = nx;
  		  Projector.coeffRef(3,3) = nx*nx;
  		  Projector.coeffRef(3,4) = 2*nx*ny;
  		  Projector.coeffRef(3,5) = ny*ny;
  		  Projector.coeffRef(4,3) = -nx*ny;
  		  Projector.coeffRef(4,4) = nx*nx-ny*ny;
  		  Projector.coeffRef(4,5) = nx*ny;
  		  Projector.coeffRef(5,3) = ny*ny;
  		  Projector.coeffRef(5,4) = -2*nx*ny;
  		  Projector.coeffRef(5,5) = nx*nx;
  		return Projector;
 
	}

	template<int dim> Sparse_matrix generate_systemA<dim>::build_InvProjector(const Point<dim> normal_vector) const
	{
		return( build_Projector( mirror(normal_vector) ) );
	}

	template<int dim> Full_matrix generate_systemA<dim>::build_Aminus( const Point<dim> normal_vector) const
	{
		MatrixXd mat;
		mat.resize(nEqn,nEqn);
		double nx = normal_vector(0);
		double ny = normal_vector(1);

		mat << sqrt(0.6),-nx,-ny,-(sqrt(0.6)*(-1 + ny*ny)),2*sqrt(0.6)*nx*ny,sqrt(0.6)*ny*ny,
		-nx,sqrt(5./3) + (-sqrt(5./3) + 1./sqrt(2))*ny*ny,nx*(sqrt(5./3)*ny - ny/sqrt(2)),-nx,-ny,0,
		-ny,nx*(sqrt(5./3)*ny - ny/sqrt(2)),1./sqrt(2) + (sqrt(5./3) - 1./sqrt(2))*ny*ny,0,-nx,-ny,
		(2 - 3*ny*ny)/sqrt(15),(-2*nx)/3.,ny/3.,((-1 + ny*ny)*(-2*sqrt(15) + 3*(-5*sqrt(2) + sqrt(15))*ny*ny))/15.,nx*ny*(-sqrt(2) + 4./sqrt(15) + (-2*sqrt(0.6) + 2*sqrt(2))*ny*ny),ny*ny*(-sqrt(2) + 2./sqrt(15) + (-sqrt(0.6) + sqrt(2))*ny*ny),
		sqrt(0.6)*nx*ny,-ny/2.,-nx/2.,nx*ny*(sqrt(0.6) - 1./sqrt(2) + (-sqrt(0.6) + sqrt(2))*ny*ny),1./sqrt(2) + (2*(-5*sqrt(2) + sqrt(15))*ny*ny)/5. + (-2*sqrt(0.6) + 2*sqrt(2))*ny*ny*ny*ny,nx*(ny/sqrt(2) + (sqrt(0.6) - sqrt(2))*ny*ny*ny),    
		(-1 + 3*ny*ny)/sqrt(15),nx/3.,(-2*ny)/3.,-((-1 + ny*ny)*(-sqrt(15) + 3*(-5*sqrt(2) + sqrt(15))*ny*ny))/15.,(nx*ny*(15*sqrt(2) - 2*sqrt(15) + 6*(-5*sqrt(2) + sqrt(15))*ny*ny))/15.,ny*ny*(sqrt(2) - 1./sqrt(15) + (sqrt(0.6) - sqrt(2))*ny*ny);

		return( mat );
	};
	template<int dim> generate_systemA<dim>::generate_systemA()
	{
		string system_dir;
		string filename;
		string filename_out;

		system_dir = "system_matrices/";

		cout << "Reading system matrices.......\n" << endl;
		for (unsigned int i = 0 ; i < dim ; i ++)
		{
			filename = system_dir + "A" + to_string(i+1) + ".txt";
			this->build_triplet(A[i],filename);

			A[i].matrix.resize(nEqn,nEqn);
			this->build_matrix_from_triplet(A[i]);
			this->print_matrix(A[i],generate_filename_to_write(system_dir,filename));
		}

		filename = system_dir + "P.txt";
		P.matrix.resize(nEqn,nEqn);
		this->build_triplet(P,filename);
		this->build_P(P);
		this->print_matrix(P,generate_filename_to_write(system_dir,filename));

		filename = system_dir + "BC.txt";
		BC.matrix.resize(nEqn,nEqn);
		this->build_BC(BC);
		this->print_matrix(BC,generate_filename_to_write(system_dir,filename));

		cout << "Done Reading Matrices......\n" << endl;
	
	}
}
