	template<int dim> class generate_systemA:public Base_EquationGenerator<dim>,protected constants_SystemA
	{
		public:
			enum Num_flux
          	{Upwind, LLF};

          	Num_flux num_flux;
			generate_systemA(const enum Num_flux num_flux);
			system_matrix A[dim];
			system_matrix P;
			system_matrix BC;

			virtual void build_BCrhs(const Tensor<1,dim,double> p,const Tensor<1,dim,double> normal_vector,Vector<double> &bc_rhs) ;
			virtual Sparse_matrix build_Projector(const Tensor<1,dim,double> normal_vector) ;
			virtual Sparse_matrix build_InvProjector(const Tensor<1,dim,double> normal_vector) ;
			virtual Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector) ;
			Full_matrix build_Aminus_boundary(const Tensor<1,dim,double> normal_vector) ;
			virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) ;		

	  protected:
	  		Full_matrix Aminus_1D_Int;
	  		Full_matrix Aminus_1D_Bound;
	  		virtual void build_Aminus1D();
	  		string generate_filename_to_write(const string folder_name,const string filename);
			virtual void build_P(system_matrix &P);
			virtual void build_BC(system_matrix &BC);
			virtual Tensor<1,dim,double> mirror(const Tensor<1,dim,double> normal_vector) ;

			/*virtual void build_Aminus1D(const Sparse_matrix Ax) const; */
	};

	template<int dim> void generate_systemA<dim>::build_Aminus1D()
	{

		assert(Aminus_1D_Bound.rows() == nEqn && Aminus_1D_Bound.cols() == nEqn);
		assert(Aminus_1D_Int.rows() == nEqn && Aminus_1D_Int.cols() == nEqn);
		
		EigenSolver<MatrixXd> ES(A[0].matrix);
  		MatrixXd vecs = ES.pseudoEigenvectors();
  		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
 		double maxEV = vals.cwiseAbs().maxCoeff();

 		Aminus_1D_Bound = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();


		switch (num_flux)
		{
			case Upwind:
			{
			 
  			  	Aminus_1D_Int = Aminus_1D_Bound;
				break;
			}

			case LLF:
			{
				 for (unsigned int i = 0 ; i < Aminus_1D_Int.rows() ; i++)
					Aminus_1D_Int(i,i) = maxEV;

				break;

			}
		}
	}

	template<int dim> void generate_systemA<dim>::source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) 
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
			filename_out = sub_directory_names[4] + "/" + filename.substr(folder_name.size());

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
	
	template<int dim> void generate_systemA<dim>::build_BCrhs(const Tensor<1,dim,double> p,const Tensor<1,dim,double> normal_vector,Vector<double> &bc_rhs) 
	{
		assert(bc_rhs.size() == nEqn);
		double norm = p.norm();

		  if( norm > 0.7 ) {
 		   bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta
  		  } else {
    		bc_rhs(1) = chi*theta0;
    		bc_rhs(4) = uW*normal_vector[1];
  		  }; 
 
	}

	template<int dim> Tensor<1,dim,double> generate_systemA<dim>::mirror(const Tensor<1,dim,double> normal_vector)
	{
		double nx = normal_vector[0], ny = normal_vector[1];
		Tensor<1,dim,double> mirrored_vector;
		mirrored_vector[0] = nx;
		mirrored_vector[1] = -ny;
		return mirrored_vector;
	}

	template<int dim> Sparse_matrix generate_systemA<dim>::build_Projector(const Tensor<1,dim,double> normal_vector) 
	{
		Sparse_matrix Projector;
		Projector.resize(nEqn,nEqn);

		double nx = normal_vector[0];
		double ny = normal_vector[1];
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

	template<int dim> Sparse_matrix generate_systemA<dim>::build_InvProjector(const Tensor<1,dim,double> normal_vector) 
	{
		return( build_Projector( mirror(normal_vector) ) );
	}

	template<int dim> Full_matrix generate_systemA<dim>::build_Aminus( const Tensor<1,dim,double> normal_vector) 
	{
		return( build_InvProjector(normal_vector) * Aminus_1D_Int * build_Projector(normal_vector) );
	}

	template<int dim> Full_matrix generate_systemA<dim>::build_Aminus_boundary(const Tensor<1,dim,double> normal_vector) 
	{
		return( build_InvProjector(normal_vector) * Aminus_1D_Bound * build_Projector(normal_vector) );
	}

	template<int dim> generate_systemA<dim>::generate_systemA(const enum Num_flux num_flux)
	:
	num_flux(num_flux)
	{
		cout << "Total Equations: " << nEqn << endl;
		string system_dir;
		string filename;
		string filename_out;

		system_dir = "system_matrices/";

		cout << "Reading system matrices.......\n" << endl;
		for (unsigned int i = 0 ; i < dim ; i ++)
		{
			filename = system_dir + to_string(nEqn) + "A" + to_string(i+1) + ".txt";
			this->build_triplet(A[i],filename);

			A[i].matrix.resize(nEqn,nEqn);
			this->build_matrix_from_triplet(A[i]);
			this->print_matrix(A[i],generate_filename_to_write(system_dir,filename));
		}

		filename = system_dir + to_string(nEqn) + "P.txt";
		P.matrix.resize(nEqn,nEqn);
		this->build_triplet(P,filename);
		this->build_P(P);
		this->print_matrix(P,generate_filename_to_write(system_dir,filename));

		filename = system_dir + to_string(nEqn) + "BC.txt";
		BC.matrix.resize(nEqn,nEqn);
		this->build_BC(BC);
		this->print_matrix(BC,generate_filename_to_write(system_dir,filename));

		Aminus_1D_Int.resize(nEqn,nEqn);
		Aminus_1D_Bound.resize(nEqn,nEqn);
		build_Aminus1D();

		cout << "Done Reading Matrices......\n" << endl;
	
	}
