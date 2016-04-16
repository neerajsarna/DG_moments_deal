	template<int dim> class generate_systemB:public Base_EquationGenerator<dim>,protected constants_SystemB
	{
		public:
			enum Num_flux
          	{Upwind, LLF};

          	Num_flux num_flux;
			generate_systemB(const enum Num_flux num_flux);
			system_matrix A[dim];
			system_matrix P;
			system_matrix BC;

			virtual void build_BCrhs(const Point<dim> p,const Point<dim> normal_vector,Vector<double> &bc_rhs) const;
			virtual Sparse_matrix build_Projector(const Point<dim> normal_vector) const;
			virtual Sparse_matrix build_InvProjector(const Point<dim> normal_vector) const;
			virtual Full_matrix build_Aminus(const Point<dim> normal_vector) const;
			Full_matrix build_Aminus_boundary(const Point<dim> normal_vector) const;
			virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) const;		

	  protected:
	  		Full_matrix Aminus_1D_Int;
	  		Full_matrix Aminus_1D_Bound;
	  		virtual void build_Aminus1D();
	  		string generate_filename_to_write(const string folder_name,const string filename);
			virtual void build_P(system_matrix &P);
			virtual void build_BC(system_matrix &BC);
			virtual Point<dim> mirror(const Point<dim> normal_vector) const ;

			/*virtual void build_Aminus1D(const Sparse_matrix Ax) const; */
	};

		template<int dim> void generate_systemB<dim>::build_Aminus1D()
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

	template<int dim> void generate_systemB<dim>::source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) const
	{
		assert(value.size() == p.size());

		for (unsigned int i = 0 ; i < value.size() ; i++)
		{
			double r = sqrt(p[i].square());

			value[i] = 0;
			value[i][0] = A0 + A2*r*r + A1*(1.0-5.0/18*r*r/(tau*tau))*p[i](0)/r;
		}

	}

	// for the present system, BC will be taken up from Mathematica assuming that chi is 1.0
	template<int dim> void generate_systemB<dim>::build_BC(system_matrix &BC)
	{
		;
	}

	template<int dim> string generate_systemB<dim>::generate_filename_to_write(const string folder_name,const string filename)
	{
		assert(sub_directory_names[4].compare("outputs/matrix_read") == 0);

		string filename_out;
		filename_out = sub_directory_names[4] + "/" +filename.substr(folder_name.size());

		return filename_out;
	}

		template<int dim> void generate_systemB<dim>::build_P(system_matrix &matrix_info)
	{
			assert(matrix_info.matrix.cols() !=0 || matrix_info.matrix.rows() !=0);

			matrix_info.matrix.coeffRef(1,1) = zeta/tau;
			for (unsigned int i = 2 ; i < matrix_info.matrix.rows() ; i++)
				matrix_info.matrix.coeffRef(i,i) = 1/tau;
	}

	template<int dim> void generate_systemB<dim>::build_BCrhs(const Point<dim> p,const Point<dim> normal_vector,Vector<double> &bc_rhs) const
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

	template<int dim> Point<dim> generate_systemB<dim>::mirror(const Point<dim> normal_vector)const
	{
		Point<dim> mirrored_vector(normal_vector(0),-normal_vector(1));
		return mirrored_vector;
	}

	template<int dim> Sparse_matrix generate_systemB<dim>::build_Projector( const Point<dim> normal_vector) const
	{
		double nx = normal_vector(0), ny = normal_vector(1);
		double nxnx = nx*nx, nyny = ny*ny;
		Sparse_matrix T(nEqn,nEqn);

		T.coeffRef(0,0) = 1.0;
		T.coeffRef(1,1) = nx;
		T.coeffRef(1,2) = ny;
		T.coeffRef(2,1) = -ny;
		T.coeffRef(2,2) = nx;
		T.coeffRef(3,3) = nxnx;
		T.coeffRef(3,4) = 2*nx*ny;
		T.coeffRef(3,5) = nyny;
		T.coeffRef(4,3) = -nx*ny;
		T.coeffRef(4,4) = nxnx-nyny;
		T.coeffRef(4,5) = nx*ny;
		T.coeffRef(5,3) = nyny;
		T.coeffRef(5,4) = -2*nx*ny;
		T.coeffRef(5,5) = nxnx;
		T.coeffRef(6,6) = nx*nxnx;
		T.coeffRef(6,7) = 3*ny*nxnx;
		T.coeffRef(6,8) = 3*nx*nyny;
		T.coeffRef(6,9) = ny*nyny;
		T.coeffRef(7,6) = -ny*nxnx;
		T.coeffRef(7,7) = nx*nxnx - 2*nx*nyny;
		T.coeffRef(7,8) = 2*ny*nxnx - ny*nyny;
		T.coeffRef(7,9) = nx*nyny;
		T.coeffRef(8,6) = nx*nyny;
		T.coeffRef(8,7) = -2*ny*nxnx + ny*nyny;
		T.coeffRef(8,8) = nx*nxnx - 2*nx*nyny;
		T.coeffRef(8,9) = ny*nxnx;
		T.coeffRef(9,6) = -ny*nyny;
		T.coeffRef(9,7) = 3*nx*nyny;
		T.coeffRef(9,8) = -3*ny*nxnx;
		T.coeffRef(9,9) = nx*nxnx;

		T.makeCompressed();
		return( T );
	};

	template<int dim> Sparse_matrix generate_systemB<dim>::build_InvProjector(const Point<dim> normal_vector) const
	{
		return( build_Projector( mirror(normal_vector) ) );
	}

	template<int dim> Full_matrix generate_systemB<dim>::build_Aminus( const Point<dim> normal_vector) const
	{
		return( build_InvProjector(normal_vector) * Aminus_1D_Int * build_Projector(normal_vector) );
	};

	template<int dim> Full_matrix generate_systemB<dim>::build_Aminus_boundary(const Point<dim> normal_vector) const
	{
		return( build_InvProjector(normal_vector) * Aminus_1D_Bound * build_Projector(normal_vector) );
	};


	template<int dim> generate_systemB<dim>::generate_systemB(const enum Num_flux num_flux)
	:
	num_flux(num_flux)
	{
		cout << "Total Equations: " << nEqn << endl;		
		string system_dir;
		string filename;
		string filename_out;

		system_dir = "system_matrices/";

		cout << "CAUTION: system running on chi = 1" << endl;
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
		this->build_triplet(BC,filename);
		this->build_matrix_from_triplet(BC);
		this->print_matrix(BC,generate_filename_to_write(system_dir,filename));

		Aminus_1D_Int.resize(nEqn,nEqn);
		Aminus_1D_Bound.resize(nEqn,nEqn);
		build_Aminus1D();

		cout << "Done Reading Matrices......\n" << endl;
	
	}



