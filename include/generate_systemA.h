	template<int dim> class generate_systemA:public Base_EquationGenerator<dim>
	{
		public:
			enum Num_flux
          	{Upwind, LLF};

          	enum System_type
          	{symmetric, un_symmetric};

          	// constants of the system
          	double tau;
          	double zeta;
          	double chi;
          	double theta0;
          	double theta1;
          	double uW;
          	double A0;
          	double A1;
          	double A2;
          	unsigned int nEqn;

          	Num_flux num_flux;
          	System_type system_type;
			generate_systemA(const enum Num_flux num_flux,const enum System_type system_type);
			system_matrix A[dim];
			system_matrix P;
			system_matrix BC;
			system_matrix S;
			system_matrix Ax;

			virtual void build_BCrhs(const Tensor<1,dim,double> p,const Tensor<1,dim,double> normal_vector,Vector<double> &bc_rhs) ;
			virtual Sparse_matrix build_Projector(const Tensor<1,dim,double> normal_vector) ;
			virtual Sparse_matrix build_InvProjector(const Tensor<1,dim,double> normal_vector) ;
			virtual Full_matrix build_Aminus(const Tensor<1,dim,double> normal_vector) ;
			Full_matrix build_Aminus_boundary(const Tensor<1,dim,double> normal_vector) ;
			virtual void source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) ;		

			struct exact_solution:public Function<dim>
			{
				exact_solution();
				double s_r(const double ,const double ) const = 0;
				double s_phi(const double r,const double ) const = 0;
				double thetaP(const double r,const double) const = 0;			
				void vector_value(const Point<dim> &p,Vector<double> &value) const = 0;
				double C1, C2, C3, C4, K1, K2, lambda1;
			};

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

		Aminus_1D_Bound.setZero();
		Aminus_1D_Int.setZero();

		EigenSolver<MatrixXd> ES(A[0].matrix);
  		MatrixXd vecs = ES.pseudoEigenvectors();
  		VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
 		double maxEV = vals.cwiseAbs().maxCoeff();

 		


		switch (num_flux)
		{
			case Upwind:
			{
				Aminus_1D_Bound = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();	 
  			  	Aminus_1D_Int = Aminus_1D_Bound;
				break;
			}

			case LLF:
			{
				 for (unsigned int i = 0 ; i < Aminus_1D_Bound.rows() ; i++)
					Aminus_1D_Bound(i,i) = fabs(maxEV);

				Aminus_1D_Int = Aminus_1D_Bound;

				break;

			}
		}
	}

	template<int dim> void generate_systemA<dim>::source_term(const vector<Point<dim>> &p,vector<Vector<double>> &value) 
	{
		assert(value.size() == p.size());

		switch(system_type)
		{
			case un_symmetric:
			{
				for (unsigned int i = 0 ; i < value.size() ; i++)
				{
					double norm = sqrt(p[i].square());

					value[i] = 0;
					value[i][0] = (A0 + A2*norm*norm + A1*p[i][0]/norm);	
				}

				break;
			}

			case symmetric:
			{	
				for (unsigned int i = 0 ; i < value.size() ; i ++)
				{
					Vector<double> force_value(this->nEqn);
					double norm = sqrt(p[i].square());

					force_value[0] = (A0 + A2*norm*norm + A1*p[i][0]/norm);		
					this->Sparse_matrix_dot_Vector(S, force_value,value[i]);
				}

				break;
			}
		}



	}

	template<int dim> string generate_systemA<dim>::generate_filename_to_write(const string folder_name,const string filename)
	{
			assert(this->sub_directory_names[4].compare("outputs/matrix_read") == 0);

			string filename_out;
			filename_out = this->sub_directory_names[4] + "/" + filename.substr(folder_name.size());

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

		switch(system_type)
		{
			case un_symmetric:
			{
				return( build_InvProjector(normal_vector) * Aminus_1D_Int * build_Projector(normal_vector) );		
				break;
			}

			case symmetric:
			{
				Eigen::MatrixXd SAn = S.matrix * build_InvProjector(normal_vector) * Ax.matrix *  build_Projector(normal_vector);

				EigenSolver<MatrixXd> ES(SAn);
				MatrixXd vecs = ES.pseudoEigenvectors();
				VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

				Eigen::MatrixXd Aminus_SAn = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();
				return(Aminus_SAn);
				break;
			}
		}
		
	}

	template<int dim> Full_matrix generate_systemA<dim>::build_Aminus_boundary(const Tensor<1,dim,double> normal_vector) 
	{
		return( build_InvProjector(normal_vector) * Aminus_1D_Bound * build_Projector(normal_vector) );
	}

	template<int dim> generate_systemA<dim>::generate_systemA(const enum Num_flux num_flux,const enum System_type system_type)
	:
	num_flux(num_flux),
	system_type(system_type)
	{
		nEqn = 6;

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

		filename = system_dir + to_string(nEqn) + "S.txt";
		S.matrix.resize(nEqn,nEqn);
		this->build_triplet(S,filename);
		this->build_matrix_from_triplet(S);
		this->print_matrix(S,generate_filename_to_write(system_dir,filename));

		switch(system_type)
		{
			case symmetric:
			{
				Ax.matrix = A[0].matrix;
				for (unsigned int i = 0 ; i < dim ; i ++)
					A[i].matrix = S.matrix * A[i].matrix;

				P.matrix = S.matrix * P.matrix;
				break;
			}

			case un_symmetric:
			{
				break;
			}

		}

		cout << "Done Reading Matrices......\n" << endl;
	
	}


	template<int dim> class generate_systemA<dim>::exact_solution::exact_solution()
	:
	Function<dim>(nEqn)
	{		

		if( fabs(tau - 0.01) < 1e-5 )
		{
			C1 = -258.6317250754188;
			C2 = 0.1403951449098079;
			C3 = -23.22459252429678;
			C4 = 6.490985249358112;
			K1 = 9.713664158637732e28;
			K2 = 1.0821414602962837e-127;	
		}


		if( fabs(tau - 0.1) < 1e-5 )
		{
			C1 = -1.9505992348859245;
			C2 = 0.143798955700584;
			C3 = 1.8448329709116922;
			C4 = 0.11342119104343287;
			K1 = 53.97166647969362;
			K2 = 1.26255082813026e-15;

		}

		if( fabs(tau - 1.0) < 1e-5 )
		{
			C1 = -0.14869287988905205;
			C2 = 0.17257205989623575;
			C3 = 1.2976971448738395;
			C4 = -0.10358631917431671;
			K1 = 0.18803597175199138;
			K2 = 0.003604064889912796;

		}

		if( fabs(tau - 10.0) < 1e-5 )
		{
			C1 = -0.1885054386029859;
			C2 = 1.4351918021177739;
			C3 = 0.11716270300600619;
			C4 = -0.004296145849336291;
			K1 = 0.18900239036185373;
			K2 = -0.6486989954482233;

		}

		lambda1 = sqrt(2 * zeta);
	}

template<int dim> double generate_systemA<dim>::exact_solution::thetaP(const double r,const double phi) const
	{
		return C3 + C4 * zeta * log(r) - (A0 * tau * zeta * pow(r,2))/4. + 
		cos(phi)*(C2 * r * zeta + C1 * zeta * pow(r,-1) - (A1 * tau * (-2 + zeta * pow(r,2)))/3.) + 
		(2 * A2 * pow(r,2) * pow(tau,3))/3. - (A2 * zeta * pow(r,4) * pow(tau,3))/16.;
	}

template<int dim> double generate_systemA<dim>::exact_solution::s_r(const double r,const double phi) const
	{
		return (A0 * r * tau)/2. - C4*pow(r,-1) + cos(phi) * (-C2 + (2 * A1 * r * tau)/3. + C1*pow(r,-2) - 
			K2*this->BI(1,r*lambda1)*pow(2,0.5)*pow(r,-1) + K1*this->BK(1,r*lambda1)*pow(2,0.5)*pow(r,-1))
		+ (A2 * pow(r,3) * pow(tau,3))/4.;

	}

template<int dim> double generate_systemA<dim>::exact_solution::s_phi(const double r,const double phi) const
	{
		return (C2 - (A1 * r * tau)/3. + C1*pow(r,-2) + 
			K2*(this->BI(0,lambda1 * r) * pow(zeta,0.5) + this->BI(2,lambda1 * r) * pow(zeta,0.5)) + 
			K1*(this->BK(0,lambda1 * r) * pow(zeta,0.5) + this->BK(2,lambda1 * r) * pow(zeta,0.5)))*sin(phi);

	}


 template<int dim> void generate_systemA<dim>::exact_solution::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		assert(value.size() != 0);

		double r = sqrt(p.square());
		double phi = atan2(p[1],p[0]);
		r /=tau;



		value[0] =  thetaP(r,phi);
											// theta
		value[1] = cos(phi) * s_r(r,phi) - sin(phi) * s_phi(r,phi);
								// qx
		value[2] = sin(phi) * s_r(r,phi) + cos(phi) * s_phi(r,phi);								// qy

	// allocating zero value for higher order moments(values not needed during error evaluation)
		for (unsigned int i = 3 ; i < this->nEqn ; i++)
			value[i] = 0; // a value for higher order moments is not needed

   }


