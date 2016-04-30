template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>::build_BC(system_matrix &matrix_info,const unsigned int system_id)
{
	assert(matrix_info.matrix.cols() == num_equations.total_nEqn[system_id] 
				|| matrix_info.matrix.rows() == num_equations.total_nEqn[system_id]);

	switch(system_id)
	{
		case 0:
		{
			assert(num_equations.total_nEqn[system_id] == 6);

			matrix_info.matrix.coeffRef(0,0) = 1.0;
  			matrix_info.matrix.coeffRef(1,0) = chi;
  			matrix_info.matrix.coeffRef(1,3) = 1.0*chi;
  			matrix_info.matrix.coeffRef(2,2) = 1.0;
  			matrix_info.matrix.coeffRef(3,3) = 1.0;
  			matrix_info.matrix.coeffRef(4,2) = chi;
		    matrix_info.matrix.coeffRef(5,5) = 1.0;			
			break;
		}

		case 1:
		{
			assert(num_equations.total_nEqn[system_id] == 10);

			matrix_info.matrix.coeffRef(0,0) = 1.0;
  			matrix_info.matrix.coeffRef(1,0) = chi;
  			matrix_info.matrix.coeffRef(2,2) = 1.0;
  			matrix_info.matrix.coeffRef(3,3) = 1.0;
  			matrix_info.matrix.coeffRef(4,2) = chi;
		    matrix_info.matrix.coeffRef(4,7) = 2 * chi;			
		    matrix_info.matrix.coeffRef(5,5) = 1.0;
		    matrix_info.matrix.coeffRef(6,3) = chi;
		    matrix_info.matrix.coeffRef(7,7) = 1;
		    matrix_info.matrix.coeffRef(8,5) = chi;
		    matrix_info.matrix.coeffRef(9,9) = 1.0;
			break;
		}

	}

}

template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>
::build_P(system_matrix &matrix_info,const unsigned int system_id)
{
	assert(matrix_info.matrix.cols() == num_equations.total_nEqn[system_id] 
			|| matrix_info.matrix.rows() == num_equations.total_nEqn[system_id]);

	matrix_info.matrix.coeffRef(1,1) = zeta/tau;

	for (unsigned int i = 2 ; i < matrix_info.matrix.rows() ; i++)
				matrix_info.matrix.coeffRef(i,i) = 1/tau;

}

template<int system_type,int num_flux,int dim> 
Tensor<1,dim,double> Base_EquationGenerator<system_type,num_flux,dim>
::mirror(const Tensor<1,dim,double> normal_vector) const
{
		double nx = normal_vector[0], ny = normal_vector[1];
		Tensor<1,dim,double> mirrored_vector;
		mirrored_vector[0] = nx;
		mirrored_vector[1] = -ny;

		return mirrored_vector;
}

template<int system_type,int num_flux,int dim> 
Sparse_matrix Base_EquationGenerator<system_type,num_flux,dim>::
build_Projector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id) const
{
		Sparse_matrix Projector;
		Projector.resize(system_data[system_id].nEqn,system_data[system_id].nEqn);

		double nx = normal_vector[0];
		double ny = normal_vector[1];
		double nxnx = nx * nx;
		double nyny = ny * ny;
		assert(Projector.rows() == system_data[system_id].nEqn 
				|| Projector.cols() == system_data[system_id].nEqn);

	switch(system_id)
	{

		// Details for system-A
		case 0:
		{
			assert(num_equations.total_nEqn[system_id] == 6);

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
			break;
		}

		// Details for system-B
		case 1:
		{
		assert(num_equations.total_nEqn[system_id] == 10);

		Projector.coeffRef(0,0) = 1.0;
		Projector.coeffRef(1,1) = nx;
		Projector.coeffRef(1,2) = ny;
		Projector.coeffRef(2,1) = -ny;
		Projector.coeffRef(2,2) = nx;
		Projector.coeffRef(3,3) = nxnx;
		Projector.coeffRef(3,4) = 2*nx*ny;
		Projector.coeffRef(3,5) = nyny;
		Projector.coeffRef(4,3) = -nx*ny;
		Projector.coeffRef(4,4) = nxnx-nyny;
		Projector.coeffRef(4,5) = nx*ny;
		Projector.coeffRef(5,3) = nyny;
		Projector.coeffRef(5,4) = -2*nx*ny;
		Projector.coeffRef(5,5) = nxnx;
		Projector.coeffRef(6,6) = nx*nxnx;
		Projector.coeffRef(6,7) = 3*ny*nxnx;
		Projector.coeffRef(6,8) = 3*nx*nyny;
		Projector.coeffRef(6,9) = ny*nyny;
		Projector.coeffRef(7,6) = -ny*nxnx;
		Projector.coeffRef(7,7) = nx*nxnx - 2*nx*nyny;
		Projector.coeffRef(7,8) = 2*ny*nxnx - ny*nyny;
		Projector.coeffRef(7,9) = nx*nyny;
		Projector.coeffRef(8,6) = nx*nyny;
		Projector.coeffRef(8,7) = -2*ny*nxnx + ny*nyny;
		Projector.coeffRef(8,8) = nx*nxnx - 2*nx*nyny;
		Projector.coeffRef(8,9) = ny*nxnx;
		Projector.coeffRef(9,6) = -ny*nyny;
		Projector.coeffRef(9,7) = 3*nx*nyny;
		Projector.coeffRef(9,8) = -3*ny*nxnx;
		Projector.coeffRef(9,9) = nx*nxnx;
		return Projector;
		break;
		}
	}

	assert( 1 == 0);
	return Projector;
}

template<int system_type,int num_flux,int dim> 
Sparse_matrix Base_EquationGenerator<system_type, num_flux,dim>::
build_InvProjector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id) const
{
			return( build_Projector( mirror(normal_vector) ,system_id) );
}

template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>::
build_BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int system_id) const
{
		assert(bc_rhs.size() == num_equations.total_nEqn[system_id]);
		double norm = p.norm();

		switch(system_id)
		{
			case 0:
			{
				assert(num_equations.total_nEqn[system_id] == 6);

				if( norm > 0.7 ) 
 		   			bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta
		  		else 
		  		{
    				bc_rhs(1) = chi*theta0;
    				bc_rhs(4) = uW*normal_vector[1];
  		  		}; 

				break;
			}

			case 1:
			{
				assert(num_equations.total_nEqn[system_id] == 10);

				if( norm > 0.7 ) 
 		   			bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta
		  		else 
		  		{
    				bc_rhs(1) = chi*theta0;
    				bc_rhs(4) = uW*normal_vector[1];
  		  		}; 

				break;
			}
		}
}

template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>::
source_term(const vector<Point<dim>> &p,
			vector<Vector<double>> &value,
			const unsigned int system_id)
{
	assert(value.size() == p.size());
	switch(system_id)
	{
		case 0:
		{
			assert(num_equations.total_nEqn[system_id] == 6);

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
						Vector<double> force_value(num_equations.total_nEqn[system_id]);
						double norm = sqrt(p[i].square());
						force_value[0] = (A0 + A2*norm*norm + A1*p[i][0]/norm);		
						Sparse_matrix_dot_Vector(system_data[system_id].S, force_value,value[i]);
					}

					break;
				}
			}
			break;
		}

		case 1:
		{
			assert(num_equations.total_nEqn[system_id] == 10);

			switch(system_type)
			{
				case un_symmetric:
				{
					for (unsigned int i = 0 ; i < value.size() ; i++)
					{
						double r = sqrt(p[i].square());
						value[i] = 0;
						value[i][0] = A0 + A2 * r * r + A1*(1.0-5.0/18*r*r/(tau*tau))*p[i][0] / r;
					}

					break;
				}

				case symmetric:
				{	
					assert(1 == 0);
					for (unsigned int i = 0 ; i < value.size() ; i ++)
					{

						Vector<double> force_value(num_equations.total_nEqn[system_id]);
						double r = sqrt(p[i].square());

						force_value[0] = A0 + A2 * r * r + A1*(1.0-5.0/18*r*r/(tau*tau))*p[i][0]/r;
						Sparse_matrix_dot_Vector(system_data[system_id].S, force_value,value[i]);
					}

					break;
				}
			}
			break;
		}
	}
}

template<int system_type,int num_flux,int dim> 
Full_matrix Base_EquationGenerator<system_type,num_flux,dim>
::build_Aminus(const Tensor<1,dim,double> normal_vector,const unsigned int system_id)
{
		switch(system_type)
		{
			case un_symmetric:
			{
				return( build_InvProjector(normal_vector,system_id) 
						* system_data[system_id].Aminus_1D_Int * build_Projector(normal_vector,system_id) );		
				break;
			}

			case symmetric:
			{
				assert(num_equations.total_nEqn[system_id] == 10);
				Eigen::MatrixXd SAn = system_data[system_id].S.matrix * build_InvProjector(normal_vector,system_id) 
										* system_data[system_id].Ax.matrix *  build_Projector(normal_vector,system_id);

				EigenSolver<MatrixXd> ES(SAn);
				MatrixXd vecs = ES.pseudoEigenvectors();
				VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

				Eigen::MatrixXd Aminus_SAn = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();
				return(Aminus_SAn);
				break;
			}
		}

}

template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>
::build_Aminus1D(Full_matrix &Aminus_1D_Int,
				Full_matrix &Aminus_1D_Bound,
				const unsigned int system_id)
{
	assert(Aminus_1D_Bound.rows() == num_equations.total_nEqn[system_id] 
			&& Aminus_1D_Bound.cols() == num_equations.total_nEqn[system_id] );

	assert(Aminus_1D_Int.rows() == num_equations.total_nEqn[system_id]  
			&& Aminus_1D_Int.cols() == num_equations.total_nEqn[system_id] );

	Aminus_1D_Bound.setZero();
	Aminus_1D_Int.setZero();

	EigenSolver<MatrixXd> ES(system_data[system_id].A[0].matrix);
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

template<int system_type,int num_flux,int dim> 
void Base_EquationGenerator<system_type,num_flux,dim>::generate_matrices(equation_data &system_data,
																		const unsigned int system_id)
{
		if (system_type == 0)
			system_data.is_symmetric = true;
		else
			system_data.is_symmetric = false;

	    system_data.nEqn = num_equations.total_nEqn[system_id];

		cout <<"System ID " << system_id 
			 << " #nEqn " << system_data.nEqn << endl;

		string system_dir;
		string filename;
		string filename_out;

		system_dir = "system_matrices/";

		cout << "Reading system matrices.......\n" << endl;
		for (unsigned int i = 0 ; i < dim ; i ++)
		{
			filename = system_dir + to_string(system_data.nEqn) + "A" + to_string(i+1) + ".txt";
			build_triplet(system_data.A[i],filename);

			system_data.A[i].matrix.resize(system_data.nEqn,system_data.nEqn);
			build_matrix_from_triplet(system_data.A[i]);
			print_matrix(system_data.A[i],generate_filename_to_write(system_dir,filename));
		}

		filename = system_dir + to_string(system_data.nEqn) + "P.txt";
		system_data.P.matrix.resize(system_data.nEqn,system_data.nEqn);
		build_triplet(system_data.P,filename);
		build_P(system_data.P,system_id);
		print_matrix(system_data.P,generate_filename_to_write(system_dir,filename));

		filename = system_dir + to_string(system_data.nEqn) + "BC.txt";
		system_data.BC.matrix.resize(system_data.nEqn,system_data.nEqn);
		build_BC(system_data.BC,system_id);
		print_matrix(system_data.BC,generate_filename_to_write(system_dir,filename));

		system_data.Aminus_1D_Int.resize(system_data.nEqn,system_data.nEqn);
		system_data.Aminus_1D_Bound.resize(system_data.nEqn,system_data.nEqn);
		build_Aminus1D(system_data.Aminus_1D_Int,
						system_data.Aminus_1D_Bound,
						system_id);

		switch(system_type)
		{
			case symmetric:
			{
				assert(num_equations.total_nEqn[system_id] == 10); // symmetric system not implemented for B type

				filename = system_dir + to_string(system_data.nEqn) + "S.txt";
				system_data.S.matrix.resize(system_data.nEqn,system_data.nEqn);
				build_triplet(system_data.S,filename);
				build_matrix_from_triplet(system_data.S);
				print_matrix(system_data.S,generate_filename_to_write(system_dir,filename));

				system_data.Ax.matrix = system_data.A[0].matrix;
				for (unsigned int i = 0 ; i < dim ; i ++)
					system_data.A[i].matrix = system_data.S.matrix * system_data.A[i].matrix;

				system_data.P.matrix = system_data.S.matrix * system_data.P.matrix;
				break;
			}

			case un_symmetric:
			{
				break;
			}

		}


		cout << "Done Reading Matrices......\n" << endl;
}