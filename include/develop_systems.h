// the following routine allocates memory for the projector
template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::init_tensor_data()
{
	
	const unsigned int no_of_tensors = tensor_info.free_index_info.size();
	tensor_project.resize(no_of_tensors);
	
	Assert(no_of_tensors == 4,ExcNotImplemented());
	Assert(tensor_project.size() != 0 , ExcNotInitialized());

	for (unsigned int i = 0 ; i < no_of_tensors ; i++)
	{
		const unsigned int n_free_indices = tensor_info.free_index_info[i];
		tensor_project[i].P.resize(n_free_indices,n_free_indices);
	}	
}


// build projectors corresponding to different degrees
template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_tensorial_projector(const double nx,const double ny)
{

	double nxnx = nx * nx;
	double nyny = ny * ny;

	tensor_project[0].P << 1;

	tensor_project[1].P << nx, ny, -ny, nx;	
	
	tensor_project[2].P << nxnx, 2*nx*ny, nyny, 
					    -nx*ny, nxnx-nyny, nx*ny, nyny, -2*nx*ny, nxnx;

		
	tensor_project[3].P << nx*nxnx, 3*ny*nxnx, 3*nx*nyny, ny*nyny, 
						  -ny*nxnx, nx*nxnx - 2*nx*nyny, 2*ny*nxnx - ny*nyny, nx*nyny, 
        				   nx*nyny, -2*ny*nxnx + ny*nyny, nx*nxnx - 2*nx*nyny, ny*nxnx,
        			       -ny*nyny, 3*nx*nyny, -3*ny*nxnx, nx*nxnx;


}

// put P in Sp at the diagonal location idx
template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::SpBlock(const unsigned int idx,const Full_matrix P,Sparse_matrix &Sp)
{
		for (int i = 0 ; i < P.rows() ; i++)
			for (int j = 0 ; j < P.cols() ; j++)
				Sp.coeffRef(idx + i,idx + j) = P(i,j);
}


template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_BC(system_matrix &matrix_info,const unsigned int system_id)
{
	assert(matrix_info.matrix.cols() == num_equations.total_nEqn[system_id] 
		   || matrix_info.matrix.rows() == num_equations.total_nEqn[system_id]);

	const unsigned int neqn_local = num_equations.total_nEqn[system_id];

	switch(neqn_local)
	{
		case 6:
		{
			assert(neqn_local == 6);

			matrix_info.matrix.coeffRef(0,0) = 1.0;
  			matrix_info.matrix.coeffRef(1,0) = chi;
  			matrix_info.matrix.coeffRef(1,3) = 1.0*chi;
  			matrix_info.matrix.coeffRef(2,2) = 1.0;
  			matrix_info.matrix.coeffRef(3,3) = 1.0;
  			matrix_info.matrix.coeffRef(4,2) = chi;
		    matrix_info.matrix.coeffRef(5,5) = 1.0;			
			break;
		}

		case 10:
		{
			assert(neqn_local == 10);

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

		// only implement characteristic boundary conditions for R13
		case 17:
		{
			Assert(1 == 0, ExcNotImplemented());
		}

	}

}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_B_tilde_inv(const system_matrix &B,
				const Full_matrix &X_minus,
				Full_matrix &B_tilde_inv)
{
	Assert(B_tilde_inv.size() != 0, ExcNotInitialized());
	B_tilde_inv = (B.matrix * X_minus).inverse();
}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_B_hat(const system_matrix &B,
			  const Full_matrix &X_minus,
			 const Full_matrix &B_tilde_inv,
			 Full_matrix &B_hat)
{
	Assert(B_hat.size() != 0,ExcNotInitialized());
	B_hat = X_minus * B_tilde_inv * B.matrix;

}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_P(system_matrix &matrix_info,const unsigned int system_id)
{
	assert(matrix_info.matrix.cols() == num_equations.total_nEqn[system_id] 
			|| matrix_info.matrix.rows() == num_equations.total_nEqn[system_id]);

	// depending upon the total number of equations we can describe the production terms
	const unsigned int neqn_local = num_equations.total_nEqn[system_id];

	// the G26A and the G26B system
	if (neqn_local == 6 || neqn_local == 10)
	{
		// only the first equations is a conservation law
		for (unsigned int i = 1 ; i < neqn_local ; i++)
				matrix_info.matrix.coeffRef(i,i) = 1.0/tau;
	}

	// all the other moment system will have dim + 2 conservation laws
	// all are BGK 
	for (unsigned int i = dim + 2 ; i  < neqn_local ; i++)
		matrix_info.matrix.coeffRef(i,i) = 1.0/tau;
		
}


// mirror of a 2d vector
template<int num_flux,int dim> 
Tensor<1,dim,double> Base_EquationGenerator<num_flux,dim>
::mirror(const Tensor<1,dim,double> normal_vector) const
{
		double nx = normal_vector[0], ny = normal_vector[1];
		Tensor<1,dim,double> mirrored_vector;
		mirrored_vector[0] = nx;
		mirrored_vector[1] = -ny;

		return mirrored_vector;
}

// builds the projector for the equations
template<int num_flux,int dim> 
Sparse_matrix Base_EquationGenerator<num_flux,dim>
::build_Projector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id) 
{
		Sparse_matrix Projector;
		Projector.resize(system_data[system_id].nEqn,system_data[system_id].nEqn);

		double nx = normal_vector[0];
		double ny = normal_vector[1];


		Assert(Projector.rows() == system_data[system_id].nEqn 
				|| Projector.cols() == system_data[system_id].nEqn,ExcNotInitialized());

		build_tensorial_projector(nx,ny);
		const unsigned int neqn_local = num_equations.total_nEqn[system_id];

	switch(neqn_local)
	{

		// Details for system-A
		// variables are: theta,qx, qy, Rxx, Rxy, Ryy
		case 6:
		{
			
			switch(system_type)
			{
				case un_symmetric:
				{
					SpBlock(0,tensor_project[0].P,Projector);
					SpBlock(1,tensor_project[1].P,Projector);
					SpBlock(3,tensor_project[2].P,Projector);				
					break;					
				}

				case symmetric:
				{
					Projector.coeffRef(0,0) = 1.0/sqrt(2);
					Projector.coeffRef(1,1) = nx/sqrt(2);
					Projector.coeffRef(1,2) = ny/sqrt(2);
					Projector.coeffRef(2,1) = -ny/sqrt(2);
					Projector.coeffRef(2,2) = nx/sqrt(2);
					Projector.coeffRef(3,3) = ((3 + sqrt(3))*nx*nx + (-3 + sqrt(3))*ny*ny)/(6 * sqrt(2));
					Projector.coeffRef(3,4) = nx*ny;
					Projector.coeffRef(3,5) = ((-3 + sqrt(3))*nx*nx + (3 + sqrt(3))*ny*ny)/(6 * sqrt(2));
					Projector.coeffRef(4,3) = -nx*ny/sqrt(2);
					Projector.coeffRef(4,4) = (nx*nx-ny*ny)/2;
					Projector.coeffRef(4,5) = nx*ny/sqrt(2);
					Projector.coeffRef(5,3) = ((-3 + sqrt(3))*nx*nx + (3 + sqrt(3))*ny*ny)/(6 * sqrt(2));
					Projector.coeffRef(5,4) = -nx*ny;
					Projector.coeffRef(5,5) = ((3 + sqrt(3))*nx*nx + (-3 + sqrt(3))*ny*ny)/(6 * sqrt(2));
					
					break;		

				}

			}
		
			break;
			

		}

		// Details for system-B
		// the variables are theta, qx, qy , Rxx,Rxy , Ryy , mxxx, mxxy, mxyy, myyy
		case 10:
		{
			
			switch(system_type)
			{
				case un_symmetric:
				{
					SpBlock(0,tensor_project[0].P,Projector);
					SpBlock(1,tensor_project[1].P,Projector);
					SpBlock(3,tensor_project[2].P,Projector);
					SpBlock(6,tensor_project[3].P,Projector);
					break;
				}

				case symmetric:
				{
					Assert(1 == 0,ExcNotImplemented());
					break;
				}

			}
			
			break;
		}

		// Details for the G26 system
		// rho, (vx, vy), theta, sigma_xx, sigma_xy, sigma_yy, qx, qy, (mxxx, mxxy, mxyy, myyy),Delta, (Rxx, Rxy, Ryy)
		case 17:
		{
			
			switch(system_type)
			{
				case un_symmetric:
				{
					  SpBlock( 0, tensor_project[0].P, Projector );	// rho
  					  SpBlock( 1, tensor_project[1].P, Projector );	// velocity
  					  SpBlock( 3, tensor_project[0].P, Projector ); // theta
  			          SpBlock( 4, tensor_project[2].P, Projector ); // sigma
  					  SpBlock( 7, tensor_project[1].P, Projector ); // q
  				      SpBlock( 9, tensor_project[3].P, Projector ); // mijk
  				      SpBlock( 13, tensor_project[0].P, Projector ); // delta
  				      SpBlock( 14, tensor_project[2].P, Projector ); // Rij

					break;
				}

				case symmetric:
				{
					Assert(1 == 0,ExcNotImplemented());
					break;
				}

			}
			
			break;
		}
	}

	return Projector;
}

// builds the inverse projector
template<int num_flux,int dim> 
Sparse_matrix Base_EquationGenerator< num_flux,dim>::
build_InvProjector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id)
{
	Sparse_matrix Inv_Projector;
	Inv_Projector.resize(system_data[system_id].nEqn,system_data[system_id].nEqn);

	double nx = normal_vector[0];
	double ny = normal_vector[1];
	double nxnx = nx * nx;
	double nyny = ny * ny;

	Assert(Inv_Projector.rows() == system_data[system_id].nEqn 
			|| Inv_Projector.cols() == system_data[system_id].nEqn,
			ExcNotInitialized());

	const unsigned int neqn_local = num_equations.total_nEqn[system_id];

	switch(neqn_local)
	{
		case 6:
		{
			switch(system_type)
			{
				case un_symmetric:
				{
					Inv_Projector = build_Projector(mirror(normal_vector),system_id);
					break;
				}

				case symmetric:
				{
					Inv_Projector.coeffRef(0,0) = sqrt(2);

					Inv_Projector.coeffRef(1,1) = sqrt(2) * nx;
					Inv_Projector.coeffRef(1,2) = -sqrt(2) * ny;
					Inv_Projector.coeffRef(2,1) = sqrt(2) * ny;
					Inv_Projector.coeffRef(2,2) = sqrt(2) * nx;

					Inv_Projector.coeffRef(3,3) = ((1 + sqrt(3))*nx*nx + (-1 + sqrt(3))*ny*ny)/sqrt(2);
					Inv_Projector.coeffRef(3,4) = -2 * sqrt(2) * nx*ny;
					Inv_Projector.coeffRef(3,5) = ((-1 + sqrt(3))*nx*nx + (1 + sqrt(3))*ny*ny)/sqrt(2);

					Inv_Projector.coeffRef(4,3) = 2 * nx * ny;
					Inv_Projector.coeffRef(4,4) = 2 * (nx*nx-ny*ny);
					Inv_Projector.coeffRef(4,5) = -2 * nx * ny;

					Inv_Projector.coeffRef(5,3) = ((-1 + sqrt(3))*nx*nx + (1 + sqrt(3))*ny*ny)/sqrt(2);
					Inv_Projector.coeffRef(5,4) = 2 * sqrt(2) * nx * ny ;
					Inv_Projector.coeffRef(5,5) = ((1 + sqrt(3))*nx*nx + (-1 + sqrt(3))*ny*ny)/sqrt(2);
					break;
				}
			}

			break;

		}
		case 10:
		case 17:
		{
			switch(system_type)
			{
				case un_symmetric:
				{
					Inv_Projector = build_Projector(mirror(normal_vector),system_id);
					break;
				}

				case symmetric:
				{
					Assert(1 == 0, ExcNotImplemented());
					break;
				}
			}
			break;
		}
	}

	return Inv_Projector;
			
}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>::
build_BCrhs(const Tensor<1,dim,double> p,
	const Tensor<1,dim,double> normal_vector,
	Vector<double> &bc_rhs,
	const unsigned int system_id) const
{

	double norm = p.norm();
	double x_cord = p[0];
	double y_cord = p[1];
	const unsigned int neqn_local = num_equations.total_nEqn[system_id];
	const unsigned int nbc_local = num_equations.nBC[system_id];

	switch(bc_type)
	{
		// the size of bc_rhs == no of negative eigen values in the system
		case characteristic:
		{
			// the boundary conditions have been implemented in the form Bu = g
			assert(bc_rhs.size() == nbc_local);
			bc_rhs = 0;
			switch(neqn_local)
			{
				case 6:
				{
					switch(mesh_info.mesh_type)
					{
						double thetaW;
						case ring:
						{
							if( norm > 0.7 ) 
		 		   			bc_rhs(0) = -chi*theta1;  // is chi \alpha and same with zeta

 					   		else 
 		   					{
 		   						bc_rhs(0) = -chi*theta0;
 		   						bc_rhs(1) = -uW*normal_vector[1];
 		   					}; 
							break;
						}

						case periodic_square:
						{
 			   				// upper edge
							if (fabs(y_cord - mesh_info.yt) < 1e-10 ) 
								thetaW = theta0;
							


 			   				// lower edge
							if ( fabs(y_cord - mesh_info.yb) < 1e-10 )
								thetaW = theta1;
							

							bc_rhs(0) = -chi * thetaW;
							bc_rhs(1) = -uW * normal_vector[1];		

							break;					

						}
						default:
						{
							Assert(1 == 0 , ExcNotImplemented());
							break;
						}
					}



 		   			break;
 		   		}

 		   		case 10:
 		   		{
 		   			Assert(1 == 0, ExcNotImplemented());
 		   			Assert(mesh_info.mesh_type == ring,ExcNotImplemented());

 		   			if( norm > 0.7 ) 
 		   				bc_rhs(0) = -chi*theta1;  // is chi \alpha and same with zeta
 		   			else 
 		   			{
 		   				bc_rhs(0) = -chi*theta0;
 		   				bc_rhs(1) = -uW*normal_vector[1];
 		   			}; 

	 		   		break;
 			   	}	

 			   	case 17:
 			   	{
 			   		// the boundary conditions obtained from the mathematica file are in the form BU + g = 0.
 			   		// But the bc_rhs in the present implementation is = -g.
 			   		Assert(mesh_info.mesh_type == periodic_square, ExcNotImplemented());
 			   		double tilde_thetaW;
 			   		// signs have been reveresed as compared to Manuel's implementation

 			   		// upper edge
 			   		// tilde_thetaW is the value of w_{wall}[1,0,0,0] as per the mathematica file
 			   		if (y_cord == mesh_info.yt)
 			   			tilde_thetaW = theta0;
 			   		

 			   		// lower edge
 			   		if (y_cord == mesh_info.yb)
 			   			tilde_thetaW = theta1;

 			   		for (unsigned int m = 0 ; m < system_data[system_id].B.matrix.outerSize() ; m++)
 			   			for (Sparse_matrix::InnerIterator n(system_data[system_id].B.matrix,m); n ; ++n)
 			   			{

 			   				if (n.col() == 3)		// only provide a boundary value for the temperature
 			   					bc_rhs(n.row()) = tilde_thetaW * n.value();
 			   			}
 			   			
 			   		break;
 			   	}	

 			   	default:
 			   	{
 			   		Assert(1 == 0,ExcInternalError());
 			   		break;
 			   	}	

		    }
		
		break;
	      }
	      // the size of bc_rhs == no of variables in the system
		case odd:
		{
			switch(neqn_local)
			{
				assert(bc_rhs.size() == num_equations.total_nEqn[system_id]);
				case 6:
				{

					if( norm > 0.7 ) 
 		   			bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta

 		   			else 
 		   			{
 		   				bc_rhs(1) = chi*theta0;
 		   				bc_rhs(4) = uW*normal_vector[1];
 		   			}; 

 		   			break;
 		   		}

 		   		case 10:
 		   		{


 		   			if( norm > 0.7 ) 
 		   				bc_rhs(1) = chi*theta1;  // is chi \alpha and same with zeta
 		   			else 
 		   			{
 		   				bc_rhs(1) = chi*theta0;
 		   				bc_rhs(4) = uW*normal_vector[1];
 		   			}; 

	 		   		break;
 			   	}

 			   	case 17:
 			   	{
 			   		Assert(1 == 0, ExcMessage("not implemented due to inconsistency of the method"));
 			   		break;
 			   	}
 			}

 		   break;
 		}
 	}

 }

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::source_term(const vector<Point<dim>> &p,
			vector<Vector<double>> &value,
			const unsigned int system_id)
{
	AssertDimension(value.size(),p.size());
	const unsigned int neqn_local = num_equations.total_nEqn[system_id];

	switch(force_type)
	{
		case type1:
		{
			for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				double norm = sqrt(p[i].square());
				value[i] = 0;
				value[i][0] = (A0 + A2*norm*norm + A1*p[i][0]/norm);	
			}
			break;
		}
		case type2:
		{
			for (unsigned int i = 0 ; i < value.size() ; i++)
			{
				double r = sqrt(p[i].square());
				value[i] = 0;
				value[i][0] = A0 + A2 * r * r + A1*(1.0-5.0/18*r*r/(tau*tau))*p[i][0] / r;
			}
			break;
		}

		// no forcing
		case type3:
		{
			for (unsigned int i = 0 ; i < value.size() ; i++)
				value[i] = 0;
			break;
		}

		// forcing for the poisson heat conduction
		case type4:
		{
			// the coefficient of the polynomial(see mathematica file for details)
			const double alpha = sqrt(2.0/3.0);

			for (unsigned int i = 0 ; i < value.size(); i++)
			{
				const double y_cord = p[i][1];
				// initialize the variable
				value[i] = 0;

				// the source term for the energy equation
				value[i][3] = -alpha * pow(y_cord,2);
			}
		}
	}

	switch(system_type)
	{
		// incase of symmetric system, the force vector gets multiplied by the symmetrizer
		case symmetric:
		{
			switch(neqn_local)
			{
				case 6:
				{

					for (unsigned int i = 0 ; i < value.size() ; i ++)
						Sparse_matrix_dot_Vector(system_data[system_id].S_half,value[i]);

					break;
				}
				case 10:
				case 17:
				{
					Assert(1 == 0,ExcNotImplemented());
					break;
				}
			}
			

			break;
		}
		case un_symmetric:
			break;
		
	}

}

template<int num_flux,int dim> 
Full_matrix Base_EquationGenerator<num_flux,dim>
::build_Aminus(const Tensor<1,dim,double> normal_vector,const unsigned int system_id)
{

	return( build_InvProjector(normal_vector,system_id) 
			* system_data[system_id].Aminus_1D_Int 
			* build_Projector(normal_vector,system_id) );		

}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::build_Aminus1D(Full_matrix &Aminus_1D_Int,
				Full_matrix &Aminus_1D_Bound,
				const unsigned int system_id)
{
	Assert(Aminus_1D_Bound.rows() == num_equations.total_nEqn[system_id] 
			&& Aminus_1D_Bound.cols() == num_equations.total_nEqn[system_id],
			ExcNotInitialized());

	Assert(Aminus_1D_Int.rows() == num_equations.total_nEqn[system_id]  
			&& Aminus_1D_Int.cols() == num_equations.total_nEqn[system_id],
			ExcNotInitialized());

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

	// collects the eigenvectors corresponding to negative eigenvalues
	Ordering_Values::Ordering order(vals);
	system_data[system_id].X_minus.resize(num_equations.total_nEqn[system_id],
										 num_equations.nBC[system_id]);

	for (unsigned int i = 0 ; i < num_equations.nBC[system_id]; i ++)
		system_data[system_id].X_minus.col(i) = vecs.col(order.index(i));


}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>::generate_matrices(equation_data &system_data,
															const unsigned int system_index)
{
		if (system_type == 0)
			system_data.is_symmetric = true;
		else
			system_data.is_symmetric = false;

		// number of equations in the system
	    system_data.nEqn = num_equations.total_nEqn[system_index];
	    // system id
	    system_data.system_id = num_equations.system_id[system_index];

	    // all the data corresponding to regularized theories will be saved in files ending with R
	    if (system_data.system_id < 0)
	    	system_data.base_filename = to_string(system_data.nEqn)+"R";
	    // else we have the normal grad's equations
	    else
	    	system_data.base_filename = to_string(system_data.nEqn);

		cout <<"System ID " << system_data.system_id
			 << " #nEqn " << system_data.nEqn << endl;

		string system_dir;
		string filename;
		string filename_out;

		system_dir = "system_matrices/";

		cout << "Reading system matrices.......\n" << endl;
		for (unsigned int i = 0 ; i < dim ; i ++)
		{
			filename = system_dir +"A" + to_string(i+1) + "_" + system_data.base_filename +   ".txt";
			build_triplet(system_data.A[i],filename);

			system_data.A[i].matrix.resize(system_data.nEqn,system_data.nEqn);
			build_matrix_from_triplet(system_data.A[i]);
			print_matrix(system_data.A[i],generate_filename_to_write(system_dir,filename));
		}

		filename = system_dir + "P_"  + system_data.base_filename +".txt";
		system_data.P.matrix.resize(system_data.nEqn,system_data.nEqn);

		// build P corresponding to the equations
		build_P(system_data.P,system_index);
		print_matrix(system_data.P,generate_filename_to_write(system_dir,filename));

		/*for any other case, BC will be generated internally using characteristic splitting*/
		if (num_equations.total_nEqn[system_index] == 6 
			|| num_equations.total_nEqn[system_index] == 10)
		{
					filename = system_dir + "BC_"  + system_data.base_filename + ".txt";
					system_data.BC.matrix.resize(system_data.nEqn,system_data.nEqn);
					build_BC(system_data.BC,system_index);
					print_matrix(system_data.BC,generate_filename_to_write(system_dir,filename));			
		}

		system_data.Aminus_1D_Int.resize(system_data.nEqn,system_data.nEqn);
		system_data.Aminus_1D_Bound.resize(system_data.nEqn,system_data.nEqn);

		build_Aminus1D(system_data.Aminus_1D_Int,
						system_data.Aminus_1D_Bound,
						system_index);


		// needs to be implemented for the B system
		// i.e the system with 10 equations (unavailability of symmetrizer and the B matrix) 
		if (num_equations.total_nEqn[system_index] != 10)
		{
			filename = system_dir + "B_" + system_data.base_filename + ".txt";
			system_data.B.matrix.resize(num_equations.nBC[system_index],system_data.nEqn);
			build_triplet(system_data.B,filename);
			build_matrix_from_triplet(system_data.B);
			print_matrix(system_data.B,generate_filename_to_write(system_dir,filename));

			system_data.B_tilde_inv.resize(num_equations.nBC[system_index],
								      num_equations.nBC[system_index]);

			system_data.B_hat.resize(system_data.nEqn,
							     system_data.nEqn);

			build_B_tilde_inv(system_data.B,
						  system_data.X_minus,
						  system_data.B_tilde_inv);

			build_B_hat(system_data.B,
			  		system_data.X_minus,
			 		system_data.B_tilde_inv,
			 		system_data.B_hat);

		}

		switch(system_type)
		{
			case symmetric:
			{
				Assert(system_data.system_id != 10,
					   ExcNotImplemented()); // symmetric system not implemented for B type

				Assert(system_data.system_id != 17,
						ExcNotImplemented());  // symmetric system not implemented for G26 equations

				filename = system_dir + "S_half_" + system_data.base_filename + ".txt";

				system_data.S_half.matrix.resize(system_data.nEqn,system_data.nEqn);
				build_triplet(system_data.S_half,filename);
				build_matrix_from_triplet(system_data.S_half);
				print_matrix(system_data.S_half,generate_filename_to_write(system_dir,filename));


				filename = system_dir + "S_half_inv_" + system_data.base_filename + ".txt";

				system_data.S_half_inv.matrix.resize(system_data.nEqn,system_data.nEqn);
				build_triplet(system_data.S_half_inv,filename);
				build_matrix_from_triplet(system_data.S_half_inv);
				print_matrix(system_data.S_half_inv,generate_filename_to_write(system_dir,filename));

				system_data.Ax.matrix = system_data.A[0].matrix;
				for (unsigned int i = 0 ; i < dim ; i ++)
					system_data.A[i].matrix = system_data.S_half.matrix * system_data.A[i].matrix * system_data.S_half_inv.matrix;

				system_data.P.matrix = system_data.S_half.matrix * system_data.P.matrix * system_data.S_half_inv.matrix;
				break;
			}

			case un_symmetric:
			{
				break;
			}

		}


		cout << "Done Reading Matrices......\n" << endl;
}
