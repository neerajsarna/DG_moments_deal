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
::fix_B_vx(system_matrix &B,
			const unsigned int system_index)
{
	const unsigned int neqn_local = num_equations.total_nEqn[system_index];
	for (unsigned int i = 0 ; i < neqn_local ; i++)
		if (i != 1)							// if the coefficient is not equal to the normal velocity
			B.matrix.coeffRef(0,i) *= epsilon;
}		

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>
::fix_BC_vx(system_matrix &BC,
			const unsigned int system_index)
{
	const unsigned int neqn_local = num_equations.total_nEqn[system_index];

	// change the first row of the BC matrix
	for (unsigned int i = 0 ; i < neqn_local ; i++)
			BC.matrix.coeffRef(1,i) *= epsilon;
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


	// all the other moment system will have dim + 2 conservation laws
	// all are BGK 
	for (unsigned int i = 4 ; i  < neqn_local ; i++)
		matrix_info.matrix.coeffRef(i,i) = 1.0/tau;
		
}


// mirror of a 2d vector
template<int num_flux,int dim> 
Tensor<1,dim,double> Base_EquationGenerator<num_flux,dim>
::mirror(const Tensor<1,dim,double> normal_vector) const
{
		double nx = normal_vector[0], ny = 0;
		

		Tensor<1,dim,double> mirrored_vector;

		if (dim == 1)
			mirrored_vector[0] = nx;

		if (dim == 2)
		{
			mirrored_vector[0] = nx;
			mirrored_vector[1] = -ny;
		}

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
		double ny = 0;


		Assert(Projector.rows() == system_data[system_id].nEqn 
				|| Projector.cols() == system_data[system_id].nEqn,ExcNotInitialized());

		build_tensorial_projector(nx,ny);
		const unsigned int neqn_local = num_equations.total_nEqn[system_id];

		Assert(neqn_local == 17,ExcNotImplemented());

					  SpBlock( 0, tensor_project[0].P, Projector );	// rho
  					  SpBlock( 1, tensor_project[1].P, Projector );	// velocity
  					  SpBlock( 3, tensor_project[0].P, Projector ); // theta
  			          SpBlock( 4, tensor_project[2].P, Projector ); // sigma
  					  SpBlock( 7, tensor_project[1].P, Projector ); // q
  				      SpBlock( 9, tensor_project[3].P, Projector ); // mijk
  				      SpBlock( 13, tensor_project[0].P, Projector ); // delta
  				      SpBlock( 14, tensor_project[2].P, Projector ); // Rij


	return Projector;
}

// builds the inverse projector
template<int num_flux,int dim> 
Sparse_matrix Base_EquationGenerator< num_flux,dim>::
build_InvProjector(const Tensor<1,dim,double> normal_vector,const unsigned int system_id)
{
	return build_Projector(mirror(normal_vector),system_id);
			
}

template<int num_flux,int dim> 
void Base_EquationGenerator<num_flux,dim>::
build_BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int system_id) 
{

	Assert(bc_rhs.size() !=0 ,ExcNotInitialized());

	double norm = p.norm();
	double x_cord = p[0];
	double y_cord = 0;
	const unsigned int neqn_local = num_equations.total_nEqn[system_id];
	const unsigned int nbc_local = num_equations.nBC[system_id];

	double thetaW;

 	// left edge
	if (fabs(x_cord - mesh_info.xl) < 1e-10)
		thetaW = theta0;


 	// right edge
	if (fabs(x_cord - mesh_info.xr) < 1e-10)
		thetaW = theta1;

	switch(bc_type)
	{
		case characteristic:
		{

			for (unsigned int m = 0 ; m < system_data[system_id].B.matrix.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(system_data[system_id].B.matrix,m); n ; ++n)
			{

 			   	// only provide a boundary value for the temperature and no condition for velocity equation
				if (n.col() == 3 && n.row() > 0)
					bc_rhs(n.row()) = -sqrt(3.0/2.0) * thetaW * n.value();
			}
			break;
		};

		// bcrhs in the odd case is the negative of g, because during the assembly we have a negative sign
		case odd:
		{
			const unsigned int num_odd = 6;

			// location of the odd variables in the variables 
			vector<unsigned int> odd_variables = {1, 5, 7, 9, 11, 15};		// caution: C-indices

			// no temperature condition for the normal velocity
			for (unsigned int i = 1 ; i < num_odd ; i ++)			
			{
				

				bc_rhs(odd_variables[i]) = -system_data[system_id].BC.matrix.coeffRef(odd_variables[i],3) * sqrt(3.0/2) * thetaW;
			}


			break;
		}

		default:
			break;
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

	for (unsigned int i = 0 ; i < value.size(); i++)
			{
				const double x_cord = p[i][0];
				// initialize the variable
				value[i] = 0;

				// the source term for the energy equation
				value[i][3] = -alpha * pow(x_cord,2);
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

		system_data.Aminus_1D_Int.resize(system_data.nEqn,system_data.nEqn);
		system_data.Aminus_1D_Bound.resize(system_data.nEqn,system_data.nEqn);

		build_Aminus1D(system_data.Aminus_1D_Int,
						system_data.Aminus_1D_Bound,
						system_index);


		// needs to be implemented for the B system
		// i.e the system with 10 equations (unavailability of symmetrizer and the B matrix) 

			filename = system_dir + "B_" + system_data.base_filename + ".txt";
			system_data.B.matrix.resize(num_equations.nBC[system_index],system_data.nEqn);
			build_triplet(system_data.B,filename);
			build_matrix_from_triplet(system_data.B);
			print_matrix(system_data.B,generate_filename_to_write(system_dir,filename));

			// build the BC matrix
			filename = system_dir + "BC_" + system_data.base_filename + ".txt";
			system_data.BC.matrix.resize(system_data.nEqn,system_data.nEqn);
			build_triplet(system_data.BC,filename);
			build_matrix_from_triplet(system_data.BC);
			print_matrix(system_data.BC,generate_filename_to_write(system_dir,filename));


			// fix the boundary conditions for normal velocity
			if (num_equations.total_nEqn[system_index] == 17)
			{
				fix_B_vx(system_data.B,system_index);
				fix_BC_vx(system_data.BC,system_index);
			}

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

		

		


		cout << "Done Reading Matrices......\n" << endl;
}
