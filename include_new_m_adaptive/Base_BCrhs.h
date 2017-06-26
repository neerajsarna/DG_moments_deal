namespace BCrhs
{
	using namespace dealii;

	// an abstract class for BCrhs
	template<int dim>
	class
	Base_BCrhs
	{
	public:
		Base_BCrhs(const constant_numerics &constants);

		const constant_numerics constants;
		virtual void BCrhs(const Tensor<1,dim,double> p,
						   const Tensor<1,dim,double> normal_vector,
						   Vector<double> &bc_rhs,
						   const unsigned int b_id) = 0;

		// prescribe the attributes of the wall. The projector matrix projects a vector into local coordinates
		// the following function is a specialization for the 2D case
		void assign_wall_properties(double &thetaW,double &vn,double &vt,double &vr,const FullMatrix<double> &proj_vector,
									const unsigned int b_id);

		// prescribe the attributes of the inflow
		// the following function is a specialization for the 2D case
		void assign_inflow_properties(double &thetaW,double &vn,double &vt,double &vr, double &rho,const FullMatrix<double> &proj_vector,
									const unsigned int b_id);

		// the following function computes the velocity of the flow in the normal coordinates of the wall
		// specialization for the 2D case
		Vector<double> local_wall_velocity(const double vx,const double vy,const double vz,const FullMatrix<double> &proj_vector);

		std::vector<Tensor<1,dim>> reinit_tangential_vectors(const Tensor<1,dim> &normal_vector);



	};

	template<int dim>
	Base_BCrhs<dim>
	::Base_BCrhs(const constant_numerics &constants)
	:
	constants(constants)
	{;}

	// reinitiliaze the tangential vectors for the 3D case 
	template<>
	std::vector<Tensor<1,3>>
	Base_BCrhs<3>
	::reinit_tangential_vectors(const Tensor<1,3,double> &normal_vector)
	{
		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = normal_vector[2];

		const double tx = -ny;
		const double ty = nx-nz;
		const double tz = ny;

		Tensor<1,3,double> t;

		t[0] = tx;
		t[1] = ty;
		t[2] = tz;

		// we now need to normalize the vector
		t /= t.norm();

		// the third orthogonal vector
		Tensor<1,3,double> r = cross_product_3d(t,normal_vector);

		// we now need to normalize the vector
		r /= r.norm();

		std::vector<Tensor<1,3>> tangential_vectors(2);

		tangential_vectors[0] = t;
		tangential_vectors[1] = r;

		return(tangential_vectors);

	}


	// same as above but for the 1D case
	template<>
	Vector<double>
	Base_BCrhs<1>
	::local_wall_velocity(const double vx,const double vy,const double vz,const FullMatrix<double> &proj_vector)
	{
		// since we are only considering the 2D case
		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),1);

		Vector<double> global_velocity(1);
		Vector<double> local_velocity(1);

		global_velocity(0) = vx;

		local_velocity(0) = proj_vector(0,0) * global_velocity(0);
		
		return(local_velocity);
	}


	//Computes the velocity of the fluid in the local coordinates of the wall
	// Input arguments
	// 1. velocity in the global coordinates
	// 2. the projector matrix which transforms the velocity into the local coordinates
	// specialization for the 2D case
	template<>
	Vector<double>
	Base_BCrhs<2>
	::local_wall_velocity(const double vx,const double vy,const double vz,const FullMatrix<double> &proj_vector)
	{
		// since we are only considering the 2D case
		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),2);

		Vector<double> global_velocity(2);
		Vector<double> local_velocity(2);

		global_velocity(0) = vx;
		global_velocity(1) = vy;

		local_velocity(0) = proj_vector(0,0) * global_velocity(0) + proj_vector(0,1) * global_velocity(1);
		local_velocity(1) = proj_vector(1,0) * global_velocity(0) + proj_vector(1,1) * global_velocity(1);

		return(local_velocity);
	}

	// same as above but for the 3D case
	template<>
	Vector<double>
	Base_BCrhs<3>
	::local_wall_velocity(const double vx,const double vy,const double vz,const FullMatrix<double> &proj_vector)
	{
		// since we are only considering the 2D case
		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),3);

		Vector<double> global_velocity(3);
		Vector<double> local_velocity(3);

		global_velocity(0) = vx;
		global_velocity(1) = vy;
		global_velocity(2) = vz;

		local_velocity(0) = proj_vector(0,0) * global_velocity(0)
							 + proj_vector(0,1) * global_velocity(1)
							 + proj_vector(0,2) * global_velocity(2);
		
		local_velocity(1) = proj_vector(1,0) * global_velocity(0)
							 + proj_vector(1,1) * global_velocity(1)
							 + proj_vector(1,2) * global_velocity(2);


		local_velocity(2) = proj_vector(2,0) * global_velocity(0)
							 + proj_vector(2,1) * global_velocity(1)
							 + proj_vector(2,2) * global_velocity(2);

		return(local_velocity);
	}



	// assign wall properties for a 1D case
	template<>
	void
	Base_BCrhs<1>
	::assign_wall_properties(double &thetaW,double &vn,double &vt,double &vr,
							const FullMatrix<double> &proj_vector, const unsigned int b_id)
	{

		bool assigned_properties = false;
		double vy_default = 0;
		double vz_default = 0.;

		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),1);

		// we can only have at the maximum two walls
		Assert(b_id <= 1,ExcNotImplemented());

		if (b_id == 0)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx0,vy_default,vz_default,proj_vector);

			thetaW = constants.theta0;
			vn = wall_velocity(0);
			assigned_properties = true;
		}


		if (b_id == 1)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx1,vy_default,vz_default,proj_vector);

			thetaW = constants.theta1;
			vn = wall_velocity(0);
			assigned_properties = true;
		}


		Assert(assigned_properties==true,ExcMessage("Wall properties not assigned, check the implementation"));

	}



	// for the 2D case
	template<>
	void
	Base_BCrhs<2>
	::assign_wall_properties(double &thetaW,double &vn,double &vt,double &vr,
							const FullMatrix<double> &proj_vector, const unsigned int b_id)
	{

		bool assigned_properties = false;
		double vz_default = 0.0;

		AssertDimension(proj_vector.n(),proj_vector.m());
		AssertDimension(proj_vector.m(),2);

		Assert(b_id <= 4,ExcNotImplemented());

		if (b_id == 0)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx0,constants.vy0,vz_default,proj_vector);

			thetaW = constants.theta0;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}


		if (b_id == 1)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx1,constants.vy1,vz_default,proj_vector);

			thetaW = constants.theta1;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 2)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx2,constants.vy2,vz_default, proj_vector);

			thetaW = constants.theta2;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 3)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx3,constants.vy3,vz_default,proj_vector);

			thetaW = constants.theta3;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 4)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx4,constants.vy4,vz_default,proj_vector);

			thetaW = constants.theta4;
			vn = wall_velocity(0);
			vt = wall_velocity(1);		
			assigned_properties = true;
		}




		Assert(assigned_properties==true,ExcMessage("Wall properties not assigned, check the implementation"));

	}


	// for the 2D case
	template<>
	void
	Base_BCrhs<3>
	::assign_wall_properties(double &thetaW,double &vn,double &vt,double &vr,
							const FullMatrix<double> &proj_vector, const unsigned int b_id)
	{

		bool assigned_properties = false;

		AssertDimension(proj_vector.n(),proj_vector.m());
		AssertDimension(proj_vector.m(),3);

		Assert(b_id <= 4,ExcNotImplemented());

		if (b_id == 0)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx0,constants.vy0,
																	constants.vz0,proj_vector);

			thetaW = constants.theta0;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			assigned_properties = true;
		}


		if (b_id == 1)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx1,constants.vy1,
																	constants.vz1,proj_vector);

			thetaW = constants.theta1;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			assigned_properties = true;
		}

		if (b_id == 2)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx2,constants.vy2,
																	constants.vz2,proj_vector);

			thetaW = constants.theta2;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			assigned_properties = true;
		}

		if (b_id == 3)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx3,constants.vy3,
																	constants.vz3,proj_vector);

			thetaW = constants.theta3;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			assigned_properties = true;
		}

		if (b_id == 4)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx4,constants.vy4,
																	constants.vz4,proj_vector);

			thetaW = constants.theta4;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			assigned_properties = true;
		}




		Assert(assigned_properties==true,ExcMessage("Wall properties not assigned, check the implementation"));

	}

	// specialization for the 1D case
	template<>
	void
	Base_BCrhs<1>
	::assign_inflow_properties(double &thetaW,double &vn,double &vt,double &vr,
							   double &rho,const FullMatrix<double> &proj_vector,
							   const unsigned int b_id)
	{

		bool assigned_properties = false;
		double vy_default = 0.0;
		double vz_default = 0.0;

		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),1);

		Assert(b_id==101 || b_id == 102 ,ExcNotImplemented());

		// the wall ids start from 101, 102
		if (b_id == 101)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx101,vy_default,vz_default,proj_vector);
		
			thetaW = constants.theta101;
			vn = wall_velocity(0);
			rho = constants.rho101;
			assigned_properties = true;
		}


		if (b_id == 102)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx102,vy_default,vz_default,proj_vector);

			thetaW = constants.theta102;
			vn = wall_velocity(0);
			rho = constants.rho102;
			assigned_properties = true;
		}

		Assert(assigned_properties==true,ExcMessage("Inflow properties not assigned, check the implementation"));

	}




	// specialization for the 2D case
	template<>
	void
	Base_BCrhs<2>
	::assign_inflow_properties(double &thetaW,double &vn,double &vt,double &vr,
								double &rho,const FullMatrix<double> &proj_vector,
								 const unsigned int b_id)
	{

		bool assigned_properties = false;
		const double vz_default = 0;

		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),2);

		Assert(b_id==101 || b_id == 102 ,ExcNotImplemented());

		// the wall ids start from 101, 102
		if (b_id == 101)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx101,constants.vy101,vz_default, proj_vector);
		
			thetaW = constants.theta101;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			rho = constants.rho101;
			assigned_properties = true;
		}


		if (b_id == 102)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx102,constants.vy102,vz_default, proj_vector);

			thetaW = constants.theta102;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			rho = constants.rho102;
			assigned_properties = true;
		}

		Assert(assigned_properties==true,ExcMessage("Inflow properties not assigned, check the implementation"));

	}


	// specialization for the 2D case
	template<>
	void
	Base_BCrhs<3>
	::assign_inflow_properties(double &thetaW,double &vn,double &vt,double &vr,
								double &rho,const FullMatrix<double> &proj_vector,
								 const unsigned int b_id)
	{

		bool assigned_properties = false;
		const double vz_default = 0;

		AssertDimension(proj_vector.m(),proj_vector.n());
		AssertDimension(proj_vector.m(),3);

		Assert(b_id==101 || b_id == 102 ,ExcNotImplemented());

		// the wall ids start from 101, 102
		if (b_id == 101)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx101,constants.vy101,
																	constants.vz101, proj_vector);
		
			thetaW = constants.theta101;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			rho = constants.rho101;
			assigned_properties = true;
		}


		if (b_id == 102)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx102,constants.vy102,
																	constants.vz102, proj_vector);
		
			thetaW = constants.theta102;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			vr = wall_velocity(2);
			rho = constants.rho102;
			assigned_properties = true;
		}

		Assert(assigned_properties==true,ExcMessage("Inflow properties not assigned, check the implementation"));

	}



	//Right hand side for a wall
	template<int dim>
	class
	BCrhs_wall:public Base_BCrhs<dim>
	{
	public:
		BCrhs_wall(const constant_numerics &constants,
				   const int nBC,
				   const Sparse_matrix &B);

		Sparse_matrix B;
		const int nBC;

		virtual void BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int b_id);
	};

	template<int dim>
	BCrhs_wall<dim>::BCrhs_wall(const constant_numerics &constants,
								const int nBC,
								const Sparse_matrix &B)
	:
	Base_BCrhs<dim>(constants),
	B(B),
	nBC(nBC)
	{;}


	// specialization for the 1D case
		// specialization for the 2D case
	template<>
	void
	BCrhs_wall<1>::BCrhs(const Tensor<1,1,double> p,
						  const Tensor<1,1,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		bc_rhs = 0;
		double vt_default = 0;
		double vr_default = 0;

		// vector normal to the wall
		const double nx = normal_vector[0];
		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));

		// three types of walls in the 1D case. Two fully accommodated wall and one specular wall
		Assert(b_id == 0 || b_id == 1 || b_id == 50,ExcMessage("Incorrect boundary id"));
				// we only compute the rhs in case of a non-specular wall
		if (b_id != 50)
		{


		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		
		FullMatrix<double> proj_vector;
		proj_vector.reinit(1,1);
		proj_vector(0,0) = nx;

		AssertDimension((int)bc_rhs.size(),nBC);

		// temperature of the wall
		double thetaW;

		// normal velocity of the wall
		double vn;


		this->assign_wall_properties(thetaW,vn,vt_default,vr_default,proj_vector,b_id);

		// we know the ID_theta in case of a 1D problem
		const unsigned int ID_theta = 2;
		
	
		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
 			   	// first prescribe the wall temperature
					if (n.col() == ID_theta && n.row() > 0)
						bc_rhs(n.row()) += sqrt(1.0/2.0) * thetaW * n.value();

				}

		// now prescribe the normal velocity
				bc_rhs(0) += vn;			
		}

	}


	// specialization for the 2D case
	template<>
	void
	BCrhs_wall<2>::BCrhs(const Tensor<1,2,double> p,
						  const Tensor<1,2,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		bc_rhs = 0;
		double vr_default = 0;

		// vector normal to the wall
		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));
		Assert(b_id == 0 || b_id == 1 || b_id == 2 || b_id ==3 || b_id == 50,ExcMessage("Incorrect boundary id"));
				// we only compute the rhs in case of a non-specular wall
		if (b_id != 50)
		{


		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		
		FullMatrix<double> proj_vector;
		proj_vector.reinit(2,2);
		proj_vector(0,0) = nx;
		proj_vector(0,1) = ny;
		proj_vector(1,0) = -ny;
		proj_vector(1,1) = nx;

		AssertDimension((int)bc_rhs.size(),nBC);

		// temperature of the wall
		double thetaW;

		// normal velocity of the wall
		double vn;

		// tangential velocity of the wall
		double vt;

		this->assign_wall_properties(thetaW,vn,vt,vr_default,proj_vector,b_id);

		const unsigned int ID_theta = this->constants.variable_map[1].find("theta")->second;
		const unsigned int ID_vy = this->constants.variable_map[1].find("vy")->second;
	
		AssertDimension(ID_theta,3);
		AssertDimension(ID_vy,2);

		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
 			   	// first prescribe the wall velocity
					if (n.col() == ID_theta && n.row() > 0)
						bc_rhs(n.row()) += -sqrt(3.0/2.0) * thetaW * n.value();

				// now prescribe the tangential velocity of the wall.
				// One to One relation with velocity
					if(n.col() == ID_vy && n.row() > 0)
						bc_rhs(n.row()) += vt * n.value();
				}

		// now prescribe the normal velocity
				bc_rhs(0) += vn;			
		}

	}



	// specialization for the 2D case
	template<>
	void
	BCrhs_wall<3>::BCrhs(const Tensor<1,3,double> p,
						  const Tensor<1,3,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		bc_rhs = 0;

		std::vector<Tensor<1,3,double>> tangential_vectors = this->reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = normal_vector[2];

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = tangential_vectors[0][2];

		const double rx = tangential_vectors[1][0];
		const double ry = tangential_vectors[1][1];
		const double rz = tangential_vectors[1][2];

		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));
		Assert(b_id == 0 || b_id == 1 || b_id == 2 || b_id ==3 || b_id == 4 || b_id == 50,ExcMessage("Incorrect boundary id"));
				// we only compute the rhs in case of a non-specular wall
		if (b_id != 50)
		{


		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		

		FullMatrix<double> proj_vector;
		proj_vector.reinit(3,3);
		proj_vector(0,0) = nx;
		proj_vector(0,1) = ny;
		proj_vector(0,2) = nz;

		proj_vector(1,0) = tx;
		proj_vector(1,1) = ty;
		proj_vector(1,2) = tz;

		proj_vector(2,0) = rx;
		proj_vector(2,1) = ry;
		proj_vector(2,2) = rz;

		AssertDimension((int)bc_rhs.size(),nBC);

		// temperature of the wall
		double thetaW;

		// normal velocity of the wall
		double vn;

		// tangential velocity of the wall
		double vt;

		// tangential velocity of the wall
		double vr;

		this->assign_wall_properties(thetaW,vn,vt,vr,proj_vector,b_id);

		const unsigned int ID_theta = this->constants.variable_map[2].find("theta")->second;
		const unsigned int ID_vy = this->constants.variable_map[2].find("vy")->second;
		const unsigned int ID_vz = this->constants.variable_map[2].find("vz")->second;
	

		AssertDimension(ID_vy,2);
		AssertDimension(ID_vz,3);
		AssertDimension(ID_theta,4);


		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
 			   	// first prescribe the wall velocity
					if (n.col() == ID_theta && n.row() > 0)
						bc_rhs(n.row()) += -sqrt(3.0/2.0) * thetaW * n.value();

				// now prescribe the tangential velocity of the wall.
				// One to One relation with velocity
					if(n.col() == ID_vy && n.row() > 0)
						bc_rhs(n.row()) += vt * n.value();

					if(n.col() == ID_vz && n.row() > 0)
						bc_rhs(n.row()) += vr * n.value();
				}

		// now prescribe the normal velocity
				bc_rhs(0) += vn;			
		}

	}



	// handling of inflow boundary conditions
	template<int dim>
	class
	BCrhs_inflow:public Base_BCrhs<dim>
	{
	public:
		BCrhs_inflow(const constant_numerics &constants,
					const int nBC,
					 const Sparse_matrix &Binflow);

		Sparse_matrix Binflow;
		const int nBC;

		virtual void BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int b_id);
	};

	template<int dim>
	BCrhs_inflow<dim>::BCrhs_inflow(const constant_numerics &constants,
									const int nBC,
									const Sparse_matrix &Binflow)
	:
	Base_BCrhs<dim>(constants),
	Binflow(Binflow),
	nBC(nBC)
	{;}


	// specialization for the 1D case
	template<>
	void
	BCrhs_inflow<1>::BCrhs(const Tensor<1,1,double> p,
						  const Tensor<1,1,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		// vector normal to the wall
		const double nx = normal_vector[0];
		double vt_default = 0;
		double vr_default = 0.;

		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));
		Assert(b_id ==101 || b_id == 102,ExcMessage("Incorrect boundary id"));

		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		FullMatrix<double> proj_vector;
		proj_vector.reinit(1,1);
		proj_vector(0,0) = nx;

		AssertDimension((int)bc_rhs.size(),nBC);
		

		//temprature of the incoming distribution function
		double thetaW;

		// normal velocity of the incoming distribution function
		double vn;

		// density of the wall
		double rho;

		this->assign_inflow_properties(thetaW,vn,vt_default,vr_default,rho,proj_vector,b_id);

		const unsigned int ID_rho = 0;
		const unsigned int ID_vx = 1;
		const unsigned int ID_theta = 2;

		bc_rhs = 0;

		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
		for (unsigned int m = 0 ; m < Binflow.outerSize() ; m++)
			for (Sparse_matrix::InnerIterator n(Binflow,m); n ; ++n)
			{
 			   	// first prescribe the temperature. In the inflow case, we do not have the epsilon in the first 
 			   	// equation and therefore take into account the value of the coefficients in the first equation also.
				if (n.col() == ID_theta)
					bc_rhs(n.row()) += sqrt(1.0/2.0) * thetaW * n.value();

				// prescribe the density of the fluid 
				if (n.col() == ID_rho)
					bc_rhs(n.row()) += rho * n.value();

				// we now prescribe the normal velocity
				if (n.row() == 0 && n.col() == ID_vx)
				{
					Assert(n.value() == 1,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vn * n.value();
				}
				
			}

	}


	// specialization for the 2D case
	template<>
	void
	BCrhs_inflow<2>::BCrhs(const Tensor<1,2,double> p,
						  const Tensor<1,2,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		// vector normal to the wall
		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		double vr_default = 0;

		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));
		Assert(b_id ==101 || b_id == 102,ExcMessage("Incorrect boundary id"));

		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		FullMatrix<double> proj_vector;
		proj_vector.reinit(2,2);
		proj_vector(0,0) = nx;
		proj_vector(0,1) = ny;
		proj_vector(1,0) = -ny;
		proj_vector(1,1) = nx;

		AssertDimension((int)bc_rhs.size(),nBC);

		//temprature of the incoming distribution function
		double thetaW;

		// normal velocity of the incoming distribution function
		double vn;

		// tangential velocity of the incoming distribution function
		double vt;

		// density of the wall
		double rho;

		this->assign_inflow_properties(thetaW,vn,vt,vr_default,rho,proj_vector,b_id);

		const unsigned int ID_rho = this->constants.variable_map[1].find("rho")->second;
		const unsigned int ID_vx = this->constants.variable_map[1].find("vx")->second;
		const unsigned int ID_vy = this->constants.variable_map[1].find("vy")->second;
		const unsigned int ID_theta = this->constants.variable_map[1].find("theta")->second;
	
		AssertDimension(ID_rho,0);
		AssertDimension(ID_vx,1);
		AssertDimension(ID_vy,2);
		AssertDimension(ID_theta,3);


		bc_rhs = 0;

		

		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
		for (unsigned int m = 0 ; m < Binflow.outerSize() ; m++)
			for (Sparse_matrix::InnerIterator n(Binflow,m); n ; ++n)
			{
 			   	// first prescribe the temperature. In the inflow case, we do not have the epsilon in the first 
 			   	// equation and therefore do not to take into account the value of the coefficients in the first equation also.
				if (n.col() == ID_theta)
					bc_rhs(n.row()) += -sqrt(3.0/2.0) * thetaW * n.value();

				// now prescribe the tangential velocity
				if(n.col() == ID_vy)
				{
					Assert(n.row() != 0,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vt * n.value();
				}

				// prescribe the density of the fluid 
				if (n.col() == ID_rho)
					bc_rhs(n.row()) += rho * n.value();

				// we now prescribe the normal velocity
				if (n.row() == 0 && n.col() == ID_vx)
				{
					Assert(n.value() == 1,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vn * n.value();
				}
				
			}

		}



	// specialization for the 2D case
	template<>
	void
	BCrhs_inflow<3>::BCrhs(const Tensor<1,3,double> p,
						  const Tensor<1,3,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{


		std::vector<Tensor<1,3,double>> tangential_vectors = this->reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = normal_vector[2];

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = tangential_vectors[0][2];

		const double rx = tangential_vectors[1][0];
		const double ry = tangential_vectors[1][1];
		const double rz = tangential_vectors[1][2];

		Assert(p.norm() >=0 ,ExcMessage("Incorrect point"));
		Assert(b_id ==101 || b_id == 102,ExcMessage("Incorrect boundary id"));

		// matrix which projects the wall velocties in 
		// the cartesian coordinates to local coordinate
		FullMatrix<double> proj_vector;
		proj_vector.reinit(3,3);
		proj_vector(0,0) = nx;
		proj_vector(0,1) = ny;
		proj_vector(0,2) = nz;

		proj_vector(1,0) = tx;
		proj_vector(1,1) = ty;
		proj_vector(1,2) = tz;

		proj_vector(2,0) = rx;
		proj_vector(2,1) = ry;
		proj_vector(2,2) = rz;

		AssertDimension((int)bc_rhs.size(),nBC);

		//temprature of the incoming distribution function
		double thetaW;

		// normal velocity of the incoming distribution function
		double vn;

		// tangential velocity of the incoming distribution function
		double vt;

		// tangential velocity
		double vr;

		// density of the wall
		double rho;

		this->assign_inflow_properties(thetaW,vn,vt,vr,rho,proj_vector,b_id);

		const unsigned int ID_rho = this->constants.variable_map[2].find("rho")->second;
		const unsigned int ID_vx = this->constants.variable_map[2].find("vx")->second;
		const unsigned int ID_vy = this->constants.variable_map[2].find("vy")->second;
		const unsigned int ID_vz = this->constants.variable_map[2].find("vz")->second;
		const unsigned int ID_theta = this->constants.variable_map[2].find("theta")->second;
	
		AssertDimension(ID_rho,0);
		AssertDimension(ID_vx,1);
		AssertDimension(ID_vy,2);
		AssertDimension(ID_vz,3);
		AssertDimension(ID_theta,4);


		bc_rhs = 0;

		

		// the original boundary conditions are of the form B.U = g. In the present function, we are prescribing
		// the vector g. The vector g can be defined with the help of the coefficients of B.
		for (unsigned int m = 0 ; m < Binflow.outerSize() ; m++)
			for (Sparse_matrix::InnerIterator n(Binflow,m); n ; ++n)
			{
 			   	// first prescribe the temperature. In the inflow case, we do not have the epsilon in the first 
 			   	// equation and therefore do not to take into account the value of the coefficients in the first equation also.
				if (n.col() == ID_theta)
					bc_rhs(n.row()) += -sqrt(3.0/2.0) * thetaW * n.value();

				// now prescribe the tangential velocity
				if(n.col() == ID_vy)
				{
					Assert(n.row() != 0,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vt * n.value();
				}


				if(n.col() == ID_vz)
				{
					Assert(n.row() != 0,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vr * n.value();
				}


				// prescribe the density of the fluid 
				if (n.col() == ID_rho)
					bc_rhs(n.row()) += rho * n.value();

				// we now prescribe the normal velocity
				if (n.row() == 0 && n.col() == ID_vx)
				{
					Assert(n.value() == 1,ExcMessage("incorrect boundary matrix"));
					bc_rhs(n.row()) += vn * n.value();
				}
				
			}

		}


}
