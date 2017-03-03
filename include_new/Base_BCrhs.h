namespace BCrhs
{
	using namespace dealii;

	// an abstract class for BCrhs
	template<int dim>
	class
	Base_BCrhs
	{
	public:
		Base_BCrhs(const constant_data &constants);

		const constant_data constants;
		virtual void BCrhs(const Tensor<1,dim,double> p,
						   const Tensor<1,dim,double> normal_vector,
						   Vector<double> &bc_rhs,
						   const unsigned int b_id) = 0;

		// prescribe the attributes of the wall. The projector matrix projects a vector into local coordinates
		void assign_wall_properties(double &thetaW,double &vn,double &vt,const FullMatrix<double> &proj_vector,
									const unsigned int b_id);

		// prescribe the attributes of the inflow
		void assign_inflow_properties(double &thetaW,double &vn,double &vt,double &rho,const FullMatrix<double> &proj_vector,
									const unsigned int b_id);

		Vector<double> local_wall_velocity(const double vx,const double vy,const FullMatrix<double> &proj_vector);
	};

	template<int dim>
	Base_BCrhs<dim>
	::Base_BCrhs(const constant_data &constants)
	:
	constants(constants)
	{;}

	template<int dim>
	Vector<double>
	Base_BCrhs<dim>
	::local_wall_velocity(const double vx,const double vy,const FullMatrix<double> &proj_vector)
	{
		Vector<double> global_velocity(2);
		Vector<double> local_velocity(2);

		global_velocity(0) = vx;
		global_velocity(1) = vy;

		local_velocity(0) = proj_vector(0,0) * global_velocity(0) + proj_vector(0,1) * global_velocity(1);
		local_velocity(1) = proj_vector(1,0) * global_velocity(0) + proj_vector(1,1) * global_velocity(1);

		return(local_velocity);
	}

	template<int dim>
	void
	Base_BCrhs<dim>
	::assign_wall_properties(double &thetaW,double &vn,double &vt,const FullMatrix<double> &proj_vector, const unsigned int b_id)
	{

		bool assigned_properties = false;

		Assert(b_id <= 4,ExcNotImplemented());

		if (b_id == 0)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx0,constants.vy0,proj_vector);

			thetaW = constants.theta0;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}


		if (b_id == 1)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx1,constants.vy1,proj_vector);

			thetaW = constants.theta1;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 2)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx2,constants.vy2,proj_vector);

			thetaW = constants.theta2;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 3)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx3,constants.vy3,proj_vector);

			thetaW = constants.theta3;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			assigned_properties = true;
		}

		if (b_id == 4)
		{
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx4,constants.vy4,proj_vector);

			thetaW = constants.theta4;
			vn = wall_velocity(0);
			vt = wall_velocity(1);		
			assigned_properties = true;
		}




		Assert(assigned_properties==true,ExcMessage("Wall properties not assigned, check the implementation"));

	}

	template<int dim>
	void
	Base_BCrhs<dim>
	::assign_inflow_properties(double &thetaW,double &vn,double &vt,double &rho,const FullMatrix<double> &proj_vector,
								 const unsigned int b_id)
	{

		bool assigned_properties = false;

		Assert(b_id==101 || b_id == 102 ,ExcNotImplemented());

		// the wall ids start from 101, 102
		if (b_id == 101)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx101,constants.vy101,proj_vector);
		
			thetaW = constants.theta101;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
			rho = constants.rho101;
			assigned_properties = true;
		}


		if (b_id == 102)
		{
			// the wall velocity vector in local coordinates
			const Vector<double> wall_velocity = local_wall_velocity(constants.vx102,constants.vy102,proj_vector);

			thetaW = constants.theta102;
			vn = wall_velocity(0);
			vt = wall_velocity(1);	
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
		BCrhs_wall(const constant_data &constants,
			const Sparse_matrix &B);

		Sparse_matrix B;

		virtual void BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int b_id);
	};

	template<int dim>
	BCrhs_wall<dim>::BCrhs_wall(const constant_data &constants,
		const Sparse_matrix &B)
	:
	Base_BCrhs<dim>(constants),
	B(B)
	{;}

	template<int dim>
	void
	BCrhs_wall<dim>::BCrhs(const Tensor<1,dim,double> p,
						  const Tensor<1,dim,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		bc_rhs = 0;

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

		AssertDimension((int)bc_rhs.size(),this->constants.nBC);
		Assert(dim > 1,ExcNotImplemented());

		// temperature of the wall
		double thetaW;

		// normal velocity of the wall
		double vn;

		// tangential velocity of the wall
		double vt;

		this->assign_wall_properties(thetaW,vn,vt,proj_vector,b_id);

		const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
		//const unsigned int ID_vx = this->constants.variable_map.find("vx")->second;
		const unsigned int ID_vy = this->constants.variable_map.find("vy")->second;
	

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


	// handling of inflow boundary conditions
	template<int dim>
	class
	BCrhs_inflow:public Base_BCrhs<dim>
	{
	public:
		BCrhs_inflow(const constant_data &constants,
			const Sparse_matrix &Binflow);

		Sparse_matrix Binflow;

		virtual void BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int b_id);
	};

	template<int dim>
	BCrhs_inflow<dim>::BCrhs_inflow(const constant_data &constants,
									const Sparse_matrix &Binflow)
	:
	Base_BCrhs<dim>(constants),
	Binflow(Binflow)
	{;}

	template<int dim>
	void
	BCrhs_inflow<dim>::BCrhs(const Tensor<1,dim,double> p,
						  const Tensor<1,dim,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		// vector normal to the wall
		const double nx = normal_vector[0];
		const double ny = normal_vector[1];

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

		AssertDimension((int)bc_rhs.size(),this->constants.nBC);
		Assert(dim > 1,ExcNotImplemented());

		//temprature of the incoming distribution function
		double thetaW;

		// normal velocity of the incoming distribution function
		double vn;

		// tangential velocity of the incoming distribution function
		double vt;

		// density of the wall
		double rho;

		this->assign_inflow_properties(thetaW,vn,vt,rho,proj_vector,b_id);

		const unsigned int ID_rho = this->constants.variable_map.find("rho")->second;
		const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
		const unsigned int ID_vx = this->constants.variable_map.find("vx")->second;
		const unsigned int ID_vy = this->constants.variable_map.find("vy")->second;
	

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



}
