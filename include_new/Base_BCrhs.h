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

		// prescribe the attributes of the wall
		void assign_wall_properties(double &thetaW,double &vn,double &vt,const unsigned int b_id);

		// prescribe the attributes of the inflow
		void assign_inflow_properties(double &thetaW,double &vn,double &vt,,double &rho,const unsigned int b_id);
	};

	template<int dim>
	Base_BCrhs<dim>
	::Base_BCrhs(const constant_data &constants)
	:
	constants(constants)
	{;}

	template<int dim>
	void
	Base_BCrhs<dim>
	::assign_wall_properties(double &thetaW,double &vn,double &vt,const unsigned int b_id)
	{

		bool assigned_properties = false;

		Assert(b_id <= 4,ExcNotImplemented());

		if (b_id == 0)
		{
			thetaW = constants.theta0;
			vn = constants.vn0;
			vt = constants.vt0;	
			assigned_properties = true;
		}


		if (b_id == 1)
		{
			thetaW = constants.theta1;
			vn = constants.vn1;
			vt = constants.vt1;	
			assigned_properties = true;
		}

		if (b_id == 2)
		{
			thetaW = constants.theta2;
			vn = constants.vn2;
			vt = constants.vt2;	
			assigned_properties = true;
		}

		if (b_id == 3)
		{
			thetaW = constants.theta3;
			vn = constants.vn3;
			vt = constants.vt3;	
			assigned_properties = true;
		}

		if (b_id == 4)
		{
			thetaW = constants.theta4;
			vn = constants.vn4;
			vt = constants.vt4;	
			assigned_properties = true;
		}


		Assert(assigned_properties==true,ExcMessage("Wall properties not assigned, check the implementation"));

	}

	template<int dim>
	void
	Base_BCrhs<dim>
	::assign_inflow_properties(double &thetaW,double &vn,double &vt,double &rho,const unsigned int b_id)
	{

		bool assigned_properties = false;

		Assert(b_id=>101 && b_id <= 102 ,ExcNotImplemented());

		// the wall ids start from 101, 102
		if (b_id == 101)
		{
			thetaW = constants.theta101;
			vn = constants.vn101;
			vt = constants.vt101;	
			rho = constants.rho101;
			assigned_properties = true;
		}


		if (b_id == 102)
		{
			thetaW = constants.theta102;
			vn = constants.vn102;
			vt = constants.vt102;	
			rho = constants.rho102;
			assigned_properties = true;
		}

		Assert(assigned_properties==true,ExcMessage("Inflow properties not assigned, check the implementation"));

	}


	// BCrhs using characteristic variables
	template<int dim>
	class
	BCrhs_char:public Base_BCrhs<dim>
	{
	public:
		BCrhs_char(const constant_data &constants,
			const Sparse_matrix &B);

		Sparse_matrix B;

		virtual void BCrhs(const Tensor<1,dim,double> p,
			const Tensor<1,dim,double> normal_vector,
			Vector<double> &bc_rhs,
			const unsigned int b_id);
	};

	template<int dim>
	BCrhs_char<dim>::BCrhs_char(const constant_data &constants,
		const Sparse_matrix &B)
	:
	Base_BCrhs<dim>(constants),
	B(B)
	{;}

	template<int dim>
	void
	BCrhs_char<dim>::BCrhs(const Tensor<1,dim,double> p,
						  const Tensor<1,dim,double> normal_vector,
						  Vector<double> &bc_rhs,
						  const unsigned int b_id)
	{
		AssertDimension((int)bc_rhs.size(),this->constants.nBC);
		Assert(dim > 1,ExcNotImplemented());

		// temprature of the wall
		double thetaW;

		// normal velocity of the wall
		double vn;

		// tangential velocity of the wall
		double vt;

		this->assign_wall_properties(thetaW,vn,vt,b_id);

		const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
		const unsigned int ID_vx = this->constants.variable_map.find("vx")->second;
		const unsigned int ID_vy = this->constants.variable_map.find("vy")->second;
	

		bc_rhs = 0;

		

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

	template<int dim>
		class
		BCrhs_odd:public Base_BCrhs<dim>
		{
		public:
			BCrhs_odd(const constant_data &constants,
							Sparse_matrix &B,	
							const MatrixUI &odd_ID);

			Sparse_matrix B;
			
			const MatrixUI odd_ID;
			virtual void BCrhs(const Tensor<1,dim,double> p,
				const Tensor<1,dim,double> normal_vector,
				Vector<double> &bc_rhs,
				const unsigned int b_id);
		};

	template<int dim>
		BCrhs_odd<dim>::BCrhs_odd(const constant_data &constants,
			Sparse_matrix &B,
			const MatrixUI &odd_ID)
		:
		Base_BCrhs<dim>(constants),
		B(B),
		odd_ID(odd_ID)
		{;}

	template<int dim>
		void
		BCrhs_odd<dim>::BCrhs(const Tensor<1,dim,double> p,
							const Tensor<1,dim,double> normal_vector,
							Vector<double> &bc_rhs,
							const unsigned int b_id)
		{
		// for this boundary implementation, we need to have a bc_rhs the size 
		// of which equals the number of equations
			AssertDimension((int)bc_rhs.size(),this->constants.nEqn);
					// temprature of the wall
			double thetaW;

		// normal velocity of the wall
			double vn;

		// tangential velocity of the wall
			double vt;

			this->assign_wall_properties(thetaW,vn,vt,b_id);

			const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
			const unsigned int ID_vx = this->constants.variable_map.find("vx")->second;
			const unsigned int ID_vy = this->constants.variable_map.find("vy")->second;
			bc_rhs = 0;

			

			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
				for (Sparse_matrix::InnerIterator n(B,m); n ; ++n)
				{
					const double odd_coeff = B.coeffRef(n.row(),odd_ID(n.row(),0));


					// the assumption in this implementation is the order in which the variables
					// appear in the odd_ID is the same order in which the boundary conditions have been specified.

					
 			   		// only provide a boundary value for the temperature and no condition for velocity equation
					if (n.col() == ID_theta && n.row() > 0)
					{
						
						bc_rhs(odd_ID(n.row(),0)) += -sqrt(3.0/2.0) * thetaW * n.value()/odd_coeff;

					}

					if (n.col() == ID_vy && n.row() > 0)
						bc_rhs(odd_ID(n.row(),0)) += vt * n.value()/odd_coeff;
				}

			bc_rhs(ID_vx) += vn;

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

		this->assign_inflow_properties(thetaW,vn,vt,rho,b_id);

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
					Assert(n.row() != 0,ExcMessage("incorrect boundary matrix. The coefficient for the tangential velocity should not 
													be present in the first equation"));
					bc_rhs(n.row()) += vt * n.value();
				}

				// prescribe the density of the fluid 
				if (n.col() == ID_rho)
					bc_rhs(n.row()) += rho * n.value();

				// we now prescribe the normal velocity
				if (n.row() == 0 && n.col() == ID_vx)
					bc_rhs(n.rows()) += vn * n.value();
				
			}

		}
}
