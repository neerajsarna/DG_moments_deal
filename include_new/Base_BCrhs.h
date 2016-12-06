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

		void assign_wall_properties(double &thetaW,double &vn,double &vt,const unsigned int b_id);
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

		Assert(b_id <= 3,ExcNotImplemented());

		if (b_id == 0)
		{
			thetaW = constants.theta0;
			vn = constants.vn0;
			vt = constants.vt0;	
		}


		if (b_id == 1)
		{
			thetaW = constants.theta1;
			vn = constants.vn1;
			vt = constants.vt1;	
		}

		if (b_id == 2)
		{
			thetaW = constants.theta2;
			vn = constants.vn2;
			vt = constants.vt2;	
		}

		if (b_id == 3)
		{
			thetaW = constants.theta3;
			vn = constants.vn3;
			vt = constants.vt3;	
		}

	}

		// BCrhs for a ring with characteristic variables
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
		}
