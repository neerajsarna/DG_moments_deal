TEST(DISABLED_BCrhsWallG20,HandlesBCrhsWallG20)
{
	const unsigned int dim = 3;

	std::string folder_name = "../system_matrices/";
	Constants::Base_Constants constants(input_file);
	std::vector<Develop_System::System<dim>> System;

	for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
			constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));


	
	for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System[i].initialize_system();




	double theta0 = 1.0;
	double theta1 = 2.0;
	double theta2 = 3.0;
	double theta3 = 0.5;
	double theta4 = 2.0;


	double vx0 = 1.0;
	double vx1 = 2.0;
	double vx2 = 0.5;
	double vx3 = 1.0;
	double vx4 = 3.0;


	double vy0 = 0.2;
	double vy1 = 0.6;
	double vy2 = 1.5;
	double vy3 = 2.8;
	double vy4 = 3.1;


	double vz0 = 0.3;
	double vz1 = 0.8;
	double vz2 = 2.0;
	double vz3 = 2.3;
	double vz4 = 1.3;


	unsigned int b_id;
	Tensor<1,dim> normal_vector;

	normal_vector[0] = 0.1;
	normal_vector[1] = 0.3;
	normal_vector[2] = 0.5;

	normal_vector /= normal_vector.norm();

	std::vector<Tensor<1,3,double>> tangential_vectors = System[0].bcrhs_wall.reinit_tangential_vectors(normal_vector);

	const double nx = normal_vector[0];
	const double ny = normal_vector[1];
	const double nz = normal_vector[2];

	const double tx = tangential_vectors[0][0];
	const double ty = tangential_vectors[0][1];
	const double tz = tangential_vectors[0][2];

	const double rx = tangential_vectors[1][0];
	const double ry = tangential_vectors[1][1];
	const double rz = tangential_vectors[1][2];

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


	Tensor<1,dim> p;

	p[0] = 0;
	p[1] = 0;
	p[2] = 0;


	// velocity of the wall in local coordinates
	Vector<double> wall_velocity = System[0].bcrhs_wall.local_wall_velocity(vx0,vy0,vz0,proj_vector);
	Vector<double> bcrhs_manuel(System[0].nBC);
	Vector<double> bcrhs(System[0].nBC);


	double vn = wall_velocity(0);
	double vt = wall_velocity(1);
	double vr = wall_velocity(2);

	// b_id = 0
	b_id = 0;


	System[0].bcrhs_wall.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = 4 *  sqrt(2/(15 * M_PI)) * theta0 * sqrt(3.0/2.0);
	bcrhs_manuel(4) = (2/15.0) * sqrt(2/M_PI) * theta0 * sqrt(3.0/2.0);
	bcrhs_manuel(5) = -(1/15.0) * sqrt(2/M_PI) * theta0 * sqrt(3.0/2.0);
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);


	wall_velocity = System[0].bcrhs_wall.local_wall_velocity(vx1,vy1,vz1,proj_vector);
	vn = wall_velocity(0);
	vt = wall_velocity(1);
	vr = wall_velocity(2);

	// b_id = 0
	b_id = 1;

	System[0].bcrhs_wall.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = 4 *  sqrt(2/(15 * M_PI)) * theta1 * sqrt(3.0/2.0);
	bcrhs_manuel(4) = (2/15.0) * sqrt(2/M_PI) * theta1* sqrt(3.0/2.0);
	bcrhs_manuel(5) = -(1/15.0) * sqrt(2/M_PI) * theta1 * sqrt(3.0/2.0);
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);


	wall_velocity = System[0].bcrhs_wall.local_wall_velocity(vx2,vy2,vz2,proj_vector);
	vn = wall_velocity(0);
	vt = wall_velocity(1);
	vr = wall_velocity(2);

	// b_id = 0
	b_id = 2;

	System[0].bcrhs_wall.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = 4 *  sqrt(2/(15 * M_PI)) * theta2 * sqrt(3.0/2.0);
	bcrhs_manuel(4) = (2/15.0) * sqrt(2/M_PI) * theta2* sqrt(3.0/2.0);
	bcrhs_manuel(5) = -(1/15.0) * sqrt(2/M_PI) * theta2 * sqrt(3.0/2.0);
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);

	wall_velocity = System[0].bcrhs_wall.local_wall_velocity(vx3,vy3,vz3,proj_vector);
	vn = wall_velocity(0);
	vt = wall_velocity(1);
	vr = wall_velocity(2);

	// b_id = 0
	b_id = 3;

	System[0].bcrhs_wall.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = 4 *  sqrt(2/(15 * M_PI)) * theta3 * sqrt(3.0/2.0);
	bcrhs_manuel(4) = (2/15.0) * sqrt(2/M_PI) * theta3* sqrt(3.0/2.0);
	bcrhs_manuel(5) = -(1/15.0) * sqrt(2/M_PI) * theta3 * sqrt(3.0/2.0);
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);

	wall_velocity = System[0].bcrhs_wall.local_wall_velocity(vx4,vy4,vz4,proj_vector);
	vn = wall_velocity(0);
	vt = wall_velocity(1);
	vr = wall_velocity(2);

	// b_id = 0
	b_id = 4;

	System[0].bcrhs_wall.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = 4 *  sqrt(2/(15 * M_PI)) * theta4 * sqrt(3.0/2.0);
	bcrhs_manuel(4) = (2/15.0) * sqrt(2/M_PI) * theta4* sqrt(3.0/2.0);
	bcrhs_manuel(5) = -(1/15.0) * sqrt(2/M_PI) * theta4 * sqrt(3.0/2.0);
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);

}


TEST(DISABLED_BCrhsInflowG20,HandlesBCrhsInflowG20)
{
	const unsigned int dim = 3;

	std::string folder_name = "../system_matrices/";
	Constants::Base_Constants constants(input_file);
	std::vector<Develop_System::System<dim>> System;

	for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
			constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));


	
	for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System[i].initialize_system();


 	double rho101 = 2.0;
 	double theta101 = 0.6;
 	double vx101 = 1;
 	double vy101 = 0.1;
 	double vz101 = 2;

 
 	double rho102 = 1.0;
 	double theta102 = 0.2;
 	double vx102 = 3.5;
 	double vy102 = 0.3;
 	double vz102 = 0.5;


	unsigned int b_id;
	Tensor<1,dim> normal_vector;

	normal_vector[0] = 0.1;
	normal_vector[1] = 0.3;
	normal_vector[2] = 0.5;

	normal_vector /= normal_vector.norm();

	std::vector<Tensor<1,3,double>> tangential_vectors = System[0].bcrhs_wall.reinit_tangential_vectors(normal_vector);

	const double nx = normal_vector[0];
	const double ny = normal_vector[1];
	const double nz = normal_vector[2];

	const double tx = tangential_vectors[0][0];
	const double ty = tangential_vectors[0][1];
	const double tz = tangential_vectors[0][2];

	const double rx = tangential_vectors[1][0];
	const double ry = tangential_vectors[1][1];
	const double rz = tangential_vectors[1][2];

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


	Tensor<1,dim> p;

	p[0] = 0;
	p[1] = 0;
	p[2] = 0;


	// velocity of the wall in local coordinates
	Vector<double> wall_velocity = System[0].bcrhs_inflow.local_wall_velocity(vx101,vy101,vz101,proj_vector);
	Vector<double> bcrhs_manuel(System[0].nBC);
	Vector<double> bcrhs(System[0].nBC);


	double vn = wall_velocity(0);
	double vt = wall_velocity(1);
	double vr = wall_velocity(2);

	// b_id = 0
	b_id = 101;


	System[0].bcrhs_inflow.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = -sqrt(2/M_PI) * rho101 - theta101 * sqrt(3.0/2.0) * (1/sqrt(3 * M_PI)) + vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = -rho101/sqrt(5 * M_PI)- theta101 * sqrt(3.0/2.0) * (-7/sqrt(30 * M_PI));
	bcrhs_manuel(4) = (2/(5 * sqrt(3 * M_PI))) * rho101 - theta101 * sqrt(3.0/2.0) * (-sqrt(2)/(5 * sqrt(M_PI)));
	bcrhs_manuel(5) = -rho101/(5 * sqrt(3 * M_PI))- theta101 * sqrt(3.0/2.0) * (1/(5 * sqrt(2 * M_PI)));
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);


	wall_velocity = System[0].bcrhs_inflow.local_wall_velocity(vx102,vy102,vz102,proj_vector);
	vn = wall_velocity(0);
	vt = wall_velocity(1);
	vr = wall_velocity(2);

	// b_id = 0
	b_id = 102;

	System[0].bcrhs_inflow.BCrhs(p,normal_vector,bcrhs,b_id);

	bcrhs_manuel(0) = -sqrt(2/M_PI) * rho102 - theta102 * sqrt(3.0/2.0) * (1/sqrt(3 * M_PI)) + vn;
	bcrhs_manuel(1) = -vt/sqrt(M_PI);
	bcrhs_manuel(2) = -vr/sqrt(M_PI);
	bcrhs_manuel(3) = -rho102/sqrt(5 * M_PI)- theta102 * sqrt(3.0/2.0) * (-7/sqrt(30 * M_PI));
	bcrhs_manuel(4) = (2/(5 * sqrt(3 * M_PI))) * rho102 - theta102 * sqrt(3.0/2.0) * (-sqrt(2)/(5 * sqrt(M_PI)));
	bcrhs_manuel(5) = -rho102/(5 * sqrt(3 * M_PI))- theta102 * sqrt(3.0/2.0) * (1/(5 * sqrt(2 * M_PI)));
	bcrhs_manuel(6) = 0;

	Compare_Float_Vec(bcrhs_manuel,bcrhs);


}

TEST(BCrhsWallKineticG20,HandlesKineticFluxG20)
{
		const unsigned int dim = 1;

	std::string folder_name = "../system_matrices/";
	Constants::Base_Constants constants(input_file);
	std::vector<Develop_System::System<dim>> System;

	for(int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System.push_back(Develop_System::System<dim>(constants.constants_num,constants.constants_sys.nEqn[i],
			constants.constants_sys.nBC[i],constants.constants_sys.Ntensors[i],folder_name));


	
	for (int i = 0 ; i < constants.constants_sys.total_systems ; i++)
		System[i].initialize_system();


	Tensor<1,dim> p;
	Tensor<1,dim> normal_vector;

	p[0] = 0.5;
	normal_vector[0] = 1.0;

	Vector<double> boundary_value(System[0].nEqn);

	System[0].bcrhs_wall_kinetic.BCrhs(p,normal_vector,boundary_value,0);

	std::cout << "Value at boundary " << boundary_value << std::endl;

}