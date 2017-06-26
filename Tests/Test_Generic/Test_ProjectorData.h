// in this file we test the projectors individually
namespace Test_ProjectorData
{

	template<int dim>
	std::vector<Tensor<1,dim>> reinit_tangential_vectors(const Tensor<1,3,double> &normal_vector)
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

	// In the following routine we check the product Proj*Proj^{-1}
	TEST(ProjectorDataIntoInvProjectorData,HandlesInvertibilityOfProjector)
	{
		const unsigned int dim = 3;
		const unsigned int Ntensors = 12;
		const unsigned int nEqn = 56;
		TensorInfo::Base_TensorInfo<dim> base_tensorinfo(nEqn,Ntensors);

		base_tensorinfo.reinit();

		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info;
		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info_inv;
		
		for (unsigned int i = 0 ; i < 5 ; i ++)
		{

			Tensor<1,dim> normal_vector;

			normal_vector[0] = i * 0.1 + pow(i,2) + 1;
			normal_vector[1] = i * 0.5 + pow(i,3) + 1;
			normal_vector[2] = i * 0.5 + pow(i,3) + 5;

			normal_vector /= normal_vector.norm();

			std::vector<Tensor<1,dim,double>> tangential_vectors = 
											  base_tensorinfo.reinit_tangential_vectors(normal_vector);

			const double nx = normal_vector[0];
			const double ny = normal_vector[1];
			const double nz = normal_vector[2];

			const double tx = tangential_vectors[0][0];
			const double ty = tangential_vectors[0][1];
			const double tz = tangential_vectors[0][2];

			const double rx = tangential_vectors[1][0];
			const double ry = tangential_vectors[1][1];
			const double rz = tangential_vectors[1][2];

			base_tensorinfo.reinit_local(nx,ny,nz,
						 				tx,ty,tz,
						 				rx,ry,rz,
										tensor_info);

			for (unsigned int degree =  0 ; degree < 10 ; degree ++ )
				tensor_info[degree].P = tensor_info[degree].P 
										* base_tensorinfo.tensor_Invsymmetrizer[degree].P;

			base_tensorinfo.reinit_local(nx,tx,rx,
										 ny,ty,ry,
										 nz,tz,rz,
										 tensor_info_inv);

			for (unsigned int degree =  0 ; degree < 10 ; degree ++ )
				tensor_info_inv[degree].P = base_tensorinfo.tensor_symmetrizer[degree].P
											* tensor_info_inv[degree].P;

			for (unsigned int j = 0 ; j < 10; j ++ )
			{
				Full_matrix prod = tensor_info[j].P * tensor_info_inv[j].P;
				Full_matrix Id;
				Id.resize(prod.rows(),prod.cols());
				Id.setZero();

				for (unsigned int k = 0 ; k < prod.rows() ; k++)
					Id.coeffRef(k,k) = 1;

				Compare_Float_Mat(Id,prod);
			}

		}	
	}

	TEST(Symmetrizer,HandlesSymmetrizer)
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


		// if the symmetrizer is right then the system matrices should be symmetric
		for (unsigned int space = 0 ; space < dim ; space ++)
			Check_Symmetricity(System[0].system_data.A[space].matrix);

		Check_Symmetricity(System[0].system_data.P.matrix);


		Check_Symmetricity(System[0].base_tensorinfo.S_half);

	}

	TEST(DirectionalFlux,HandlesDirectionalFlux)
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


		Tensor<1,dim> normal_vector;

		normal_vector[0] = 0.1;
		normal_vector[1] = 0.5;
		normal_vector[2] = 0.3;

		normal_vector /= normal_vector.norm();


		Full_matrix An = System[0].base_tensorinfo.reinit_Invglobal(normal_vector) 
						 * System[0].system_data.Ax.matrix * 
						 System[0].base_tensorinfo.reinit_global(normal_vector);


		Full_matrix An_manuel;

		An_manuel.resize(System[0].nEqn,System[0].nEqn);
		An_manuel.setZero();

		for (unsigned int space = 0 ; space < dim ; space++)
			An_manuel += System[0].system_data.A[space].matrix * normal_vector[space];


		Compare_Float_Mat(An_manuel,An);



	}
}
