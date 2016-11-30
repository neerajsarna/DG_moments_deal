namespace Test_SymmetrizerData
{

	// We check the following product, symmetrizer * Inv symmetrizer
	TEST(SymmetrizerDataIntoInvSymmetrizerData,HandlesInvertibilityOfSymmetrizer)
	{
		const unsigned int dim = 2;
		TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
		const unsigned int Ntensors = 12;


		base_tensorinfo.reinit(Ntensors);

		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info;
		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info_inv;
		

		tensor_info.resize(base_tensorinfo.max_tensorial_degree+1);
		tensor_info_inv.resize(base_tensorinfo.max_tensorial_degree+1);

		base_tensorinfo.allocate_tensor_memory(tensor_info);
		base_tensorinfo.allocate_tensor_memory(tensor_info_inv);

		base_tensorinfo.reinit_symmetrizer(tensor_info);
		base_tensorinfo.reinit_Invsymmetrizer(tensor_info_inv);


		for (unsigned int j = 0 ; j < base_tensorinfo.max_tensorial_degree + 1 ; j ++ )
		{

			Full_matrix prod = tensor_info[j].P * tensor_info_inv[j].P;
			Full_matrix Id;
			Id.resize(prod.rows(),prod.cols());
			Id.setZero();

			for (unsigned int k = 0 ; k < prod.rows() ; k++)
				Id.coeffRef(k,k) = 1;

			Compare_Float_Mat(Id,prod);
		}

		

		std::cout << "Tested tensors from " << 0 << " to: " << base_tensorinfo.max_tensorial_degree << std::endl;	
	}

}

