// in this file we test the projectors individually
namespace Test_ProjectorData
{
	#include "Projector_Data.h"
	#include "InvProjector_Data.h"


	TEST(ProjectorData,HandlesProjectorData)
	{
		const unsigned int dim = 2;
		TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
		const unsigned int Ntensors = 12;


		base_tensorinfo.reinit(Ntensors);

		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info;
		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info_manuel;
		
		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			double nx = i * 0.1 + pow(i,2) + 1;
			double ny = i * 0.5 + pow(i,3) + 1;
			const double norm = sqrt(nx * nx + ny * ny);

			nx /= norm;
			ny /=norm;



			base_tensorinfo.reinit_local(nx,ny,tensor_info);
			build_Proj_manuel(tensor_info_manuel,
                        	nx,ny,base_tensorinfo.max_tensorial_degree,
                      		base_tensorinfo);

			for (unsigned int j = 0 ; j < base_tensorinfo.max_tensorial_degree + 1 ; j++)
			{

				Compare_Float_Mat(tensor_info[j].P,tensor_info_manuel[j].P);
			}


		}

		std::cout << "Tested tensors from " << 0 << " to: " << base_tensorinfo.max_tensorial_degree << std::endl;

	}

	TEST(InvProjectorData,HandlesInvProjectorData)
	{
		const unsigned int dim = 2;
		TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
		const unsigned int Ntensors = 12;


		base_tensorinfo.reinit(Ntensors);

		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info;
		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info_manuel;
		
		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			double nx = i * 0.1 + pow(i,2) + 1;
			double ny = i * 0.5 + pow(i,3) + 1;
			const double norm = sqrt(nx * nx + ny * ny);

			nx /= norm;
			ny /=norm;

			base_tensorinfo.reinit_Invlocal(nx,ny,tensor_info);
			build_InvProj_manuel(tensor_info_manuel,
                        	nx,ny,base_tensorinfo.max_tensorial_degree,
                      		base_tensorinfo);

			for (unsigned int j = 0 ; j < base_tensorinfo.max_tensorial_degree + 1 ; j++)
			{

				Compare_Float_Mat(tensor_info[j].P,tensor_info_manuel[j].P);
			}

		}

		std::cout << "Tested tensors from " << 0 << " to: " << base_tensorinfo.max_tensorial_degree << std::endl;

	}

	// We check the product of all the projectors. So we check the Product = Mat * Inv(Mat)
	TEST(ProjectorDataIntoInvProjectorData,HandlesInvertibilityOfProjector)
	{
		const unsigned int dim = 2;
		TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
		const unsigned int Ntensors = 12;


		base_tensorinfo.reinit(Ntensors);

		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info;
		std::vector<typename TensorInfo::Base_TensorInfo<dim>::projector_data> tensor_info_inv;
		
		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			double nx = i * 0.1 + pow(i,2) + 1;
			double ny = i * 0.5 + pow(i,3) + 1;
			const double norm = sqrt(nx * nx + ny * ny);

			nx /= norm;
			ny /=norm;

			base_tensorinfo.reinit_local(nx,ny,tensor_info);
			base_tensorinfo.reinit_Invlocal(nx,ny,tensor_info_inv);


			for (unsigned int j = 0 ; j < 12; j ++ )
			{

				Full_matrix prod = tensor_info[j].P * tensor_info_inv[j].P;
				Full_matrix Id;
				Id.resize(prod.rows(),prod.cols());
				Id.setZero();

				// std::cout << "<<<<<<<<<<<<<<Printing Projector>>>>>>>>>>>>>" << std::endl;
				// std::cout << tensor_info[j].P << std::endl;

				// std::cout << "<<<<<<<<<<<<<<Printing InvProjector>>>>>>>>>>>>>" << std::endl;
				// std::cout << tensor_info_inv[j].P << std::endl;
				for (unsigned int k = 0 ; k < prod.rows() ; k++)
					Id.coeffRef(k,k) = 1;

				Compare_Float_Mat(Id,prod);
			}

		}	

		std::cout << "Tested tensors from " << 0 << " to: " << base_tensorinfo.max_tensorial_degree << std::endl;	
	}
}