namespace Test_TensorInfo
{
	using namespace dealii;

	#include "develop_Aminus.h"

	Full_matrix Projector_manuel(const double nx,const double ny,
								   const unsigned int nEqn)
	{
		Full_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		global_Projector << 1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.*nx,1.*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   -1.*ny,1.*nx,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0.7886751345948129*pow(nx,2) - 0.21132486540518713*pow(ny,2),
   1.4142135623730951*nx*ny,-0.21132486540518713*pow(nx,2) + 
    0.7886751345948129*pow(ny,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.*nx*ny,
   0.7071067811865476*(pow(nx,2) - pow(ny,2)),1.*nx*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   -0.21132486540518713*pow(nx,2) + 0.7886751345948129*pow(ny,2),
   -1.4142135623730951*nx*ny,0.7886751345948129*pow(nx,2) - 
    0.21132486540518713*pow(ny,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.*nx,1.*ny,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,-1.*ny,1.*nx,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0.6051359354609693*pow(nx,3) - 0.5516289482287868*nx*pow(ny,2),
   1.4476551742303836*pow(nx,2)*ny - 0.18387631607626226*pow(ny,3),
   -0.1838763160762623*pow(nx,3) + 1.4476551742303838*nx*pow(ny,2),
   -0.5516289482287868*pow(nx,2)*ny + 0.6051359354609696*pow(ny,3),0,0,0,0,0,0,0,0,
   0,0,0,0,0,-0.9728885676134938*pow(nx,2)*ny + 0.18387631607626226*pow(ny,3),
   0.4825517247434612*pow(nx,3) - 1.1489797655631846*nx*pow(ny,2),
   1.1489797655631848*pow(nx,2)*ny - 0.48255172474346125*pow(ny,3),
   -0.1838763160762623*pow(nx,3) + 0.9728885676134943*nx*pow(ny,2),0,0,0,0,0,0,0,0,
   0,0,0,0,0,-0.18387631607626226*pow(nx,3) + 0.9728885676134937*nx*pow(ny,2),
   -1.1489797655631846*pow(nx,2)*ny + 0.4825517247434612*pow(ny,3),
   0.48255172474346125*pow(nx,3) - 1.1489797655631848*nx*pow(ny,2),
   0.9728885676134942*pow(nx,2)*ny - 0.1838763160762623*pow(ny,3),0,0,0,0,0,0,0,0,0,
   0,0,0,0,0.5516289482287868*pow(nx,2)*ny - 0.6051359354609693*pow(ny,3),
   -0.18387631607626226*pow(nx,3) + 1.4476551742303836*nx*pow(ny,2),
   -1.4476551742303838*pow(nx,2)*ny + 0.1838763160762623*pow(ny,3),
   0.6051359354609696*pow(nx,3) - 0.5516289482287868*nx*pow(ny,2),0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0.7886751345948129*pow(nx,2) - 0.21132486540518713*pow(ny,2),
   1.4142135623730951*nx*ny,-0.21132486540518713*pow(nx,2) + 
    0.7886751345948129*pow(ny,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.*nx*ny,
   0.7071067811865476*(pow(nx,2) - pow(ny,2)),1.*nx*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   -0.21132486540518713*pow(nx,2) + 0.7886751345948129*pow(ny,2),
   -1.4142135623730951*nx*ny,0.7886751345948129*pow(nx,2) - 
    0.21132486540518713*pow(ny,2);

		return(global_Projector);
	}

	Full_matrix InvProjector_manuel(const double nx,const double ny,
									  const unsigned int nEqn)
	{
		Full_matrix Inv_global_Projector;
		Inv_global_Projector.resize(nEqn,nEqn);

		Inv_global_Projector <<1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.*nx,-1.*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   1.*ny,1.*nx,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   1.3660254037844386*pow(nx,2) + 0.3660254037844386*pow(ny,2),-2.*nx*ny,
   0.3660254037844386*pow(nx,2) + 1.3660254037844386*pow(ny,2),0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,1.4142135623730951*nx*ny,1.4142135623730951*(pow(nx,2) - pow(ny,2)),
   -1.4142135623730951*nx*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0.3660254037844386*pow(nx,2) + 1.3660254037844386*pow(ny,2),2.*nx*ny,
   1.3660254037844386*pow(nx,2) + 0.3660254037844386*pow(ny,2),0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,1.*nx,-1.*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.*ny,1.*nx,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,1.8689147936150907*pow(nx,3) + 0.7121499099253351*nx*pow(ny,2),
   -4.182444560994602*pow(nx,2)*ny - 0.7121499099253351*pow(ny,3),
   0.7121499099253351*pow(nx,3) + 4.182444560994602*nx*pow(ny,2),
   -0.7121499099253351*pow(nx,2)*ny - 1.8689147936150907*pow(ny,3),0,0,0,0,0,0,0,0,
   0,0,0,0,0,2.3436814002319806*pow(nx,2)*ny + 0.7121499099253349*pow(ny,3),
   2.3436814002319806*pow(nx,3) - 2.5509130706879564*nx*pow(ny,2),
   -2.5509130706879564*pow(nx,2)*ny + 2.3436814002319806*pow(ny,3),
   0.7121499099253349*pow(nx,3) + 2.3436814002319806*nx*pow(ny,2),0,0,0,0,0,0,0,0,0,
   0,0,0,0,0.7121499099253349*pow(nx,3) + 2.343681400231981*nx*pow(ny,2),
   2.5509130706879573*pow(nx,2)*ny - 2.343681400231981*pow(ny,3),
   2.343681400231981*pow(nx,3) - 2.5509130706879573*nx*pow(ny,2),
   -2.343681400231981*pow(nx,2)*ny - 0.7121499099253349*pow(ny,3),0,0,0,0,0,0,0,0,0,
   0,0,0,0,0.7121499099253349*pow(nx,2)*ny + 1.868914793615091*pow(ny,3),
   0.7121499099253349*pow(nx,3) + 4.182444560994603*nx*pow(ny,2),
   4.182444560994603*pow(nx,2)*ny + 0.7121499099253349*pow(ny,3),
   1.868914793615091*pow(nx,3) + 0.7121499099253349*nx*pow(ny,2),0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   1.3660254037844386*pow(nx,2) + 0.3660254037844386*pow(ny,2),-2.*nx*ny,
   0.3660254037844386*pow(nx,2) + 1.3660254037844386*pow(ny,2),0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,1.4142135623730951*nx*ny,1.4142135623730951*(pow(nx,2) - pow(ny,2)),
   -1.4142135623730951*nx*ny,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0.3660254037844386*pow(nx,2) + 1.3660254037844386*pow(ny,2),2.*nx*ny,
   1.3660254037844386*pow(nx,2) + 0.3660254037844386*pow(ny,2);

		return(Inv_global_Projector);
	}

	// TEST(iFullDegree,HandlesiFullDegree)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 20;
	// 	VectorXd ntensors(tensors) ;
		 

	// 	for (unsigned int i = 0 ; i < tensors ; i ++)
	// 		ntensors(i) = i + 1;

	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(3);

	// 	VectorXd result(tensors);
	// 	VectorXd result_manuel(ntensors);


	// 	result_manuel<< 0, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7;

	// 	for (unsigned int i = 0 ; i < ntensors.size() ; i ++)
	// 	{
	// 		result(i) = base_tensorinfo.iFullDegree(ntensors(i));
	// 		EXPECT_EQ(result(i),result_manuel(i));
	// 	}

	// }

	// TEST(iTraces,HandlesiTraces)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 20;
	// 	VectorXd ntensors(tensors) ;
		 

	// 	for (unsigned int i = 0 ; i < tensors ; i ++)
	// 		ntensors(i) = i + 1;

	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(3);

	// 	VectorXd result(tensors);
	// 	VectorXd result_manuel(ntensors);


	// 	result_manuel<< 0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0;

	// 	for (unsigned int i = 0 ; i < ntensors.size() ; i ++)
	// 	{
	// 		result(i) = base_tensorinfo.iTraces(ntensors(i));
	// 		EXPECT_EQ(result(i),result_manuel(i));
	// 	}		
	// }

	// TEST(iTensorDegree,HandlesiTensorDegree)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 20;
	// 	VectorXd ntensors(tensors) ;


	// 	for (unsigned int i = 0 ; i < tensors ; i ++)
	// 		ntensors(i) = i + 1;

	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(3);

	// 	VectorXd result(tensors);
	// 	VectorXd result_manuel(ntensors);


	// 	result_manuel<< 0, 1, 0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 2, 4, 6, 1, 3, 5, 7;

	// 	for (unsigned int i = 0 ; i < ntensors.size() ; i ++)
	// 	{
	// 		result(i) = base_tensorinfo.iTensorDegree(ntensors(i));
	// 		EXPECT_EQ(result(i),result_manuel(i));
	// 	}		
	// }

	// TEST(varIdx,HandlesvarIdx)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 4;
	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(tensors);

	// 	MatrixUI varIdx_result(20,2);

	// 	varIdx_result << 0, 0, 0, 1, 1, 0, 0, 2, 1, 1, 0, 3, 2, 0, 1, 2, 0, 4, 2, 1, 1, 3, 0, 
	// 						5, 3, 0, 2, 2, 1, 4, 0, 6, 3, 1, 2, 3, 1, 5, 0, 7;

	// 	for (unsigned int i = 0 ; i < varIdx_result.rows(); i ++)
	// 	{
	// 		EXPECT_EQ(base_tensorinfo.varIdx(i,0),varIdx_result(i,0)) << "Traces do not match";
	// 		EXPECT_EQ(base_tensorinfo.varIdx(i,1),varIdx_result(i,1)) << "Free indices do not match"; 
	// 	}
	// }

	// TEST(FreeIndices,HandlesFreeIndices)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 20;
	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(tensors);

	// 	MatrixUI free_indices(tensors,1);

	// 	free_indices << 1, 2, 1, 3, 2, 4, 1, 3, 5, 2, 4, 6, 1, 3, 5, 7, 2, 4, 6, 8;

	// 	for (unsigned int i = 0 ; i < free_indices.size() ; i++)
	// 		EXPECT_EQ(base_tensorinfo.free_indices(i),free_indices(i))<<"free indices do not match";
	// }

	// TEST(FreeIndicesCumilative,HandlesFreeIndicesCumilative)
	// {
	// 	const unsigned int dim = 2;
	// 	const unsigned int tensors = 20;
	// 	TensorInfo::Base_TensorInfo<dim> base_tensorinfo;
	// 	base_tensorinfo.reinit(tensors);

	// 	MatrixUI free_indices_cumilative(tensors,1);		
	// 	free_indices_cumilative << 0, 1, 3, 4, 7, 9, 13, 14, 17, 22, 24, 28, 34, 35, 38, 43, 50, 52, 56,62;

	// 	for (unsigned int i = 0 ; i < free_indices_cumilative.size(); i ++)
	// 		EXPECT_EQ(base_tensorinfo.free_indices_cumilative(i),free_indices_cumilative(i))<< "Cumilative free indices do not match";
	// }

	TEST(Projector, HandlesG26Projector)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

		std::vector<Tensor<1,dim>> normal_vectors(5);

		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			normal_vectors[i][0] = i * 0.1 + pow(i,2);
			normal_vectors[i][1] = i * 0.5 + pow(i,3);

			Sparse_matrix Projector = G26.build_Projector(normal_vectors[i]);
			Full_matrix result = Projector_manuel(normal_vectors[i][0],normal_vectors[i][1],G26.constants.nEqn);
			
			Compare_Float_Mat(Projector,result);
		}
	}

	TEST(InvProjector,HandlesG26InvProjector)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

		std::vector<Tensor<1,dim>> normal_vectors(5);

		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			normal_vectors[i][0] = i * 0.1 + pow(i,2);
			normal_vectors[i][1] = i * 0.5 + pow(i,3);

			Sparse_matrix Inv_Projector = G26.build_InvProjector(normal_vectors[i]);
			Full_matrix result = InvProjector_manuel(normal_vectors[i][0],normal_vectors[i][1],G26.constants.nEqn);
			
			Compare_Float_Mat(Inv_Projector,result);
		}

	}

	TEST(AminusAnG26,HandlesAminusAnG26)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

		std::vector<Tensor<1,dim>> normal_vectors(5);

		for (unsigned int i = 0 ; i < 5 ; i ++)
		{
			normal_vectors[i][0] = i * 0.1 + pow(i,2) +1 ;
			normal_vectors[i][1] = i * 0.5 + pow(i,3) + 1;

			normal_vectors[i] = normal_vectors[i]/normal_vectors[i].norm();

			Full_matrix Aminus = G26.build_Aminus(normal_vectors[i]);
			Full_matrix result = Aminus_manuel(normal_vectors[i][0],normal_vectors[i][1],G26.constants.nEqn);
			
			Compare_Float_Mat(Aminus,result);
		}	
	}

	}