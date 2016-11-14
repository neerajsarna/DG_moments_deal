// In the following file we genereate tensor info for different sytems
// along with symmetrizer

namespace TensorInfo
{
	using namespace dealii;

	template<int dim>
	class
	Base_TensorInfo
	{
		public:
			
			struct projector_data
			{
				// matrix which stores the local projector matrix
				Full_matrix P;
			};	

			// In all the normal grad's theories, we know the number of tensors we will be using
			Base_TensorInfo();
			void reinit(const unsigned int n_tensors);

			// For certain systems, for ex the A system we do not have a direct correlation between the various
			// Ntensors and varIdx. For such system we directly provide varIdx
			void reinit(const MatrixUI &var_idx,const unsigned int n_tensors);

			MatrixUI varIdx;
			// we have the max tensorial degree because only then the projectors for individual tensors can be 
			// constructed
			const unsigned int max_tensorial_degree = 11;
			unsigned int Ntensors;

			// total number of equations in the system
			unsigned int nEqn;
			MatrixUI free_indices;
			MatrixUI free_indices_cumilative;

			// terminology belongs to the Mathematica file
			unsigned int iFullDegree(const unsigned int tensor_degree);
			unsigned int iTraces(const unsigned int tensor_degree);
			unsigned int iTensorDegree(const unsigned int tensor_degree);

			// computes the maximum tensorial degree depending upon varIdx
			unsigned int compute_max_tensorial_degree();

			// computes the number of 2D components for a particular tensorial degree
			unsigned int components_2D(const unsigned int tensorial_degree);

			// computes the number of 3D components for a particular tensorial degree
			unsigned int components_3D(const unsigned int tensor_degree);

			// computes the total number of equations for a given Ntensors
			unsigned int compute_nEqn(MatrixUI &free_indices);

			MatrixUI generate_varIdx();
			MatrixUI generate_num_free_indices(MatrixUI &varIdx);
			MatrixUI generate_cumilative_index(MatrixUI &free_indices);
			void allocate_tensor_memory(std::vector<projector_data> &tensor_project);

			// put a block of a full_matrix at a particular location in a sparse_matrix
			// idx is the diagonal position at which P has to be placed at Sp
			void SpBlock(const unsigned int idx,const Full_matrix &P,Sparse_matrix &Sp);

			void reinit_local(const double nx,const double ny,
							  std::vector<projector_data> &tensor_project);

			// same as above but for inverse projector
			void reinit_Invlocal(const double nx,const double ny,
				   				std::vector<projector_data> &tensor_project);

			Sparse_matrix reinit_global(const double nx,const double ny);

			// inverse of the projector matrix
			Sparse_matrix reinit_Invglobal(const double nx,const double ny);

			// develops the mirror of a normal vector
			Tensor<1,dim> mirror(const Tensor<1,dim,double> normal_vector) const;

			// routines to generate the symmetrizer
			Sparse_matrix S_half;

			Sparse_matrix S_half_inv;

			// create symmetrizer for the system
			void create_Symmetrizer();
			void reinit_symmetrizer(std::vector<projector_data> &tensor_Invsymmetrizer);

			// create inverse symmetrizer for the system
			void create_InvSymmetrizer();
			void reinit_Invsymmetrizer(std::vector<projector_data> &tensor_Invsymmetrizer);
			

			// compute the number of equations with the help of Ntensors
			// // functions not presently needed
			// std::vector<unsigned int> odd_id_global;
			// Full_matrix IDXtrf(const unsigned int tensor_degree);
			// bool is_Odd(const unsigned int value);
			// bool is_Odd(const double value);

			// void odd_id_local(const Full_matrix &IDXtrf,std::vector<unsigned int> &odd_id_local);
			// std::vector<unsigned int>  generate_ID_odd();
			

	};

	template<int dim>
	Base_TensorInfo<dim>::Base_TensorInfo()
	{;}

	// initialize the class depending upon the system
	template<int dim>
	void
	Base_TensorInfo<dim>::reinit(const unsigned int n_tensors)
	{
		Ntensors = n_tensors;

		varIdx.resize(Ntensors,2);
		free_indices.resize(Ntensors,1);
		Assert(max_tensorial_degree == 11,ExcNotImplemented());
		
		// dim == 1 has not been implemented
		Assert(dim > 1,ExcNotImplemented());

		// generate the varIdx, similar to the Mathematica file
		varIdx = generate_varIdx();
		Assert(compute_max_tensorial_degree() <= max_tensorial_degree,ExcMessage("Projector and Symmetrizer data not available for this tensorial degree"));

		// now we generate the number of components corresponding to all the free indices
		free_indices = generate_num_free_indices(varIdx);

		// we also generate the cumilative free indices for the development of projector
		free_indices_cumilative = generate_cumilative_index(free_indices);
		// id_odd_global = generate_ID_odd();

		// compute the total number of equations in the system
		nEqn = compute_nEqn(free_indices);

		create_Symmetrizer();
		create_InvSymmetrizer();
	}

	// same as above but with different parameters
	template<int dim>
	void 
	Base_TensorInfo<dim>::reinit(const MatrixUI &var_idx,const unsigned int n_tensors)
	{
		varIdx = var_idx;
		Ntensors = n_tensors;

		free_indices.resize(varIdx.rows(),1);
		Assert(varIdx.rows() != 0 || varIdx.cols() != 0,ExcNotInitialized());
		Assert(free_indices.size() != 0,ExcNotInitialized());
		Assert(max_tensorial_degree == 7,ExcNotImplemented());


		Assert(compute_max_tensorial_degree() <= max_tensorial_degree,ExcMessage("Projector and Symmetrizer data not available for this tensorial degree"));

		// the situation for dim ==1 has not be implemented yet
		Assert(dim > 1,ExcNotImplemented());

		// now we generate the number of components corresponding to all the free indices
		free_indices = generate_num_free_indices(varIdx);

		// we also generate the cumilative free indices for the development of projector
		free_indices_cumilative = generate_cumilative_index(free_indices);

		// id_odd_global = generate_ID_odd();

		// compute the total number of equations in the system
		nEqn = compute_nEqn(free_indices);

		// create the symmetrizer for this particular ssytem
		create_Symmetrizer();

		// create the symmetrizer for this particular system
		create_InvSymmetrizer();
	}

	template<int dim>
	unsigned int
	Base_TensorInfo<dim>
	::iFullDegree(const unsigned int tensor_degree)
	{
		return((unsigned int)(sqrt(4*tensor_degree-3)-1));
	}

	template<int dim>
	unsigned int
	Base_TensorInfo<dim>
	::iTraces(const unsigned int tensor_degree)
	{
		return((unsigned int)(pow((1 + 1 * iFullDegree(tensor_degree)/2.0),2)-tensor_degree));
	}

	template<int dim>
	unsigned int
	Base_TensorInfo<dim>
	::iTensorDegree(const unsigned int tensor_degree)
	{
		return(iFullDegree(tensor_degree)-2*iTraces(tensor_degree));
	}

	// the following function computes the maximum tensorial degree using varIdx.
	template<int dim>
	unsigned int 
	Base_TensorInfo<dim>
	::compute_max_tensorial_degree()
	{
		return(varIdx.col(1).maxCoeff());
	}

	template<int dim>
	MatrixUI
	Base_TensorInfo<dim>
	::generate_varIdx()
	{
		MatrixUI result(Ntensors,2);

		for (unsigned int i = 0 ; i < Ntensors ; i ++)
		{
			// first store the number of traces. Due to C indices we add a plus one.
			result(i,0) = iTraces(i + 1);

			// now store the number of free indices
			result(i,1) = iTensorDegree(i + 1);
		}

		return(result);
	}

	// computes the number of free components of a tensor of a particular degree in 2D
	template<int dim>
	unsigned int 
	Base_TensorInfo<dim>
	::components_2D(const unsigned int tensor_degree)
	{
		return(tensor_degree + 1);
	}

	// computes the number of free components of a tensor of a particular degree in 3D
	template<int dim>
	unsigned int
	Base_TensorInfo<dim>
	::components_3D(const unsigned int tensor_degree)
	{
		return(2 * tensor_degree + 1);
	}


	template<int dim>
	unsigned int 
	Base_TensorInfo<dim>
	::compute_nEqn(MatrixUI &free_indices)
	{
		return(free_indices.sum());
	}
	// template<int dim>
	// Full_matrix
	// Base_TensorInfo<dim>
	// ::IDXtrf(const unsigned int tensor_degree)
	// {
	// 	// total number of components in a trace free tensor of tensor_degree
	// 	const unsigned int num_components = components_3D(tensor_degree);

	// 	// we store the indices of the tensor
	// 	Full_matrix result(num_components,3);

	// 	// the loop is from o to tensor_degree
	// 	unsigned int count = 0;

	// 	for (unsigned int i = 0 ; i <= tensor_degree  ; i ++)
	// 		for (unsigned int j = 0 ; j <= std::min(i,1) ; j++ )
	// 		{
	// 			result(count) << tensor_degree - i, i - j, j ;
	// 			count ++;
	// 		}
	// }

	// // check whether a number is even or odd
	// template<int dim>
	// bool
	// Base_TensorInfo<dim>
	// ::is_Odd(const unsigned int value)
	// {
	// 	if (value %2.0 == 0)
	// 		return(false);

	// 	if ( value %2.0 != 0)
	// 		return(true);
	// }

	// // check whether a double is even or odd
	// template<int dim>
	// bool
	// Base_TensorInfo<dim>
	// ::is_Odd(const double value)
	// {
	// 	if (value %2 == 0)
	// 		return(false);

	// 	if ( value %2 != 0)
	// 		return(true);
	// }
	// // the following function accepts the IDXtrf corresponding to a particular tensorial degree
	// // and then return the id of local odd variables
	// // Description of the input parameters:
	// // 1. IDXtrf
	// // 2. vector to be assigned
	// template<int dim>
	// void 
	// Base_TensorInfo<dim>
	// ::odd_id_local(const Full_matrix &IDXtrf,std::vector<unsigned int> &odd_id_local)
	// {
	// 	const unsigned int num_components = IDXtrf.rows();

	// 	for (unsigned int i = 0 ; i < num_components ; i ++)
	// 		if(is_Odd(IDXtrf(i,0)))
	// 			odd_id_local.push_back(i);

	// }

	// // we compute the indices of the odd variables for 3D
	// template<3>
	// std::vector<unsigned int> 
	// Base_TensorInfo<3>
	// ::generate_ID_odd()
	// {
	// 	// is dependent upon the vaule of varIdx therefore need the following assertion
	// 	Assert(varIdx.rows() != 0,ExcNotInitialized());
	// 	Assert(varIdx.cols() != 0,ExcNotInitialized());
		
	// 	// in the following variable we store the global id of the odd variables
	// 	std::vector<unsigned int> odd_id_global;
	// 	// since we always need the total number of components which have been computed in the past 
	// 	const unsigned int cumilative_component = 0;

	// 	// now we loop over all the tensors we have in the system
	// 	for (unsigned int i = 0 ; i < Ntensors ; i ++)
	// 	{
	// 		// the tensorial degree of the present tensor being considered.
	// 		const unsigned int tensorial_degree = varIdx(i,1);

	// 		// the total number of components of the present tensor
	// 		const unsigned int num_components_local = components_3D(tensor_degree);

	// 		// the local id of the odd variables i.e. id of the odd variables in the present tensor
	// 		std::vector<unsigned int> odd_id_local;

	// 		// computation of the present local odd id
	// 		odd_id_local(IDXtrf(tensorial_degree),odd_id_local);

	// 		// loop over all the components of odd_id_local
	// 		for (unsigned int j = 0 ; j < odd_id_local.size() ; j++)
	// 			odd_id_global.push_back(cumilative_component + odd_id_global[j]);

	// 		// update the cumilative component
	// 		cumilative_component += num_components_local;
	// 	}

	// 	return(odd_id_global);

	// }

	// computes the total number of free indices corresponding to a tensorial degree.
	// components = (dim - 1) * n + 1
	template<int dim>
	MatrixUI
	Base_TensorInfo<dim>
	::generate_num_free_indices(MatrixUI &varIdx)
	{
		Assert(varIdx.cols() != 0,ExcNotInitialized());
		Assert(varIdx.rows() !=0 ,ExcNotInitialized());
		MatrixUI result;
		result.resize(varIdx.rows(),1);

		AssertDimension(result.size(),varIdx.rows());

		for (unsigned int i = 0 ; i < result.size() ; i ++)
			result(i) = (dim - 1)  * varIdx(i,1) + 1;
		

		return(result);
	}

	// we convert the free indices to cumilative index for the projector
	template<int dim>
	MatrixUI
	Base_TensorInfo<dim>
	::generate_cumilative_index(MatrixUI &free_indices)
	{
		MatrixUI result;
		result.resize(free_indices.size(),1);
		unsigned int cumilative_index = 0;

		// the first location will be zero
		for (unsigned int i = 0 ; i < result.size() ; i ++)
		{
			result(i) = cumilative_index;
			cumilative_index += free_indices(i);
		}

		return(result);
	}

	template<int dim>
	void
	Base_TensorInfo<dim>
	::SpBlock(const unsigned int idx,const Full_matrix &P,Sparse_matrix &Sp)
	{
		for (int i = 0 ; i < P.rows() ; i++)
			for (int j = 0 ; j < P.cols() ; j++)
				Sp.coeffRef(idx + i,idx + j) = P(i,j);

	}

	// allocate memory for individual tensors
	template<>
	void 
	Base_TensorInfo<2>
	::allocate_tensor_memory(std::vector<projector_data> &tensor_project)
	{
		// to make sure that tensor project has been initialized
		AssertDimension(tensor_project.size(),max_tensorial_degree + 1);

		for (unsigned int i = 0 ; i < max_tensorial_degree + 1 ; i ++)
			tensor_project[i].P.resize(components_2D(i),components_2D(i));
	}

	// now we allocate values to the various tensor project depending upon the normal vector

// same as above but for the inverse of the projector



	template<>
	Tensor<1,2>
	Base_TensorInfo<2>
	::mirror(const Tensor<1,2,double> normal_vector)const
	{
		double nx = normal_vector[0], ny = normal_vector[1];
		Tensor<1,2,double> mirrored_vector;
		mirrored_vector[0] = nx;
		mirrored_vector[1] = -ny;

		return mirrored_vector;
	}

	// constains the projector data for particular symmetrizers
	#include "Tensorial_Projector.h"
	#include "Tensorial_InvProjector.h"
	#include "Tensorial_Symmetrizer.h"
	#include "Tensorial_InvSymmetrizer.h"


	template<>
	Sparse_matrix
	Base_TensorInfo<2>
	::reinit_global(const double nx,const double ny)
	{

		// the block matrices to be used in the global projector
		std::vector<projector_data> tensor_project;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		// so we first initialize the tensorial projectors
		reinit_local(nx,ny,tensor_project);

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
		{

			SpBlock(free_indices_cumilative(i),tensor_project[varIdx(i,1)].P,global_Projector);
		}

		global_Projector.makeCompressed();

		return(global_Projector);
	}

	// same as above but for the inverse of projector
	template<>
	Sparse_matrix
	Base_TensorInfo<2>
	::reinit_Invglobal(const double nx,const double ny)
	{

		// the block matrices to be used in the global projector
		std::vector<projector_data> tensor_Invproject;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		// so we first initialize the tensorial projectors
		reinit_Invlocal(nx,ny,tensor_Invproject);

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
		{

			SpBlock(free_indices_cumilative(i),tensor_Invproject[varIdx(i,1)].P,global_Projector);
		}

		global_Projector.makeCompressed();

		return(global_Projector);
	}


	template<>
	void 
	Base_TensorInfo<2>
	::create_Symmetrizer()
	{
		S_half.resize(nEqn,nEqn);

		// symmetrizer for a particular tensor
		std::vector<projector_data> tensor_symmetrizer;

		tensor_symmetrizer.resize(max_tensorial_degree + 1);

		AssertDimension(tensor_symmetrizer.size(),max_tensorial_degree + 1);
		allocate_tensor_memory(tensor_symmetrizer);

		reinit_symmetrizer(tensor_symmetrizer);

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
		{

			SpBlock(free_indices_cumilative(i),tensor_symmetrizer[varIdx(i,1)].P,S_half);
		}

		S_half.makeCompressed();
	}

	template<>
	void
	Base_TensorInfo<2>
	::create_InvSymmetrizer()
	{
		S_half_inv.resize(nEqn,nEqn);

		// symmetrizer for a particular tensor
		std::vector<projector_data> tensor_Invsymmetrizer;

		tensor_Invsymmetrizer.resize(max_tensorial_degree + 1);

		AssertDimension(tensor_Invsymmetrizer.size(),max_tensorial_degree + 1);
		allocate_tensor_memory(tensor_Invsymmetrizer);

		reinit_Invsymmetrizer(tensor_Invsymmetrizer);

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_Invsymmetrizer[varIdx(i,1)].P,S_half_inv);
		

		S_half_inv.makeCompressed();

	}


}