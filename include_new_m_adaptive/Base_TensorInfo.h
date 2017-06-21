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
			Base_TensorInfo(const int nEqn,const int Ntensors);
			void reinit();

			// For certain systems, for ex the A system we do not have a direct correlation between the various
			// Ntensors and varIdx. For such system we directly provide varIdx
			void reinit_systemA(const MatrixUI &var_idx,const unsigned int n_tensors);

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
			unsigned int components(const unsigned int tensorial_degree);


			// computes the total number of equations for a given Ntensors
			unsigned int compute_nEqn(MatrixUI &free_indices);

			MatrixUI generate_varIdx();
			MatrixUI generate_varIdx_systemA();
			MatrixUI generate_num_free_indices(MatrixUI &varIdx);
			MatrixUI generate_cumilative_index(MatrixUI &free_indices);
			void allocate_tensor_memory(std::vector<projector_data> &tensor_project);

			// put a block of a full_matrix at a particular location in a sparse_matrix
			// idx is the diagonal position at which P has to be placed at Sp
			void SpBlock(const unsigned int idx,const Full_matrix &P,Sparse_matrix &Sp);

			void reinit_local(const double nx,const double ny,const double nz,
							  const double tx,const double ty,const double tz,
							  const double rx,const double ry,const double rz,
							  std::vector<projector_data> &tensor_project);

			// development of projector matrix for the 2D case
			Sparse_matrix reinit_global(const Tensor<1,dim> &normal_vector);
			
			// inverse of the projector matrix
			Sparse_matrix reinit_Invglobal(const Tensor<1,dim> &normal_vector);
			
			// computes the tangential vectors for a particular given normal vector
			std::vector<Tensor<1,dim>> reinit_tangential_vectors(const Tensor<1,dim> &normal_vector);

			// routines to generate the symmetrizer
			Sparse_matrix S_half;

			Sparse_matrix S_half_inv;

			// create symmetrizer for the system in the 2D setting
			void create_Symmetrizer();
			void reinit_symmetrizer(std::vector<projector_data> &tensor_symmetrizer);

			// create inverse symmetrizer for the system in the 2D setting
			void create_InvSymmetrizer();
			void reinit_Invsymmetrizer(std::vector<projector_data> &tensor_Invsymmetrizer);

			// symmetrizer for a particular tensor
			std::vector<projector_data> tensor_symmetrizer;
			std::vector<projector_data> tensor_Invsymmetrizer;
			
	};

	template<int dim>
	Base_TensorInfo<dim>::Base_TensorInfo(const int nEqn,const int Ntensors)
	:
	Ntensors(Ntensors),
	nEqn(nEqn)
	{

		// in the 1D case, number of tensors = number of equations
		if (dim == 1)
			AssertDimension(Ntensors,nEqn);

	}

	// initialize the class depending upon the system
	// specialization for the 2D case.
	template<int dim>
	void
	Base_TensorInfo<dim>::reinit()
	{
	
		// this variable is similar to the mathematica file and is computed in a similar way
		varIdx.resize(Ntensors,2);
		free_indices.resize(Ntensors,1);

		Assert(max_tensorial_degree == 11,ExcNotImplemented());

		// generate the varIdx, similar to the Mathematica file
		if (Ntensors == 3)
			varIdx = generate_varIdx_systemA();

		else
			varIdx = generate_varIdx();

		//Assert(compute_max_tensorial_degree() <= max_tensorial_degree,ExcMessage("Projector and Symmetrizer data not available for this tensorial degree"));

		// now we generate the number of components corresponding to all the free indices
		free_indices = generate_num_free_indices(varIdx);

		// we also generate the cumilative free indices for the development of projector
		free_indices_cumilative = generate_cumilative_index(free_indices);
		
		create_Symmetrizer();
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

	template<int dim>
	MatrixUI
	Base_TensorInfo<dim>
	::generate_varIdx_systemA()
	{
		MatrixUI result;
		result.resize(this->Ntensors,2);

		Assert(result.rows() !=0 ,ExcNotInitialized());
		Assert(result.cols() !=0 ,ExcNotInitialized());
		
		// theta
		result(0,0) = 1;
		result(0,1) = 0;

		// q
		result(1,0) = 0;
		result(1,1) = 1;

		// R
		result(2,0) = 1;
		result(2,1) = 2;

		return(result);
	}


	// computes the number of free components of a tensor of a particular degree in 2D
	template<int dim>
	unsigned int 
	Base_TensorInfo<dim>
	::components(const unsigned int tensor_degree)
	{
		return((dim-1) * tensor_degree + 1);
	}

	template<int dim>
	unsigned int 
	Base_TensorInfo<dim>
	::compute_nEqn(MatrixUI &free_indices)
	{
		Assert(dim > 1,ExcMessage("Incorrect dimension"));
		return(free_indices.sum());
	}

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

	// give a location defined by idx, Place the full matrix P at the location 
	// (idx,idx) in the sparse matrix Sp
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
	template<int dim>
	void 
	Base_TensorInfo<dim>
	::allocate_tensor_memory(std::vector<projector_data> &tensor_project)
	{
		// to make sure that tensor project has been initialized
		AssertDimension(tensor_project.size(),max_tensorial_degree + 1);

		// for every tensorial degree we allocate the memory
		for (unsigned int i = 0 ; i < max_tensorial_degree + 1 ; i ++)
			tensor_project[i].P.resize(components(i),components(i));
	}

	template<>
	std::vector<Tensor<1,1>>
	Base_TensorInfo<1>
	::reinit_tangential_vectors(const Tensor<1,1,double> &normal_vector )
	{
		// for any given situation there are two tangential vectors
		std::vector<Tensor<1,1>> tangential_vectors(2);

		// for the 1D case there is no sense of a tangential direction
		for (int i = 0 ; i < 2 ; i ++)
			tangential_vectors[i][0] = 0.0;

		return(tangential_vectors);
	}

	template<>
	std::vector<Tensor<1,2>>
	Base_TensorInfo<2>
	::reinit_tangential_vectors(const Tensor<1,2,double> &normal_vector)
	{
		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double tx = -ny;
		const double ty = nx;

		std::vector<Tensor<1,2>> tangential_vectors(2);

		tangential_vectors[0][0] = tx;
		tangential_vectors[0][1] = ty;

		// for the 2D case, there is no notion of a third normal direction
		tangential_vectors[1][0] = 0.0;
		tangential_vectors[1][1] = 0.0;

		return(tangential_vectors);
	}

	template<>
	std::vector<Tensor<1,3>>
	Base_TensorInfo<3>
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

		t /= t.norm();

		// the third orthogonal vector
		Tensor<1,3,double> r = cross_product_3d(t,normal_vector);

		r /= r.norm();

		std::vector<Tensor<1,3>> tangential_vectors(2);

		tangential_vectors[0] = t;
		tangential_vectors[1] = r;

		return(tangential_vectors);

	}


	// constains the projector data for particular symmetrizers
	#include "Tensorial_Projector.h"
	//#include "Tensorial_InvProjector.h"
	#include "Tensorial_Symmetrizer.h"
	#include "Tensorial_InvSymmetrizer.h"

	// development of the Projector matrix for the 1D case
	template<>
	Sparse_matrix
	Base_TensorInfo<1>
	::reinit_global(const Tensor<1,1> &normal_vector)
	{

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);
		const double nx = normal_vector[0];

		// in the 1D case, the normal vector can either be plus or minus one
		Assert(fabs(nx - 1.0) < 1e-10 || fabs(nx + 1.0) < 1e-10,ExcMessage("Incorrect normal vector"));


		for (unsigned int i = 0 ; i < Ntensors ; i++)
		{
			// if the moment is even then it is unaffected by the rotation of the axis
			if (i % 2 == 0)
				global_Projector.coeffRef(i,i) = pow(nx,2);

			// if the moment is odd then it simply changes sign 
			else
				global_Projector.coeffRef(i,i) = nx;
		}


		return(global_Projector);
	}



	// development of the Projector matrix for the 2D case
	template<>
	Sparse_matrix
	Base_TensorInfo<2>
	::reinit_global(const Tensor<1,2> &normal_vector)
	{

		// the block matrices to be used in the global projector
		std::vector<projector_data> tensor_project;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		std::vector<Tensor<1,2,double>> tangential_vectors = reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = 0;

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = 0;

		const double rx = 0;
		const double ry = 0;
		const double rz = 0;

		// initialize the projector for individual tensors
		// all the other components have been set to zero
		reinit_local(nx,ny,nz,
					 tx,ty,tz,
					 rx,ry,rz,
					tensor_project);

		for (unsigned int degree = 0 ; degree <= max_tensorial_degree ; degree++)
			tensor_project[degree].P = tensor_project[degree].P * tensor_Invsymmetrizer[degree].P;

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_project[varIdx(i,1)].P,global_Projector);
		

		global_Projector.makeCompressed();

		return(global_Projector);
	}

	// same as above but for the 3D case
	template<>
	Sparse_matrix
	Base_TensorInfo<3>
	::reinit_global(const Tensor<1,3> &normal_vector)
	{

		// the block matrices to be used in the global projector
		std::vector<projector_data> tensor_project;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		// the first entry of the vector contains the perpendicular vector t and the second 
		// entry contains the perpendicular entry r
		std::vector<Tensor<1,3,double>> tangential_vectors = reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = normal_vector[2];

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = tangential_vectors[0][2];

		const double rx = tangential_vectors[1][0];
		const double ry = tangential_vectors[1][1];
		const double rz = tangential_vectors[1][2];

		// initialize the projector for individual tensors
		// all the other components have been set to zero
		reinit_local(nx,ny,nz,
					 tx,ty,tz,
					 rx,ry,rz,
					tensor_project);

		// now we multiply by the symmetrizer
		for (unsigned int degree = 0 ; degree <= max_tensorial_degree ; degree++)
			tensor_project[degree].P = tensor_project[degree].P * tensor_Invsymmetrizer[degree].P;

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_project[varIdx(i,1)].P,global_Projector);
		

		global_Projector.makeCompressed();

		return(global_Projector);
	}



	template<>
	Sparse_matrix
	Base_TensorInfo<1>
	::reinit_Invglobal(const Tensor<1,1,double> &normal_vector)
	{
		// in the 1D case, the Projector and the Inv projector are the same.
		return(reinit_global(normal_vector));
	}


	// same as above but for the inverse of projector, specialization for the 2D case
	template<>
	Sparse_matrix
	Base_TensorInfo<2>
	::reinit_Invglobal(const Tensor<1,2,double> &normal_vector)
	{

		// the block matrices to be used in the global projector
		std::vector<projector_data> tensor_Invproject;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		std::vector<Tensor<1,2,double>> tangential_vectors = reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = 0;

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = 0;

		const double rx = 0;
		const double ry = 0;
		const double rz = 0;
		

		// so we first initialize the tensorial projectors
		reinit_local(nx,tx,rx,
						ny,ty,ry,
						nz,tz,rz,
						tensor_Invproject);

		// now we multiply with appropriate powers of the symmetrizer
		for (unsigned int degree = 0 ; degree <= max_tensorial_degree ; degree++)
			tensor_Invproject[degree].P = tensor_symmetrizer[degree].P * tensor_Invproject[degree].P ;

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_Invproject[varIdx(i,1)].P,global_Projector);
		

		global_Projector.makeCompressed();

		return(global_Projector);
	}

	template<>
	Sparse_matrix
	Base_TensorInfo<3>
	::reinit_Invglobal(const Tensor<1,3,double> &normal_vector)
	{

		// the block matrices to be used in the global projector
		// Caution: since the assembly could be done in a parallel way, it is necessary
		// to make the following entry local
		std::vector<projector_data> tensor_Invproject;

		//the global projector for the equations
		Sparse_matrix global_Projector;
		global_Projector.resize(nEqn,nEqn);

		std::vector<Tensor<1,3,double>> tangential_vectors = reinit_tangential_vectors(normal_vector);

		const double nx = normal_vector[0];
		const double ny = normal_vector[1];
		const double nz = normal_vector[2];

		const double tx = tangential_vectors[0][0];
		const double ty = tangential_vectors[0][1];
		const double tz = tangential_vectors[0][2];

		const double rx = tangential_vectors[1][0];
		const double ry = tangential_vectors[1][1];
		const double rz = tangential_vectors[1][2];
		

		// so we first initialize the tensorial projectors
		reinit_local(nx,tx,rx,
						ny,ty,ry,
						nz,tz,rz,
						tensor_Invproject);

		// now we multiply with appropriate powers of the symmetrizer
		for (unsigned int degree = 0 ; degree <= max_tensorial_degree ; degree++)
			tensor_Invproject[degree].P = tensor_symmetrizer[degree].P * tensor_Invproject[degree].P ;

		Assert(global_Projector.rows() != 0,ExcNotInitialized());
		Assert(global_Projector.cols() != 0,ExcNotInitialized());

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_Invproject[varIdx(i,1)].P,global_Projector);
		

		global_Projector.makeCompressed();

		return(global_Projector);
	}


	// specialization for the 1D case. Develops the inverse of the Projector


	template<>
	void 
	Base_TensorInfo<1>
	::create_Symmetrizer()
	{
		S_half.resize(nEqn,nEqn);


		// the moment system is already symmetric
		for (unsigned int i =0 ; i < nEqn ; i++)
			S_half.coeffRef(i,i) = 1.0;

		S_half.makeCompressed();
	}


	template<>
	void 
	Base_TensorInfo<2>
	::create_Symmetrizer()
	{
		S_half.resize(nEqn,nEqn);


		tensor_symmetrizer.resize(max_tensorial_degree + 1);

		AssertDimension(tensor_symmetrizer.size(),max_tensorial_degree + 1);
		allocate_tensor_memory(tensor_symmetrizer);

		reinit_symmetrizer(tensor_symmetrizer);

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_symmetrizer[varIdx(i,1)].P,S_half);

		S_half.makeCompressed();
	}

	template<>
	void 
	Base_TensorInfo<3>
	::create_Symmetrizer()
	{
		S_half.resize(nEqn,nEqn);


		tensor_symmetrizer.resize(max_tensorial_degree + 1);

		AssertDimension(tensor_symmetrizer.size(),max_tensorial_degree + 1);
		allocate_tensor_memory(tensor_symmetrizer);

		reinit_symmetrizer(tensor_symmetrizer);

		// then using the tensorial projectors we initialize the global projector
		for (unsigned int i = 0 ; i < Ntensors ; i++)
			SpBlock(free_indices_cumilative(i),tensor_symmetrizer[varIdx(i,1)].P,S_half);

		S_half.makeCompressed();
	}


	template<>
	void
	Base_TensorInfo<1>
	::create_InvSymmetrizer()
	{
		S_half_inv.resize(nEqn,nEqn);

		
		// the symmetrizer is an identity matrix
		for (unsigned int i = 0 ; i < nEqn ; i++)
			S_half_inv.coeffRef(i,i) = 1;

		S_half_inv.makeCompressed();

	}


	template<int dim>
	void
	Base_TensorInfo<dim>
	::create_InvSymmetrizer()
	{
		S_half_inv.resize(nEqn,nEqn);

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
