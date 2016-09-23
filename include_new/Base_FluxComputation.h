// The following namespace is devoted to the computation of fluxes
namespace Base_FluxComputation
{
	using namespace dealii;

	template<int dim>
	class
	Base_FluxComputation
	{
		protected:
			Base_FluxComputation();

			// build the flux matrix corresponding to the x-direction
			Full_matrix build_Aminus1D(Sparse_matrix &Ax);

	}
}