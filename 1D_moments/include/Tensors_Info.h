namespace Tensor_info
{
	using namespace std;
	using namespace dealii;

	template<int dim> class tensor_development
	{
		public:
			tensor_development() {};
			void map_free_indices_to_variable(tensor_data &tensor_info);
	};

	template<int dim>
	void 
	tensor_development<dim>
	::map_free_indices_to_variable(tensor_data &tensor_info)
	{
				// allocate for the maximum number of free indices possible
				tensor_info.free_index_info.resize(4);

				// scalar quantity, one variable
				tensor_info.free_index_info[0] = 1;

				// vector quantity, two variables
				// sx, sy
				tensor_info.free_index_info[1] = 2;

				/*Now we will allocate higher order tensor. All of them are symmetric and trace free*/
				/*second order tensor*/
				// rxx, rxy , ryy
				tensor_info.free_index_info[2] = 3;

				/*third order tensor*/
				// mxxx, mxxy, mxyy, myyy
				tensor_info.free_index_info[3] = 4;
			
	}
}
