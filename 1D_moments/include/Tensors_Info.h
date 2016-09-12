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
				// allocate for the maximum number of different tensorial degree
				tensor_info.free_index_info.resize(5);

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

				/*fourthe order tensor*/
				//w[0, 4, 0, 0], w[0, 3, 1, 0], w[0, 2, 2, 0], w[0, 1, 3, 0], w[0, 0, 4,0]
				tensor_info.free_index_info[4] = 5;
			
	}
}
