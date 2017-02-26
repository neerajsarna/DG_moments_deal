// in the following function we test the value of the dof at the vertices
namespace TestVertexDof
{
	using namespace dealii;
	
	TEST(VertexDof,HandlesVertexDof)
	{
		const unsigned int dim = 2;
		Assert(dim == 2,ExcNotImplemented());

		std::string folder_name = "../system_matrices/";
	
		const unsigned int p = 1;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		Triangulation<dim> triangulation;

		Point<dim> p1;
		Point<dim> p2;
		std::vector<unsigned int > repetitions(dim);

		p1(0) = 0;
		p1(1) = 0;

		p2(0) = 2;
		p2(1) = 2;

		repetitions[0] = 1;
		repetitions[1] = 1;

            //The diagonal of the rectangle is the time joining p1 and p2
		GridGenerator::subdivided_hyper_rectangle(triangulation,
													repetitions,
													p1,
													p2);

       		
        std::cout << "#Cells " << triangulation.n_active_cells() << std::endl;
        

		
		FE_DGQ<dim> finite_element(p);
		DoFHandler<dim> dof_handler(triangulation);
		
		dof_handler.distribute_dofs(finite_element);

		std::cout << "#Dofs" << dof_handler.n_dofs() << std::endl;

		typename DoFHandler<dim>::active_cell_iterator cell =	dof_handler.begin_active(),
													  end_c = 	dof_handler.end();


		std::cout << "#Dofs per vertex " << finite_element.dofs_per_vertex <<std::endl;
		std::vector<Point<dim>> support_points = finite_element.get_generalized_support_points();

		Vector<double> trial_solution(dof_handler.n_dofs());
		Vector<double> solution_value(1);

		trial_solution = 1;



		for (unsigned int i = 0 ; i < support_points.size(); i++)
			

		 for (; cell != end_c ; cell++)
		 	for (unsigned int vertex = 0 ; vertex < GeometryInfo<dim>::vertices_per_cell ; vertex ++)
		 	{
		 		dealii::VectorTools::point_value(dof_handler,
												trial_solution,
												cell->vertex(vertex),
												solution_value);
		 		std::cout << "value at the vertex " << solution_value << std::endl;
		 	}
		// 		std::cout << cell->vertex_dof_index(vertex,1) << std::endl;

	}
}
