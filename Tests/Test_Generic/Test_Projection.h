using namespace dealii;

TEST(Projection,HandlesProjection)
{
	    const unsigned int dim = 1;
		const unsigned int p = 1;
        ASSERT_EQ(dim,1) << "3D not implemented" << std::endl;

		Triangulation<dim> triangulation;
	        
	     // we create a line
		GridGenerator::hyper_cube(triangulation,-0.5,0.5);


		hp::FECollection<dim> finite_element;
		hp::FECollection<dim> finite_element_max_moments;

		hp::DoFHandler<dim> dof_handler(triangulation);
		hp::DoFHandler<dim> dof_handler_max_moments(triangulation);

		finite_element.push_back(FESystem<dim>(FE_DGQ<dim>(1),1,FE_Nothing<dim>(),1));
		finite_element.push_back(FESystem<dim>(FE_DGQ<dim>(1),1,FE_DGQ<dim>(1),1));
		finite_element_max_moments.push_back(FESystem<dim>(FE_DGQ<dim>(1),2));



		typename hp::DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(),
													     endc = dof_handler.end();


		typename hp::DoFHandler<1>::active_cell_iterator cell_2 = dof_handler_max_moments.begin_active(),
													     endc_2 = dof_handler_max_moments.end();



		for(; cell != endc ; cell++)
			cell->set_active_fe_index(0);

		for(; cell_2 != endc_2 ; cell_2++)
			cell_2->set_active_fe_index(0);

		dof_handler.distribute_dofs(finite_element);
		dof_handler_max_moments.distribute_dofs(finite_element_max_moments);

		QGauss<dim> quad(2);
		hp::QCollection<dim> hp_quad;

		hp_quad.push_back(quad);
		hp_quad.push_back(quad);

		Vector<double> solution(dof_handler.n_dofs());
		Vector<double> solution_max_moments(dof_handler_max_moments.n_dofs());

		ConstraintMatrix constraints;
		constraints.close();


		// now we project a function onto the first space

		VectorTools::project(dof_handler,
							 constraints,
							 hp_quad,
							 ConstantFunction<dim>(1,2),
							 solution);



		// now we interpolate between the two dof handlers
		FETools::interpolate	(dof_handler,
								 solution,
								dof_handler_max_moments,
								solution_max_moments);


		std::cout << "Dofs1 " << std::endl;
		for (unsigned int i = 0 ; i < solution.size(); i++)
			std::cout << solution(i) << std::endl;


		std::cout << "Dofs2 " << std::endl;
		for (unsigned int i = 0 ; i < solution_max_moments.size(); i++)
			std::cout << solution_max_moments(i) << std::endl;



}