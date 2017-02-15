// in the following file we test our residual computation
namespace TestPostProc
{
	using namespace dealii;

	template<int dim> class dummy_function:public Function<dim>
	{
		public:
			dummy_function(const unsigned int nEqn):Function<dim>(nEqn),
													nEqn(nEqn){};

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;
			const unsigned int nEqn;
	};

	template<int dim> void dummy_function<dim>::vector_value(const Point<dim> &p,Vector<double> &value) const
	{
		
  		const double x = p[0];
  		const double y = p[1];
  		

  		Assert(dim == 2,ExcNotImplemented());

  		for(unsigned int eq = 0 ; eq < nEqn ; eq++)
  			value(eq) = sin(M_PI *x ) * sin(M_PI*y); 
  		
	}

	TEST(SolvingSystemA,HandlesSolvingSystemA)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		SystemA::SystemA<dim> systemA(constants.constants,folder_name);

		ExactSolution::ExactSolution_SystemA_ring<dim>  exact_solution_systemA(constants.constants,systemA.base_tensorinfo.S_half);

		ConvergenceTable convergence_table;

		Triangulation<dim> triangulation;
		

		SphericalManifold<dim> boundary;

		Point<dim> center;
		const double inner_radius = 0.5;
		const double outer_radius = 2.0;
		const unsigned int parts = 10;

		center(0) = 0.0;
		center(1) = 0.0;

		GridGenerator::hyper_shell (triangulation,
        							center, inner_radius, outer_radius,
              						parts);

		triangulation.set_all_manifold_ids_on_boundary(0);
        triangulation.set_manifold(0,boundary);

			
            // Point<dim> p1;
            // Point<dim> p2;
            // std::vector<unsigned int > repetitions(dim);

            // p1(0) = 0.5;
            // p1(1) = 0.5;

            // p2(0) = 1.0;
            // p2(1) = 1.0;

            // repetitions[0] = 10;
            // repetitions[1] = 10;

            // //The diagonal of the rectangle is the time joining p1 and p2
            // GridGenerator::subdivided_hyper_rectangle(triangulation,
            //                           	              repetitions,
            //                             	            p1,
            //                                 	        p2);


			FESystem<dim> finite_element(FE_Q<dim>(1),constants.constants.nEqn);
			DoFHandler<dim> dof_handler(triangulation);
			MappingQ<dim,dim> mapping(1);

			// in the following routine we test the l2 norm of the residual
			for (unsigned int i = 0 ; i < 5 ; i ++)
			{
			dof_handler.distribute_dofs(finite_element);
			Vector<double> solution(dof_handler.n_dofs());
			ConstraintMatrix constraints;


			std::cout << "number of active cells " << triangulation.n_active_cells() << std::endl;
			std::cout << "number of Dofs " << dof_handler.n_dofs() << std::endl;

			constraints.close();

			// project the solution to our finite element space
			VectorTools::project(dof_handler,
								constraints,
								QGauss<dim>(3),
								exact_solution_systemA,
								solution);

				

			PostProc::Base_PostProc<dim> postproc(constants.constants,&exact_solution_systemA,&dof_handler,&mapping);				
			double residual_strong_form = postproc.compute_residual(solution,finite_element,&systemA,triangulation.n_active_cells());
			double error_per_itr;
			ComponentSelectFunction<dim> weight(0,constants.constants.nEqn); 
			Vector<double> error_per_cell(triangulation.n_active_cells()) ;

			VectorTools::integrate_difference (mapping,dof_handler,solution,
          									exact_solution_systemA,
          									error_per_cell,
          									QMidpoint<dim>(),
          									VectorTools::H1_seminorm,
          									&weight);  

			// FILE *fp;
			// std::string filename = "error_variation" + std::to_string(dof_handler.n_dofs());
			// fp = fopen(filename.c_str(),"w+");


			// typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), end_c = dof_handler.end();
			// unsigned int counter = 0;

			// for (; cell != end_c ; cell++)
			// {
			// 	fprintf(fp, "%f %f %f\n",cell->center()(0),cell->center()(1),error_per_cell(counter));
			// 	counter ++;
			// }

			// fclose(fp);



			std::cout << "error in H1 norm: " << error_per_cell.l2_norm() << std::endl;
			std::cout << "residual: " << residual_strong_form << std::endl;
			triangulation.refine_global(1);	

			}


		// in the following routine we compute the residual at a particular point
		/*for (unsigned int i = 0 ; i < 7 ; i ++)
		{
			Point<dim> p;
			p(0) = 0.6;
			p(1) = 0.6;

			std::cout << "#Cells " << triangulation.n_active_cells() << std::endl;

			dof_handler.distribute_dofs(finite_element);
			Vector<double> solution(dof_handler.n_dofs());
			ConstraintMatrix constraints;
			constraints.close();

			VectorTools::project(dof_handler,
								constraints,
								QGauss<dim>(3),
								exact_solution_systemA,
								solution);

			// the results returned from the point_gradient function
			std::vector<Tensor<1,dim,double>> value_gradient(constants.constants.nEqn);
			Vector<double> value(constants.constants.nEqn);
			std::vector<Vector<double>> value_per_quad(dim,Vector<double>(constants.constants.nEqn));	
			Vector<double> difference_per_quad(constants.constants.nEqn);
			Vector<double> exact_solution_value(constants.constants.nEqn);
			difference_per_quad = 0;

			std::vector<Vector<double>> source_term_value(1,Vector<double>(constants.constants.nEqn));
			MatrixOpt::Base_MatrixOpt matrix_opt;

			// the following value has been taken from mathematica
			std::vector<double> exact_gradient_value = {0.143286, 0.197513, -0.0895235, 0.1506, -0.193412, -0.289674};
			exact_solution_systemA.vector_value(p,exact_solution_value);
			



			// the gradient value has to be multiplied by the system matrices
			VectorTools::point_gradient	(dof_handler,
										solution,
												p,
												value_gradient);

			// the point value has to be multiplied by the Production term
			VectorTools::point_value(dof_handler,
											solution,
											p,
											value);

			std::cout << "error in value " << fabs(exact_solution_value(0)-value(0)) << std::endl;


			std::vector<Point<dim>> p_vector(1);
			p_vector[0](0) = 0.6;
			p_vector[0](1) = 0.6;

			systemA.source_term(p_vector,source_term_value);

					//base_exactsolution->vector_value(q_points[q],exact_solution);

					// now we take the transpose of values so that we can use the sparse vector product already developed
					
		
					for (unsigned int eq = 0 ; eq < constants.constants.nEqn ; eq++)
						for (unsigned int space = 0 ; space < dim ; space ++)
							value_per_quad[space](eq) = value_gradient[eq][space];

					//the residual from the convective term
					for (unsigned int space = 0 ; space < dim ; space++)
						difference_per_quad += matrix_opt.Sparse_matrix_dot_Vector(systemA.system_data.A[space].matrix,
																					value_per_quad[space]);

					// the residual from the right hand side of the equation
					difference_per_quad += matrix_opt.Sparse_matrix_dot_Vector(systemA.system_data.P.matrix,value);
										   
					// computation of the l2 norm of the solution in this cell
					for (unsigned int eq = 0 ; eq < constants.constants.nEqn ; eq++)
						// contribution from the source term
						difference_per_quad(eq) -= source_term_value[0](eq);


					std::cout << "Error in all the equations " << difference_per_quad << std::endl;


			triangulation.refine_global(1);



	}*/

 }

}

