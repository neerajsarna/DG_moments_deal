template<int dim> void Solver_DG<dim>::h_adapt()
{
  switch(refinement)
  {
    case global:
    {
      triangulation.refine_global(1);
      break;
    }
    case adaptive:
    {
      FEValues<dim> fe_midpoint(mapping,finite_element,QMidpoint<dim>(),update_values|update_quadrature_points);
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      Vector<double> estimated_error_per_cell(triangulation.n_active_cells());
      int counter = 0;

      for(; cell != endc ; cell++)
      {

        Tensor<2,dim> Y;
        std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;

        active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
          GeometryInfo<dim>::max_children_per_face);

        fe_midpoint.reinit(cell);

        Point<dim> this_center = fe_midpoint.quadrature_point(0);

        vector<Vector<double>> this_midpoint_value(1);
        this_midpoint_value[0].reinit(generate_systemA<dim>::nEqn);

        fe_midpoint.get_function_values (solution,
          this_midpoint_value);

        Tensor<1,dim> projected_gradient;

        for (unsigned int face_no = 0 ; face_no < GeometryInfo<dim>::faces_per_cell ; ++face_no)
          if (!cell->face(face_no)->at_boundary())
          {
            const typename DoFHandler<dim>::face_iterator itr_face = cell->face(face_no);             
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);            


            if (neighbor->active())
              active_neighbors.push_back(neighbor);

            else
            {
              for (unsigned int subface_no = 0 ; subface_no < itr_face->n_children() ; subface_no ++)
                active_neighbors.push_back(cell->neighbor_child_on_subface(face_no,subface_no));
            }
          }

          typename vector<typename DoFHandler<dim>::active_cell_iterator>::const_iterator  neighbor_ptr = active_neighbors.begin();
          for(; neighbor_ptr != active_neighbors.end() ; neighbor_ptr++)
          {
            const typename DoFHandler<dim>::active_cell_iterator
            neighbor = *neighbor_ptr;

            fe_midpoint.reinit(neighbor);

            vector<Vector<double>> neighbor_midpoint_value(1);
            neighbor_midpoint_value[0].reinit(generate_systemA<dim>::nEqn);

            const Point<dim> neighbor_center = fe_midpoint.quadrature_point(0);

            fe_midpoint.get_function_values (solution,
              neighbor_midpoint_value);

            Tensor<1,dim> y = neighbor_center - this_center;
            double distance = y.norm();
            y /= distance;

            for(unsigned int i = 0 ; i < dim ; i ++)
              for (unsigned int j = 0 ; j < dim ; j ++)
                Y[i][j] += y[i] * y[j];

                projected_gradient += (neighbor_midpoint_value[0][0] - this_midpoint_value[0][0])/distance * y; //adjusting the mesh as per the temperature
              }

              assert(determinant(Y) != 0);
              const Tensor<2,dim> Y_inverse = invert(Y);
              Tensor<1,dim> gradient = Y_inverse * projected_gradient;

              estimated_error_per_cell[counter] = pow(cell->diameter(),1+dim * 1.0/2) * sqrt(gradient.norm_square());
              counter++;
            }



            GridRefinement::refine_and_coarsen_fixed_number (triangulation,
             estimated_error_per_cell,
             0.7, 0.1);
          
          cell = triangulation.begin_active();

          for (; cell != endc ; cell ++)
            for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face ++)
              if (cell->face(face)->at_boundary())
              {
                cell->clear_coarsen_flag();
                cell->set_refine_flag();
              }

            triangulation.execute_coarsening_and_refinement();        

            break;
          }

          case adaptive_kelly:
          {
          Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
          KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(ngp),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);

          GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.7, 0.1);

          // in the analytical treatment of the equations, we found certain instabilities for 
          // small curvature and thus it is necessary to always execute refinement in the interior cells
          typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

          for (; cell != endc ; cell ++)
            for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face ++)
              if (cell->face(face)->at_boundary())
              {
                cell->clear_coarsen_flag();
                cell->set_refine_flag();
              }
          

          triangulation.execute_coarsening_and_refinement();      
            break;
          }

        }
        
      }
      