// In the following routine we generate the mesh
namespace MeshGenerator
{
	using namespace dealii;

	template<int dim>
	class
	Base_MeshGenerator
	{
		public:
			Base_MeshGenerator(const std::string &mesh_file_name,
							   const std::string &output_file_name,
							   const constant_data &constants);
        
        	SphericalManifold<dim> boundary;
        	Triangulation<dim> triangulation;
        	GridIn<dim> gridin;

			const std::string mesh_file_name;

			const constant_data constants;

			const std::string output_file_name;
			// print the info contained in the input triangulation
			void print_mesh_info() const;

			// print the grid to an eps file
			// The index is just used to manipulate the name of the file
			void print_grid(const unsigned int index) const;

			// generate ring shaped triangulation internally
			void mesh_internal_ring();
			void mesh_internal_square(const unsigned int part_x,const unsigned int part_y);
			void mesh_internal_periodic_square(const unsigned int part_x,const unsigned int part_y);
			void set_periodic_bid()const;
			void set_square_bid()const;

			// read a ring shaped triangulation generated from gmsh
			void mesh_gmsh_ring();
			void mesh_gmsh_periodic_square();

			// the following routine handles the refinement of the grid
			void refinement_handling(const unsigned int present_cycle,
									 const unsigned int total_cycles,
									 const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution);

			Vector<double> compute_L2_manuel(const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution);

			void adapt_apriori(const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution);
	};

	template<int dim>
	Base_MeshGenerator<dim>::Base_MeshGenerator(const std::string &mesh_file_name,
											    const std::string &output_file_name,
											    const constant_data &constants)
	:
	mesh_file_name(mesh_file_name),
	constants(constants),
	output_file_name(output_file_name)
	{
		switch(constants.mesh_options)
		{
			case read_msh:
			{

				switch(constants.mesh_type)
				{
					case ring:
					{
						mesh_gmsh_ring();
						break;
					}
					case periodic_square:
					{
						mesh_gmsh_periodic_square();
						break;
					}

					default:
					{
						Assert(1 == 0,ExcMessage("should not have reached here"));
						break;
					}
				}
				
				break;
				
			}
			case generate_internal:
			{

				switch(constants.mesh_type)
				{
					case ring:
					{
						mesh_internal_ring();
						break;
					}
					case periodic_square:
					{
						mesh_internal_periodic_square(constants.part_x,constants.part_y);
						break;
					}

					case Mesh_type::square:
					{
						mesh_internal_square(constants.part_x,constants.part_y);
						break;
					}
					default:
					{
						Assert(1 == 0,ExcMessage("Should not have reached here"));
						break;
					}
				}
				break;
			}

			default:
			{
				Assert(1 == 0,ExcMessage("Should not have reached here "));
				break;
			}
		}

	}

	// the following routine prints the mesh info on the screen
	template<int dim>
	void
	Base_MeshGenerator<dim>::print_mesh_info()const
	{
    std::cout << "****************Mesh Info*************" << std::endl;

    std::string mesh_info;
    mesh_info = " #dim " + std::to_string(dim) +  ", #Cells " + std::to_string(triangulation.n_active_cells());


    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    for (; cell!=endc; ++cell)
    {
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      {
        if (cell->face(face)->at_boundary())
          boundary_count[cell->face(face)->boundary_id()]++;
      }
    }
    mesh_info += ", boundary indicators: ";
    
    for (std::map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
     it!=boundary_count.end();
     ++it)
     mesh_info += std::to_string(it->first) + "("+std::to_string(it->second) +"times)";

   if (triangulation.has_hanging_nodes())
      mesh_info += ", hanging nodes: True";
   else
      mesh_info += ", hanging nodes: False";

   mesh_info += ", vertices: " + std::to_string(triangulation.n_vertices());

   std::cout << mesh_info << std::endl;
   std::cout << "************************************" << std::endl ;
	}

	// the following routine prints the mesh to an eps file
	template<int dim>
	void
	Base_MeshGenerator<dim>::print_grid(const unsigned int index)const
	{
	  std::string file_for_grid = output_file_name + std::to_string(index);
	  std::ofstream out (file_for_grid.c_str());
      GridOut grid_out;
      grid_out.write_eps (triangulation, out);
	}

	template<int dim>
	void
	Base_MeshGenerator<dim>::refinement_handling(const unsigned int present_cycle,
												 const unsigned int total_cycles,
												 const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution)
	{
		// if its not the last cycle then we refine
		switch(constants.refinement)
		{
			case global:
			{
				if (present_cycle != total_cycles - 1)
					triangulation.refine_global(1);				

				break;
			}

			case apriori:
			{
						if (present_cycle != total_cycles - 1)
								adapt_apriori(active_cells, mapping, base_exactsolution, dof_handler, solution);

						break;
			}

			default:
			{
				Assert(1 == 0, ExcMessage("Should not have reached here"));
				break;
			}
		}

	}

	template<int dim>
	void 
	Base_MeshGenerator<dim>::adapt_apriori(const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution)
	{

        unsigned int component = constants.variable_map.find(constants.error_variable)->second;

        const unsigned int ngp = constants.p + 1;
        // error per cell of the domain
        Vector<double> error_per_cell(active_cells);      

        ComponentSelectFunction<dim> weight(component,constants.nEqn);                              // used to compute only the error in theta

        // computation of L2 error
        VectorTools::integrate_difference (*mapping,*dof_handler,solution,
          									*base_exactsolution,
          									error_per_cell,
          									QGauss<dim>(ngp),
          									VectorTools::Linfty_norm,
          									&weight);  

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
             error_per_cell,
             0.3, 0.0);

        triangulation.execute_coarsening_and_refinement();       


	}

	// same as the integrate_difference function but manuel
	template<int dim>
	Vector<double> 
	Base_MeshGenerator<dim>::compute_L2_manuel(const unsigned int active_cells,
										   const MappingQ<dim> *mapping,
										    ExactSolution::Base_ExactSolution<dim> *base_exactsolution,
									  		const DoFHandler<dim> *dof_handler,
									  		const Vector<double> &solution)
	{
		  const UpdateFlags update_flags  = update_values | update_JxW_values | update_quadrature_points;
		  QGauss<dim> quadrature(constants.p+1);

  		  FEValues<dim>  fe_v(*mapping,dof_handler->get_fe(),quadrature, update_flags);
  		  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
  											      	end_c = dof_handler->end();

  		  const int total_ngp = quadrature.size();
  		  std::vector<double> Jacobians_interior(total_ngp);
  		  Vector<double> error_per_cell(active_cells);

  		  std::vector<Vector<double>> solution_value(total_ngp);
  		  Vector<double> exact_solution_value(constants.nEqn);

  		  for (unsigned int i = 0 ; i < total_ngp ; i ++)
  		  	solution_value[i].reinit(constants.nEqn);

  		  error_per_cell = 0;
  		  unsigned int counter = 0;


  		  //std::cout << "Manuel computation **************" << std::endl;
  		  for(; cell != end_c ; cell++)
  		  {
  		  	 fe_v.reinit(cell);
  			 Jacobians_interior = fe_v.get_JxW_values();

  			 fe_v.get_function_values(solution,solution_value);

  			 for (unsigned int q = 0 ; q < total_ngp ; q++)
  			 {
  			 	base_exactsolution->vector_value(fe_v.quadrature_point(q),exact_solution_value); 
  			 	error_per_cell(counter) += pow(exact_solution_value(0)-solution_value[q](0),2) * Jacobians_interior[q];
  			 }

  			 error_per_cell(counter) = sqrt(error_per_cell(counter));
  			
  			 counter ++;
  		  }

  		  return(error_per_cell);
	}



	#include "MeshGenerator_Ring.h"
	#include "MeshGenerator_PeriodicSquare.h"
	#include "MeshGenerator_Square.h"
}
