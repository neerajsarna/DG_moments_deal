// In the following routine we generate the mesh
namespace MeshGenerator
{
	using namespace dealii;

	template<int dim>
	class
	Base_MeshGenerator
	{
		public:
			Base_MeshGenerator(
							   const std::string &output_file_name,
							   const std::string &mesh_filename,
							   const Mesh_type mesh_type,
							   const Problem_type problem_type,
							   const double part_x,
							   const double part_y,
							   const int initial_refinement);
        
        	const Mesh_type mesh_type;
        	const Problem_type problem_type;
        	const double part_x;
        	const double part_y;
        	const int initial_refinement;

        	SphericalManifold<dim> boundary;
        	Triangulation<dim> triangulation;
        	GridIn<dim> gridin;

			const std::string output_file_name;
			const std::string mesh_filename;
			// print the info contained in the input triangulation
			void print_mesh_info() const;
			void print_grid_points() const;

			// print the grid to an eps file
			// The index is just used to manipulate the name of the file
			void print_grid(const unsigned int index) const;

			// routine for developing the mesh
			void develop_mesh();

			// read mesh from .msh file
			void read_gmsh();

			// generate ring shaped triangulation internally
			void mesh_internal_ring();
			void mesh_internal_square(const unsigned int parts_x,const unsigned int parts_y);
			void mesh_internal_square_circular_cavity();
			void mesh_internal_line();
			void mesh_internal_box_cylinder();
			void mesh_internal_box();
			void set_bid_cube_inflow();
			
			void set_periodic_bid()const;
			void set_square_bid()const;
			void set_square_circular_cavity_bid()const;
			void set_bid_line_inflow_outflow();
			void set_bid_line_heat_conduction();
			void set_bid_box_cylinder();

			// the following routine handles the refinement of the grid
			void refinement_handling(const unsigned int present_cycle,
									 const unsigned int total_cycles);
	};

		// specialization for the 2D case
	template<>
	void 
	Base_MeshGenerator<1>::develop_mesh()
	{
		mesh_internal_line();
	}

	
	// specialization for the 2D case
	template<>
	void 
	Base_MeshGenerator<2>::develop_mesh()
	{
		bool loaded_mesh = false;

		// if a problem involves inflow and outflow, we will simply use gmsh for meshing
		if (problem_type == inflow_outflow)
		{
			loaded_mesh = true;
			read_gmsh();
		}

		else 
		{
			if (mesh_type == ring)
			{
				fflush(stdout);
				std::cout << "Developing a ring " << std::endl;
				loaded_mesh = true;
				mesh_internal_ring();
			}

			if (mesh_type == square_domain)
			{
				loaded_mesh = true;
				mesh_internal_square(part_x,part_y);
			}

			if (mesh_type == square_circular_cavity)
			{
				loaded_mesh = true;
				mesh_internal_square_circular_cavity();
			}

		}

		AssertThrow(loaded_mesh,ExcMessage("Mesh not loaded"));

	}

	template<>
	void 
	Base_MeshGenerator<3>::develop_mesh()
	{

		switch(mesh_type)
		{
			case box_cylinder:
			{
				read_gmsh();

				switch(problem_type)
				{
					case inflow_outflow:
					{
						set_bid_box_cylinder();
						break;
					}

					case heat_conduction:
					case lid_driven_cavity:
					case periodic:
					default:
					{
						Assert(1 == 0 ,ExcMessage("Test Case not implemented for this boundary"));
						break;
					}

				}
				break;
			}

			case box:
			{
				read_gmsh();
				switch(problem_type)
				{
					case inflow_outflow:
					{
						set_bid_cube_inflow();
						break;
					}

					case heat_conduction:
					case lid_driven_cavity:
					case periodic:
					default:
					{
						Assert(1 == 0 ,ExcMessage("Test Case not implemented for this boundary"));
						break;
					}

				}
				break;
			}

			case ring:
			case square_domain:
			case square_circular_cavity:
			case NACA5012:
			case line:
			default:
			{
				Assert(1 == 0 , ExcMessage("Should not have reached here"));
				break;
			}
		}

		//triangulation.refine_global(initial_refinement);
	}

	
	template<int dim>
	Base_MeshGenerator<dim>::Base_MeshGenerator(
											    const std::string &output_file_name,
											    const std::string &mesh_filename,
							   					const Mesh_type mesh_type,
							   					const Problem_type problem_type,
							   					const double part_x,
							   					const double part_y,
							   					const int initial_refinement)
	:
	mesh_type(mesh_type),
	problem_type(problem_type),
	part_x(part_x),
	part_y(part_y),
	initial_refinement(initial_refinement),
	output_file_name(output_file_name),
	mesh_filename(mesh_filename)
	{
		// first we develop the mesh
		develop_mesh();
		
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
												 const unsigned int total_cycles)
	{
		// we do not refine if we are on the last cycle 
	   if (present_cycle != total_cycles - 1)
		switch (mesh_type)
		{
			case ring:
			{
				switch(problem_type)
				{
					// heat conduction between the outer and the inner ring
					case heat_conduction:
					{
						// since we have set the manifolds therefore we can perform the refinement
						if (present_cycle != total_cycles - 1)
							triangulation.refine_global(1);

						break;
					}

					// in case of inflow and outflow the boundary ids will be prescribed by gmsh
					case inflow_outflow:
					{
						AssertThrow(total_cycles == 1,ExcMessage("Could not perform grid refinement"));
						break;
					}

					case periodic:
					case lid_driven_cavity:
					default:
					{
						AssertThrow(1 == 0,ExcNotImplemented());
						break;
					}		

				}

				break;
			}

			case square_domain:
			{

				switch(problem_type)
				{
					case periodic:
					{
                	// Since for this case the boundary ids have to be prescribed again therefore we need
                	// to generate the triangulation again.
						mesh_internal_square(part_x,part_y + 100 * (present_cycle + 1));
						break;
					}

					case heat_conduction:
					case lid_driven_cavity:
					case inflow_outflow:
					{
                	// The boundary indicators stay intact even after the refinement
						triangulation.refine_global(1);
						break;
					}

					default:
					{
						AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
						break;
					}

				}
				break;
			}

			case square_circular_cavity:
			{
				switch(problem_type)
			    {
				    case heat_conduction:
                    {
                    	// since we have provided the manifolds therefore grid refinement can be conducted.
                        triangulation.refine_global(1);
                        break;
                    }

                    case inflow_outflow:
                    {
                    	AssertThrow(total_cycles == 1,ExcMessage("Could not perform grid refinement"));
                        break;
                    }

                    case lid_driven_cavity:
                    case periodic:
                    default:
                    {
                        AssertThrow(1 == 0 ,ExcNotImplemented());
                        break;
                    } 
                }
					
				break;
			}

			// we do not generate such geometries internalies
			case box_cylinder:
			case NACA5012:
			case box:
			{
				AssertThrow(total_cycles == 1,ExcMessage("Perform the Grid Refinement in Gmsh"));
				break;
			}

			case line:
			{
				triangulation.refine_global(1);
				break;
			}
		}

	}


	template<int dim>
	void 
	Base_MeshGenerator<dim>::read_gmsh()
	{
			std::cout << "Reading Gmsh file from: " << mesh_filename << std::endl;
		    triangulation.clear();
  			gridin.attach_triangulation(triangulation);
  			std::ifstream f(mesh_filename);
  			gridin.read_msh(f);
	}

	template<>
	void
	Base_MeshGenerator<3>::print_grid_points()const
	{
		typename Triangulation<3>::active_cell_iterator cell = triangulation.begin_active(), 
														endc = triangulation.end();

		FILE *fp;
		fp = fopen("grid_points","w+");

		// print all the points corresponding to a particular plane
		for (; cell != endc ; cell++)
			for (unsigned int vertex = 0 ; vertex  < GeometryInfo<3>::vertices_per_cell ; vertex++)
			{
				Point<3> coords_vertex = cell->vertex(vertex);
				// if the cells lie on the z plane then we print the coordinates
				const double z_coord = coords_vertex(2);
				if (fabs(z_coord) < 1e-16)
					fprintf(fp, "%f %f \n",coords_vertex(0),coords_vertex(1));

			}

		fclose(fp);
			
		
	}

	// mesh generator for a ring shaped geometry
	#include "MeshGenerator_Ring.h"

	// mesh generator for a square
	#include "MeshGenerator_Square.h"

	// mesh generator for a square with a ring
	#include "MeshGenerator_Square_CircularCavity.h"

	// mesh generator periodic square
	#include "MeshGenerator_PeriodicSquare.h"

	// mesh generator line
	#include "MeshGenerator_Line.h"


	// develops a box over a sphere
	#include "MeshGenerator_BoxCylinder.h"

	// develops a cube
	#include "MeshGenerator_Box.h"

}
