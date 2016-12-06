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
									 const unsigned int total_cycles);
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
												 const unsigned int total_cycles)
	{
		// if its not the last cycle then we refine
		if (present_cycle != total_cycles - 1)
			triangulation.refine_global(1);
	}

	#include "MeshGenerator_Ring.h"
	#include "MeshGenerator_PeriodicSquare.h"
	#include "MeshGenerator_Square.h"
}
