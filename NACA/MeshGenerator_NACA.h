	// generate a ring shaped mesh internally
	template<int dim>
	void
	MeshGenerator::Base_MeshGenerator<dim>::mesh_gmsh_NACA_channel()
	{
            triangulation.clear();
  			gridin.attach_triangulation(triangulation);
  			std::ifstream f("../../NACA/NACA5012.msh");
  			gridin.read_msh(f);
	}



