	// generate a ring shaped mesh internally
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_internal_ring()
	{
		std::cout << "Generating the mesh " <<std::endl;
		fflush(stdout);
		
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
	}

	// read a ring shaped mesh generated from gmsh
	template<int dim>
	void
	Base_MeshGenerator<dim>::mesh_gmsh_ring()
	{
		std::cout << "Reading the mesh " << std::endl;
		fflush(stdout);

		gridin.attach_triangulation(triangulation);
        std::ifstream f(mesh_file_name);
        gridin.read_msh(f);
        triangulation.set_all_manifold_ids_on_boundary(0);
        triangulation.set_manifold(0,boundary);
	}

