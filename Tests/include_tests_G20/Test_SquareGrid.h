namespace SquareGrid
{
	using namespace dealii;

	TEST(SquareGrid,HandlesSquareGrid)
	{
		const unsigned int dim = 2;
		ASSERT_EQ(dim,2) << "3D not implemented" << std::endl;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		std::string mesh_file = "mesh_file";
		std::string output_file = "grid";

		MeshGenerator::Base_MeshGenerator<dim> mesh_generator(mesh_file,output_file,constants.constants);

		mesh_generator.mesh_internal_square(50,50);


		typename Triangulation<dim>::active_cell_iterator cell = mesh_generator.triangulation.begin_active(),
														   endc = mesh_generator.triangulation.end();





		for (; cell != endc ; cell++)
			for (unsigned int face = 0 ;face < GeometryInfo<dim>::faces_per_cell  ; face++)
				if (cell->face(face)->at_boundary())
				{
					double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);
                    
                    if (x_cord == constants.constants.xr)
                      EXPECT_EQ(cell->face(face)->boundary_id(),2);

                    // left edge
                    if (x_cord == constants.constants.xl)
                      EXPECT_EQ(cell->face(face)->boundary_id(),0);

                    // This is the first wall
                    if (y_cord == constants.constants.yb)
                      EXPECT_EQ(cell->face(face)->boundary_id(),1);

                    // top edge, This is the second wall
                    if (y_cord == constants.constants.yt)
                      EXPECT_EQ(cell->face(face)->boundary_id(),3);
				}



	}
}