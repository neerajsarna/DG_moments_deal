namespace Mesh_Handler
{
  #include <map>
  using namespace dealii;
  using namespace std;

  template<int dim> class Base_MeshHandler
  {
  public:
    virtual void print_mesh_info(const Triangulation<dim> &tria) const = 0;
    virtual void print_grid(const Triangulation<dim> &tria,const string filename) const = 0;
    virtual void generate_mesh(Triangulation<dim> &triangulation,
                               SphericalManifold<dim> &boundary,
                               GridIn<dim> &gridin) const = 0;      

  };

  template<int dim> class mesh_generation:public Base_MeshHandler<dim>
  {
  public:
    mesh_generation(const mesh_data &mesh_info);
    virtual void print_mesh_info(const Triangulation<dim> &tria) const;
    virtual void print_grid(const Triangulation<dim> &tria,const string filename) const;
    virtual void generate_mesh(Triangulation<dim> &triangulation,
                              SphericalManifold<dim> &boundary,
                              GridIn<dim> &gridin)const ;     
    void set_periodic_bid(Triangulation<dim> &tria)const;
    const mesh_data mesh_info;

  };

  template<int dim> mesh_generation<dim>::mesh_generation(const mesh_data &mesh_info)
  :
  mesh_info(mesh_info)
  {

  }

  template<int dim> 
  void 
  mesh_generation<dim>
  ::set_periodic_bid(Triangulation<dim> &tria) const
  {
        typename Triangulation<dim>::cell_iterator cell = tria.begin(),
                               endc = tria.end();

        for (; cell != endc ; cell++)
        {
                for (unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // right edge
                    if (x_cord == mesh_info.xr)
                      cell->face(face)->set_boundary_id(2);

                    // left edge
                    if (x_cord == mesh_info.xl)
                      cell->face(face)->set_boundary_id(0);

                    // bottom edge
                    if (y_cord == mesh_info.yb)
                      cell->face(face)->set_boundary_id(1);

                    // top edge
                    if (y_cord == mesh_info.yt)
                      cell->face(face)->set_boundary_id(3);
                   }
        }
  }

  template<int dim> void mesh_generation<dim>::generate_mesh(Triangulation<dim> &triangulation,
                                                            SphericalManifold<dim> &boundary,
                                                            GridIn<dim> &gridin) const
  {
    triangulation.clear();

    const Point<dim> center (0,0);
    const double inner_radius = 0.5,
    outer_radius = 2.0;

    switch (mesh_info.mesh_options)
    {
      case read_msh:
      {
        gridin.attach_triangulation(triangulation);
        ifstream f(mesh_info.mesh_filename);
        gridin.read_msh(f);
        triangulation.set_all_manifold_ids_on_boundary(0);
        triangulation.set_manifold(0,boundary);
        break;
      }

      case generate_internal:
      {

        switch(mesh_info.mesh_type)
        {
          case ring:
          {
            GridGenerator::hyper_shell (triangulation,
              center, inner_radius, outer_radius,
              10);

            triangulation.set_all_manifold_ids_on_boundary(0);
            triangulation.set_manifold(0,boundary);
            break;
          }

          case periodic_square:
          {
            const unsigned int division_per_dim = 10;
            const double xl = mesh_info.xl;
            const double xr = mesh_info.xr;

            GridGenerator::subdivided_hyper_cube  (triangulation,
                                                    division_per_dim,
                                                    mesh_info.xl,
                                                    mesh_info.xr);

            set_periodic_bid(triangulation);
            break;
          }
        }


        break;
      }

    }

  }

  template<int dim> void mesh_generation<dim>::print_grid(const Triangulation<dim> &triangulation,const string filename) const
    {
      ofstream out (filename.c_str());
      GridOut grid_out;
      grid_out.write_eps (triangulation, out);
    } 

  template<int dim> void mesh_generation<dim>::print_mesh_info(const Triangulation<dim> &triangulation) const
  {
    cout << "****************Mesh Info*************" << endl;

    string mesh_info;
    mesh_info = " #dim " + to_string(dim) +  ", #Cells " + to_string(triangulation.n_active_cells());


    map<unsigned int, unsigned int> boundary_count;
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
    
    for (map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
     it!=boundary_count.end();
     ++it)
     mesh_info += to_string(it->first) + "("+to_string(it->second) +"times)";

   if (triangulation.has_hanging_nodes())
      mesh_info += ", hanging nodes: True";
   else
      mesh_info += ", hanging nodes: False";

   mesh_info += ", vertices: " + to_string(triangulation.n_vertices());

   cout << mesh_info << endl;
   cout << "************************************" << endl ;
 }
}
