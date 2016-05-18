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
    virtual void generate_mesh(Triangulation<dim> &triangulation,SphericalManifold<dim> &boundary,const string mesh_file_name)= 0;      
  };

  template<int dim> class mesh_generation:public Base_MeshHandler<dim>
  {
  public:
    enum Meshing_Options {read_msh, generate_internal};
    mesh_generation(const enum Meshing_Options mesh_options);
    virtual void print_mesh_info(const Triangulation<dim> &tria) const;
    virtual void print_grid(const Triangulation<dim> &tria,const string filename) const;
    virtual void generate_mesh(Triangulation<dim> &triangulation,SphericalManifold<dim> &boundary,const string mesh_file_name);            
    Meshing_Options mesh_options;
  private:
    GridIn<dim> gridin;
  };

  template<int dim> mesh_generation<dim>::mesh_generation(const enum Meshing_Options mesh_options)
  :
  mesh_options(mesh_options)
  {

  }

  template<int dim> void mesh_generation<dim>::generate_mesh(Triangulation<dim> &triangulation,SphericalManifold<dim> &boundary,const string mesh_file_name)
  {
    triangulation.clear();

    const Point<dim> center (0,0);
    const double inner_radius = 0.5,
    outer_radius = 2.0;

    switch (mesh_options)
    {
      case read_msh:
      {
        gridin.attach_triangulation(triangulation);
        ifstream f(mesh_file_name);
        gridin.read_msh(f);
        triangulation.set_all_manifold_ids_on_boundary(0);
        triangulation.set_manifold(0,boundary);
        break;
      }

      case generate_internal:
      {
        GridGenerator::hyper_shell (triangulation,
          center, inner_radius, outer_radius,
          10);

        triangulation.set_all_manifold_ids_on_boundary(0);
        triangulation.set_manifold(0,boundary);
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
