  template<int num_flux,int dim>
  void 
  Solver_DG<num_flux,dim>::add_periodic_sparsity(DynamicSparsityPattern &dsp,
            vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>> &periodicity_vector)
{
    const unsigned int n_pairs = periodicity_vector.size();
    vector<types::global_dof_index> local_dof_cell1;
    vector<types::global_dof_index> local_dof_cell2;

  
    for (unsigned int i = 0 ; i < n_pairs ; i ++)
    {
      const typename DoFHandler<dim>::cell_iterator cell_1 = periodicity_vector[i].cell[0];
      const typename DoFHandler<dim>::cell_iterator cell_2 = periodicity_vector[i].cell[1];

      const unsigned int dof_cell1 = cell_1->get_fe().dofs_per_cell;
      const unsigned int dof_cell2 = cell_2->get_fe().dofs_per_cell;

      local_dof_cell1.resize(dof_cell1);
      local_dof_cell2.resize(dof_cell2);

      cell_1->get_dof_indices(local_dof_cell1);
      cell_2->get_dof_indices(local_dof_cell2);


      for (unsigned int i = 0 ; i < dof_cell1 ; i++)
        for (unsigned int j = 0 ; j < dof_cell2 ; j++)
        {
          dsp.add(local_dof_cell1[i],local_dof_cell2[j]);
          dsp.add(local_dof_cell2[j],local_dof_cell1[i]);
        }

      }


}

  template<int num_flux,int dim>
  void 
  Solver_DG<num_flux,dim>
  ::divide_periodicity(vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > &periodicity_vector,
             map<double, typename DoFHandler<dim>::cell_iterator> &xminus1_set,
             map<double, typename DoFHandler<dim>::cell_iterator> &xplus1_set,
             map<double, unsigned int> &set_xl_face,
             map<double, unsigned int> &set_xr_face,
             const double xl,
             const double xr)
{
      Assert(xl < xr,ExcMessage("provided min of boundary greater than provided max"));
      Assert(periodicity_vector.size() !=0,ExcNotInitialized());

      unsigned int cells_cord_plus1 = 0;
      unsigned int cells_cord_minus1 = 0;

      const unsigned int n_pairs = periodicity_vector.size();

      for (unsigned int i = 0 ; i < n_pairs ; i ++)
      {
        const typename DoFHandler<dim>::cell_iterator cell_1 = periodicity_vector[i].cell[0];
        const typename DoFHandler<dim>::cell_iterator cell_2 = periodicity_vector[i].cell[1];

        pair<double,typename DoFHandler<dim>::cell_iterator> pair_plus;
        pair<double,typename DoFHandler<dim>::cell_iterator> pair_minus;

        pair<double,unsigned int> pair_plus_face;
        pair<double,unsigned int> pair_minus_face;

        for (unsigned int face = 0 ; 
          face < GeometryInfo<dim>::faces_per_cell ;
          face++)
        {
          const typename DoFHandler<dim>::face_iterator face_itr_1 = cell_1->face(face);
          const typename DoFHandler<dim>::face_iterator face_itr_2 = cell_2->face(face);

          const double x_cord1 = face_itr_1->center()(0);
          const double x_cord2 = face_itr_2->center()(0);

          const double y_center_cord1 = cell_1->center()(1);
          const double y_center_cord2 = cell_2->center()(1);

          Assert(y_center_cord2 == y_center_cord1,
                 ExcMessage("periodic cells have different y coordinates"));



          // computation for the first cell
          if (face_itr_1->at_boundary())
          {
            if (x_cord1 == mesh_info.xr)       //Right edge division on the basis of cell center
            {
              cells_cord_plus1 ++;
              pair_plus = make_pair(y_center_cord1,cell_1);
              pair_plus_face = make_pair(y_center_cord1,face);

              xplus1_set.insert(pair_plus);
              set_xr_face.insert(pair_plus_face);
            }

            if (x_cord1 == mesh_info.xl)        // Left edge
            {
              cells_cord_minus1 ++;
              pair_minus = make_pair(y_center_cord1,cell_1);
              pair_minus_face = make_pair(y_center_cord1,face);

              xminus1_set.insert(pair_minus);
              set_xl_face.insert(pair_minus_face);
            }
          }

          // computation for the secod cell
          if (face_itr_2->at_boundary())
          {
            if (x_cord2 == mesh_info.xr)       // division on the basis of cell center
            {
              cells_cord_plus1 ++;
              pair_plus = make_pair(y_center_cord2,cell_2);
              pair_plus_face = make_pair(y_center_cord2,face);

              xplus1_set.insert(pair_plus);
              set_xr_face.insert(pair_plus_face);
            }

            if (x_cord2 == mesh_info.xl)        // 
            {
              cells_cord_minus1 ++;
              pair_minus = make_pair(y_center_cord2,cell_2);
              pair_minus_face = make_pair(y_center_cord2,face);

              xminus1_set.insert(pair_minus);
              set_xl_face.insert(pair_minus_face);
            }
          }
        }
      }

      Assert(xplus1_set.size() == xminus1_set.size(),ExcMessage("unequal splitting"));
}

template<int num_flux,int dim>
  typename DoFHandler<dim>::cell_iterator 
  Solver_DG<num_flux,dim>::get_periodic_neighbor(const double xcord_face,
                                                 const double ycord_center) const
  {
      // if cell at the left then find cell at plus1
      if (fabs(xcord_face - mesh_info.xl) < 1e-10) 
      {
        auto it = set_xr.find(ycord_center);
        Assert(it != set_xr.end(),
              ExcMessage("cannot find the neighbor on the right boundary"));

        return(it->second);
      }

      if (fabs(xcord_face - mesh_info.xr) < 1e-10)
      {
        auto it = set_xl.find(ycord_center);
        Assert(it != set_xl.end(),
              ExcMessage("cannot find the neighbor on the left boundary"));
        return(it->second);
      }

      Assert(1 == 0,ExcMessage("could not find periodic neighbor"));

      // the following statment has been added to avoid warning. The code will anyhow crash
      // due to the above statement.
      return set_xr.begin()->second;

  }

  template<int num_flux,int dim>
  unsigned int
  Solver_DG<num_flux,dim>::get_periodic_neighbor_face(const double xcord_face,
                                                 const double ycord_center) const
  {
      // if cell at the left then find cell at plus1
      if (xcord_face == mesh_info.xl)
      {
        auto it = set_xr_face.find(ycord_center);
        Assert(it != set_xr_face.end(),
              ExcMessage("cannot find the face number on the right boundary"));
        return(it->second);
      }

      else
      {
        auto it = set_xl_face.find(ycord_center);
        Assert(it != set_xr_face.end(),
              ExcMessage("cannot find the face number on the left boundary"));
        return(it->second);
      }

      Assert(1 == 0,ExcMessage("could not find periodic neighbor face"));
      return 0;

  }