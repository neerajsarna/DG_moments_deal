// in the following 
namespace PeriodicityHandler
{
	using namespace dealii;

	template<int dim>
	class 
	Base_PeriodicityHandler
	{
		public:
			Base_PeriodicityHandler(const double xl,const double xr);

		const double xl;
		const double xr;

        // stores the y coord of the cell center and cell iterator
        std::map<double, typename DoFHandler<dim>::cell_iterator> xminus1_set;
        std::map<double, typename DoFHandler<dim>::cell_iterator> xplus1_set;

        // stores the y coord of the cell center and the local face number
        std::map<double, unsigned int> set_xl_face;
        std::map<double, unsigned int> set_xr_face;

		std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vector;

		void develop_periodic_faces(DoFHandler<dim> &dof_handler);

        void add_periodic_sparsity(DynamicSparsityPattern &dsp);

        void divide_periodicity();

        // given the x coord of the face and the y cordinate of the center, the following argument,
        // return the iterator to the neighbouring cell.
        typename DoFHandler<dim>::cell_iterator get_periodic_neighbor(const double xcord_face,
                                                					 const double ycord_center) const;

        // given the xcord face and the ycord of the center of the cell, it provides the face number of
        // the neighbour
        unsigned int get_periodic_neighbor_face(const double xcord_face,
                                                 const double ycord_center)const;
	};

	template<int dim>
	Base_PeriodicityHandler<dim>::Base_PeriodicityHandler(const double xl,const double xr)
	:
	xl(xl),xr(xr)
	{;}

	template<int dim>
	void 
	Base_PeriodicityHandler<dim>::develop_periodic_faces(DoFHandler<dim> &dof_handler)
	{
		// collect the periodic faces of the mesh
		// The ID of the periodic faces is 100 and 101
		GridTools::collect_periodic_faces(dof_handler, 100,101,0,periodicity_vector);
	}

	// in the following routine we add entries to a given sparsity pattern so as to accommodate periodicity
	template<int dim>
	void
	Base_PeriodicityHandler<dim>::add_periodic_sparsity(DynamicSparsityPattern &dsp)
	{

			Assert(periodicity_vector.size() != 0,ExcNotInitialized());

      // we first compute the total number of pairs which deal has found
		    const unsigned int n_pairs = periodicity_vector.size();

      // the global degrees of freedom of cell1
		    std::vector<types::global_dof_index> local_dof_cell1;

      // the global degrees of freedom of cell 2
    		std::vector<types::global_dof_index> local_dof_cell2;

  
    		for (unsigned int i = 0 ; i < n_pairs ; i ++)
    		{

          // iterator for cell 1 in the pair
    			const typename DoFHandler<dim>::cell_iterator cell_1 = periodicity_vector[i].cell[0];

          // iterator for cell 2 in the pair
    			const typename DoFHandler<dim>::cell_iterator cell_2 = periodicity_vector[i].cell[1];

          // number of dofs in cell1
    			const unsigned int dof_cell1 = cell_1->get_fe().dofs_per_cell;

          // number of dofs in cell2
    			const unsigned int dof_cell2 = cell_2->get_fe().dofs_per_cell;

          // resize the desired vectors
    			local_dof_cell1.resize(dof_cell1);
    			local_dof_cell2.resize(dof_cell2);


    			cell_1->get_dof_indices(local_dof_cell1);
    			cell_2->get_dof_indices(local_dof_cell2);

          // add new locations in the dynamic sparsity pattern
    			for (unsigned int i = 0 ; i < dof_cell1 ; i++)
    				for (unsigned int j = 0 ; j < dof_cell2 ; j++)
    				{
    					dsp.add(local_dof_cell1[i],local_dof_cell2[j]);
    					dsp.add(local_dof_cell2[j],local_dof_cell1[i]);
    				}

      }


	}

  // The function develop_periodic_faces generates the periodicity vector, for details see dealii documentation. 
  // In divide periodicity we divide the periodicity vector depending upon the x and y coordinates of the cell. In xplus1_set
  // we store the cells which have a face at x = 1. In the pair xplus1_set we store the ycoord and the cell iterator. In xminus1_set we do
  // the same thing but for xplus1_set.
	template<int dim>
	void
	Base_PeriodicityHandler<dim>::divide_periodicity()
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

        std::pair<double,typename DoFHandler<dim>::cell_iterator> pair_plus;
        std::pair<double,typename DoFHandler<dim>::cell_iterator> pair_minus;

        std::pair<double,unsigned int> pair_plus_face;
        std::pair<double,unsigned int> pair_minus_face;

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
            if (x_cord1 == xr)       //Right edge division on the basis of cell center
            {
              cells_cord_plus1 ++;
              pair_plus = std::make_pair(y_center_cord1,cell_1);
              pair_plus_face = std::make_pair(y_center_cord1,face);

              xplus1_set.insert(pair_plus);
              set_xr_face.insert(pair_plus_face);
            }

            if (x_cord1 == xl)        // Left edge
            {
              cells_cord_minus1 ++;
              pair_minus = std::make_pair(y_center_cord1,cell_1);
              pair_minus_face = std::make_pair(y_center_cord1,face);

              xminus1_set.insert(pair_minus);
              set_xl_face.insert(pair_minus_face);
            }
          }

          // computation for the second cell
          if (face_itr_2->at_boundary())
          {
            if (x_cord2 == xr)       // division on the basis of cell center
            {
              cells_cord_plus1 ++;
              pair_plus = std::make_pair(y_center_cord2,cell_2);
              pair_plus_face = std::make_pair(y_center_cord2,face);

              xplus1_set.insert(pair_plus);
              set_xr_face.insert(pair_plus_face);
            }

            if (x_cord2 == xl)        // 
            {
              cells_cord_minus1 ++;
              pair_minus = std::make_pair(y_center_cord2,cell_2);
              pair_minus_face = std::make_pair(y_center_cord2,face);

              xminus1_set.insert(pair_minus);
              set_xl_face.insert(pair_minus_face);
            }
          }
        }
      }

      Assert(xplus1_set.size() == xminus1_set.size(),ExcMessage("unequal splitting"));

	}

	template<int dim>
    typename DoFHandler<dim>::cell_iterator
    Base_PeriodicityHandler<dim>::get_periodic_neighbor(const double xcord_face,const double ycord_center) const
    {
    // if cell at the left then find cell at plus1
      if (fabs(xcord_face - xl) < 1e-10) 
      {
        auto it = xplus1_set.find(ycord_center);
        Assert(it != xplus1_set.end(),
              ExcMessage("cannot find the neighbor on the right boundary"));

        return(it->second);
      }

      if (fabs(xcord_face - xr) < 1e-10)
      {
        auto it = xminus1_set.find(ycord_center);
        Assert(it != xminus1_set.end(),
              ExcMessage("cannot find the neighbor on the left boundary"));
        return(it->second);
      }

      Assert(1 == 0,ExcMessage("could not find periodic neighbor"));

      // the following statment has been added to avoid warning. The code will anyhow crash
      // due to the above statement.
      return xplus1_set.begin()->second;

    }

    template<int dim>
    unsigned int 
    Base_PeriodicityHandler<dim>::get_periodic_neighbor_face(const double xcord_face,
    														 const double ycord_center)const
    {
    // if cell at the left then find neighbor at plus 1
      if (xcord_face == xl)
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

}