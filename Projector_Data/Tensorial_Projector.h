
// In the following file we store the tensorial projector for a symmetric system. That is a normal projector * S^{-0.5}
 template<>
 void 
 Base_TensorInfo<2>
 ::reinit_local(const double nx,const double ny,const double nz,
                 const double tx,const double ty,const double tz,
                 const double rx,const double ry,const double rz,
               std::vector<projector_data> &tensor_project)
 {
   // first we allocate the required memory
   tensor_project.resize(max_tensorial_degree + 1);

   // now we allocate memory for every individual projector
   allocate_tensor_memory(tensor_project);

   AssertDimension(tensor_project.size(),max_tensorial_degree + 1);

   // we first check whether P has been properly allocated or not
   for (unsigned int i = 0 ; i < max_tensorial_degree ; i++)
   {
     AssertDimension(tensor_project[i].P.rows(),components(i));
     AssertDimension(tensor_project[i].P.cols(),components(i));
   }

    unsigned int degree = 0;

   tensor_project[0].P << 1;

    degree = 1;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
   tensor_project[degree].P << nx,ny,tx,ty;  

    degree = 2;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
   tensor_project[degree].P << pow(nx,2),2*nx*ny,pow(ny,2),nx*tx,ny*tx + nx*ty,ny*ty,pow(tx,2),2*tx*ty,pow(ty,2); // end of the third row

    degree = 3;
   AssertDimension(tensor_project[degree].P.rows(),components(degree));
   tensor_project[degree].P << pow(nx,3),3*ny*pow(nx,2),3*nx*pow(ny,2),pow(ny,3),tx*pow(nx,2),
                              2*nx*ny*tx + ty*pow(nx,2),2*nx*ny*ty + tx*pow(ny,2),ty*pow(ny,2),nx*pow(tx,2),
                              2*nx*tx*ty + ny*pow(tx,2),2*ny*tx*ty + nx*pow(ty,2),ny*pow(ty,2),pow(tx,3),
                              3*ty*pow(tx,2),3*tx*pow(ty,2),pow(ty,3); // fourth row


    degree = 4;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
       tensor_project[degree].P <<pow(nx,4),4*ny*pow(nx,3),6*pow(nx,2)*pow(ny,2),4*nx*pow(ny,3),pow(ny,4),
   tx*pow(nx,3),3*ny*tx*pow(nx,2) + ty*pow(nx,3),3*ny*ty*pow(nx,2) + 3*nx*tx*pow(ny,2),
   3*nx*ty*pow(ny,2) + tx*pow(ny,3),ty*pow(ny,3),pow(nx,2)*pow(tx,2),
   2*tx*ty*pow(nx,2) + 2*nx*ny*pow(tx,2),
   4*nx*ny*tx*ty + pow(ny,2)*pow(tx,2) + pow(nx,2)*pow(ty,2),
   2*tx*ty*pow(ny,2) + 2*nx*ny*pow(ty,2),pow(ny,2)*pow(ty,2),nx*pow(tx,3),
   3*nx*ty*pow(tx,2) + ny*pow(tx,3),3*ny*ty*pow(tx,2) + 3*nx*tx*pow(ty,2),
   3*ny*tx*pow(ty,2) + nx*pow(ty,3),ny*pow(ty,3),pow(tx,4),4*ty*pow(tx,3),
   6*pow(tx,2)*pow(ty,2),4*tx*pow(ty,3),pow(ty,4);

    degree = 5;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
     tensor_project[degree].P <<pow(nx,5),5*ny*pow(nx,4),10*pow(nx,3)*pow(ny,2),10*pow(nx,2)*pow(ny,3),
   5*nx*pow(ny,4),pow(ny,5),tx*pow(nx,4),4*ny*tx*pow(nx,3) + ty*pow(nx,4),
   4*ny*ty*pow(nx,3) + 6*tx*pow(nx,2)*pow(ny,2),
   6*ty*pow(nx,2)*pow(ny,2) + 4*nx*tx*pow(ny,3),4*nx*ty*pow(ny,3) + tx*pow(ny,4),
   ty*pow(ny,4),pow(nx,3)*pow(tx,2),2*tx*ty*pow(nx,3) + 3*ny*pow(nx,2)*pow(tx,2),
   6*ny*tx*ty*pow(nx,2) + 3*nx*pow(ny,2)*pow(tx,2) + pow(nx,3)*pow(ty,2),
   6*nx*tx*ty*pow(ny,2) + pow(ny,3)*pow(tx,2) + 3*ny*pow(nx,2)*pow(ty,2),
   2*tx*ty*pow(ny,3) + 3*nx*pow(ny,2)*pow(ty,2),pow(ny,3)*pow(ty,2),pow(nx,2)*pow(tx,3),
   3*ty*pow(nx,2)*pow(tx,2) + 2*nx*ny*pow(tx,3),
   6*nx*ny*ty*pow(tx,2) + pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,2)*pow(ty,2),
   3*ty*pow(ny,2)*pow(tx,2) + 6*nx*ny*tx*pow(ty,2) + pow(nx,2)*pow(ty,3),
   3*tx*pow(ny,2)*pow(ty,2) + 2*nx*ny*pow(ty,3),pow(ny,2)*pow(ty,3),nx*pow(tx,4),
   4*nx*ty*pow(tx,3) + ny*pow(tx,4),4*ny*ty*pow(tx,3) + 6*nx*pow(tx,2)*pow(ty,2),
   6*ny*pow(tx,2)*pow(ty,2) + 4*nx*tx*pow(ty,3),4*ny*tx*pow(ty,3) + nx*pow(ty,4),
   ny*pow(ty,4),pow(tx,5),5*ty*pow(tx,4),10*pow(tx,3)*pow(ty,2),10*pow(tx,2)*pow(ty,3),
   5*tx*pow(ty,4),pow(ty,5);

    degree = 6;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P <<pow(nx,6),6*ny*pow(nx,5),15*pow(nx,4)*pow(ny,2),20*pow(nx,3)*pow(ny,3),
   15*pow(nx,2)*pow(ny,4),6*nx*pow(ny,5),pow(ny,6),tx*pow(nx,5),
   5*ny*tx*pow(nx,4) + ty*pow(nx,5),5*ny*ty*pow(nx,4) + 10*tx*pow(nx,3)*pow(ny,2),
   10*ty*pow(nx,3)*pow(ny,2) + 10*tx*pow(nx,2)*pow(ny,3),
   10*ty*pow(nx,2)*pow(ny,3) + 5*nx*tx*pow(ny,4),5*nx*ty*pow(ny,4) + tx*pow(ny,5),
   ty*pow(ny,5),pow(nx,4)*pow(tx,2),2*tx*ty*pow(nx,4) + 4*ny*pow(nx,3)*pow(tx,2),
   8*ny*tx*ty*pow(nx,3) + 6*pow(nx,2)*pow(ny,2)*pow(tx,2) + pow(nx,4)*pow(ty,2),
   12*tx*ty*pow(nx,2)*pow(ny,2) + 4*nx*pow(ny,3)*pow(tx,2) + 4*ny*pow(nx,3)*pow(ty,2),
   8*nx*tx*ty*pow(ny,3) + pow(ny,4)*pow(tx,2) + 6*pow(nx,2)*pow(ny,2)*pow(ty,2),
   2*tx*ty*pow(ny,4) + 4*nx*pow(ny,3)*pow(ty,2),pow(ny,4)*pow(ty,2),pow(nx,3)*pow(tx,3),
   3*ty*pow(nx,3)*pow(tx,2) + 3*ny*pow(nx,2)*pow(tx,3),
   9*ny*ty*pow(nx,2)*pow(tx,2) + 3*nx*pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,3)*pow(ty,2),
   9*nx*ty*pow(ny,2)*pow(tx,2) + pow(ny,3)*pow(tx,3) + 9*ny*tx*pow(nx,2)*pow(ty,2) + 
    pow(nx,3)*pow(ty,3),3*ty*pow(ny,3)*pow(tx,2) + 9*nx*tx*pow(ny,2)*pow(ty,2) + 
    3*ny*pow(nx,2)*pow(ty,3),3*tx*pow(ny,3)*pow(ty,2) + 3*nx*pow(ny,2)*pow(ty,3),
   pow(ny,3)*pow(ty,3),pow(nx,2)*pow(tx,4),4*ty*pow(nx,2)*pow(tx,3) + 2*nx*ny*pow(tx,4),
   8*nx*ny*ty*pow(tx,3) + pow(ny,2)*pow(tx,4) + 6*pow(nx,2)*pow(tx,2)*pow(ty,2),
   4*ty*pow(ny,2)*pow(tx,3) + 12*nx*ny*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,2)*pow(ty,3),
   6*pow(ny,2)*pow(tx,2)*pow(ty,2) + 8*nx*ny*tx*pow(ty,3) + pow(nx,2)*pow(ty,4),
   4*tx*pow(ny,2)*pow(ty,3) + 2*nx*ny*pow(ty,4),pow(ny,2)*pow(ty,4),nx*pow(tx,5),
   5*nx*ty*pow(tx,4) + ny*pow(tx,5),5*ny*ty*pow(tx,4) + 10*nx*pow(tx,3)*pow(ty,2),
   10*ny*pow(tx,3)*pow(ty,2) + 10*nx*pow(tx,2)*pow(ty,3),
   10*ny*pow(tx,2)*pow(ty,3) + 5*nx*tx*pow(ty,4),5*ny*tx*pow(ty,4) + nx*pow(ty,5),
   ny*pow(ty,5),pow(tx,6),6*ty*pow(tx,5),15*pow(tx,4)*pow(ty,2),20*pow(tx,3)*pow(ty,3),
   15*pow(tx,2)*pow(ty,4),6*tx*pow(ty,5),pow(ty,6);

    degree = 7;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P << pow(nx,7),7*ny*pow(nx,6),21*pow(nx,5)*pow(ny,2),35*pow(nx,4)*pow(ny,3),
   35*pow(nx,3)*pow(ny,4),21*pow(nx,2)*pow(ny,5),7*nx*pow(ny,6),pow(ny,7),tx*pow(nx,6),
   6*ny*tx*pow(nx,5) + ty*pow(nx,6),6*ny*ty*pow(nx,5) + 15*tx*pow(nx,4)*pow(ny,2),
   15*ty*pow(nx,4)*pow(ny,2) + 20*tx*pow(nx,3)*pow(ny,3),
   20*ty*pow(nx,3)*pow(ny,3) + 15*tx*pow(nx,2)*pow(ny,4),
   15*ty*pow(nx,2)*pow(ny,4) + 6*nx*tx*pow(ny,5),6*nx*ty*pow(ny,5) + tx*pow(ny,6),
   ty*pow(ny,6),pow(nx,5)*pow(tx,2),2*tx*ty*pow(nx,5) + 5*ny*pow(nx,4)*pow(tx,2),
   10*ny*tx*ty*pow(nx,4) + 10*pow(nx,3)*pow(ny,2)*pow(tx,2) + pow(nx,5)*pow(ty,2),
   20*tx*ty*pow(nx,3)*pow(ny,2) + 10*pow(nx,2)*pow(ny,3)*pow(tx,2) + 
    5*ny*pow(nx,4)*pow(ty,2),20*tx*ty*pow(nx,2)*pow(ny,3) + 5*nx*pow(ny,4)*pow(tx,2) + 
    10*pow(nx,3)*pow(ny,2)*pow(ty,2),
   10*nx*tx*ty*pow(ny,4) + pow(ny,5)*pow(tx,2) + 10*pow(nx,2)*pow(ny,3)*pow(ty,2),
   2*tx*ty*pow(ny,5) + 5*nx*pow(ny,4)*pow(ty,2),pow(ny,5)*pow(ty,2),pow(nx,4)*pow(tx,3),
   3*ty*pow(nx,4)*pow(tx,2) + 4*ny*pow(nx,3)*pow(tx,3),
   12*ny*ty*pow(nx,3)*pow(tx,2) + 6*pow(nx,2)*pow(ny,2)*pow(tx,3) + 
    3*tx*pow(nx,4)*pow(ty,2),18*ty*pow(nx,2)*pow(ny,2)*pow(tx,2) + 
    4*nx*pow(ny,3)*pow(tx,3) + 12*ny*tx*pow(nx,3)*pow(ty,2) + pow(nx,4)*pow(ty,3),
   12*nx*ty*pow(ny,3)*pow(tx,2) + pow(ny,4)*pow(tx,3) + 
    18*tx*pow(nx,2)*pow(ny,2)*pow(ty,2) + 4*ny*pow(nx,3)*pow(ty,3),
   3*ty*pow(ny,4)*pow(tx,2) + 12*nx*tx*pow(ny,3)*pow(ty,2) + 
    6*pow(nx,2)*pow(ny,2)*pow(ty,3),3*tx*pow(ny,4)*pow(ty,2) + 4*nx*pow(ny,3)*pow(ty,3),
   pow(ny,4)*pow(ty,3),pow(nx,3)*pow(tx,4),
   4*ty*pow(nx,3)*pow(tx,3) + 3*ny*pow(nx,2)*pow(tx,4),
   12*ny*ty*pow(nx,2)*pow(tx,3) + 3*nx*pow(ny,2)*pow(tx,4) + 
    6*pow(nx,3)*pow(tx,2)*pow(ty,2),
   12*nx*ty*pow(ny,2)*pow(tx,3) + pow(ny,3)*pow(tx,4) + 
    18*ny*pow(nx,2)*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,3)*pow(ty,3),
   4*ty*pow(ny,3)*pow(tx,3) + 18*nx*pow(ny,2)*pow(tx,2)*pow(ty,2) + 
    12*ny*tx*pow(nx,2)*pow(ty,3) + pow(nx,3)*pow(ty,4),
   6*pow(ny,3)*pow(tx,2)*pow(ty,2) + 12*nx*tx*pow(ny,2)*pow(ty,3) + 
    3*ny*pow(nx,2)*pow(ty,4),4*tx*pow(ny,3)*pow(ty,3) + 3*nx*pow(ny,2)*pow(ty,4),
   pow(ny,3)*pow(ty,4),pow(nx,2)*pow(tx,5),5*ty*pow(nx,2)*pow(tx,4) + 2*nx*ny*pow(tx,5),
   10*nx*ny*ty*pow(tx,4) + pow(ny,2)*pow(tx,5) + 10*pow(nx,2)*pow(tx,3)*pow(ty,2),
   5*ty*pow(ny,2)*pow(tx,4) + 20*nx*ny*pow(tx,3)*pow(ty,2) + 
    10*pow(nx,2)*pow(tx,2)*pow(ty,3),
   10*pow(ny,2)*pow(tx,3)*pow(ty,2) + 20*nx*ny*pow(tx,2)*pow(ty,3) + 
    5*tx*pow(nx,2)*pow(ty,4),10*pow(ny,2)*pow(tx,2)*pow(ty,3) + 10*nx*ny*tx*pow(ty,4) + 
    pow(nx,2)*pow(ty,5),5*tx*pow(ny,2)*pow(ty,4) + 2*nx*ny*pow(ty,5),
   pow(ny,2)*pow(ty,5),nx*pow(tx,6),6*nx*ty*pow(tx,5) + ny*pow(tx,6),
   6*ny*ty*pow(tx,5) + 15*nx*pow(tx,4)*pow(ty,2),
   15*ny*pow(tx,4)*pow(ty,2) + 20*nx*pow(tx,3)*pow(ty,3),
   20*ny*pow(tx,3)*pow(ty,3) + 15*nx*pow(tx,2)*pow(ty,4),
   15*ny*pow(tx,2)*pow(ty,4) + 6*nx*tx*pow(ty,5),6*ny*tx*pow(ty,5) + nx*pow(ty,6),
   ny*pow(ty,6),pow(tx,7),7*ty*pow(tx,6),21*pow(tx,5)*pow(ty,2),35*pow(tx,4)*pow(ty,3),
   35*pow(tx,3)*pow(ty,4),21*pow(tx,2)*pow(ty,5),7*tx*pow(ty,6),pow(ty,7);


    degree = 8;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,8);

    tensor_project[degree].P << pow(nx,8),8*ny*pow(nx,7),28*pow(nx,6)*pow(ny,2),56*pow(nx,5)*pow(ny,3),70*pow(nx,4)*pow(ny,4),56*pow(nx,3)*pow(ny,5),28*pow(nx,2)*pow(ny,6),
   8*nx*pow(ny,7),pow(ny,8),tx*pow(nx,7),7*ny*tx*pow(nx,6) + ty*pow(nx,7),7*ny*ty*pow(nx,6) + 21*tx*pow(nx,5)*pow(ny,2),
   21*ty*pow(nx,5)*pow(ny,2) + 35*tx*pow(nx,4)*pow(ny,3),35*ty*pow(nx,4)*pow(ny,3) + 35*tx*pow(nx,3)*pow(ny,4),
   35*ty*pow(nx,3)*pow(ny,4) + 21*tx*pow(nx,2)*pow(ny,5),21*ty*pow(nx,2)*pow(ny,5) + 7*nx*tx*pow(ny,6),7*nx*ty*pow(ny,6) + tx*pow(ny,7),ty*pow(ny,7),
   pow(nx,6)*pow(tx,2),2*tx*ty*pow(nx,6) + 6*ny*pow(nx,5)*pow(tx,2),12*ny*tx*ty*pow(nx,5) + 15*pow(nx,4)*pow(ny,2)*pow(tx,2) + pow(nx,6)*pow(ty,2),
   30*tx*ty*pow(nx,4)*pow(ny,2) + 20*pow(nx,3)*pow(ny,3)*pow(tx,2) + 6*ny*pow(nx,5)*pow(ty,2),
   40*tx*ty*pow(nx,3)*pow(ny,3) + 15*pow(nx,2)*pow(ny,4)*pow(tx,2) + 15*pow(nx,4)*pow(ny,2)*pow(ty,2),
   30*tx*ty*pow(nx,2)*pow(ny,4) + 6*nx*pow(ny,5)*pow(tx,2) + 20*pow(nx,3)*pow(ny,3)*pow(ty,2),
   12*nx*tx*ty*pow(ny,5) + pow(ny,6)*pow(tx,2) + 15*pow(nx,2)*pow(ny,4)*pow(ty,2),2*tx*ty*pow(ny,6) + 6*nx*pow(ny,5)*pow(ty,2),pow(ny,6)*pow(ty,2),
   pow(nx,5)*pow(tx,3),3*ty*pow(nx,5)*pow(tx,2) + 5*ny*pow(nx,4)*pow(tx,3),
   15*ny*ty*pow(nx,4)*pow(tx,2) + 10*pow(nx,3)*pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,5)*pow(ty,2),
   30*ty*pow(nx,3)*pow(ny,2)*pow(tx,2) + 10*pow(nx,2)*pow(ny,3)*pow(tx,3) + 15*ny*tx*pow(nx,4)*pow(ty,2) + pow(nx,5)*pow(ty,3),
   30*ty*pow(nx,2)*pow(ny,3)*pow(tx,2) + 5*nx*pow(ny,4)*pow(tx,3) + 30*tx*pow(nx,3)*pow(ny,2)*pow(ty,2) + 5*ny*pow(nx,4)*pow(ty,3),
   15*nx*ty*pow(ny,4)*pow(tx,2) + pow(ny,5)*pow(tx,3) + 30*tx*pow(nx,2)*pow(ny,3)*pow(ty,2) + 10*pow(nx,3)*pow(ny,2)*pow(ty,3),
   3*ty*pow(ny,5)*pow(tx,2) + 15*nx*tx*pow(ny,4)*pow(ty,2) + 10*pow(nx,2)*pow(ny,3)*pow(ty,3),3*tx*pow(ny,5)*pow(ty,2) + 5*nx*pow(ny,4)*pow(ty,3),
   pow(ny,5)*pow(ty,3),pow(nx,4)*pow(tx,4),4*ty*pow(nx,4)*pow(tx,3) + 4*ny*pow(nx,3)*pow(tx,4),
   16*ny*ty*pow(nx,3)*pow(tx,3) + 6*pow(nx,2)*pow(ny,2)*pow(tx,4) + 6*pow(nx,4)*pow(tx,2)*pow(ty,2),
   24*ty*pow(nx,2)*pow(ny,2)*pow(tx,3) + 4*nx*pow(ny,3)*pow(tx,4) + 24*ny*pow(nx,3)*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,4)*pow(ty,3),
   16*nx*ty*pow(ny,3)*pow(tx,3) + pow(ny,4)*pow(tx,4) + 36*pow(nx,2)*pow(ny,2)*pow(tx,2)*pow(ty,2) + 16*ny*tx*pow(nx,3)*pow(ty,3) + 
    pow(nx,4)*pow(ty,4),4*ty*pow(ny,4)*pow(tx,3) + 24*nx*pow(ny,3)*pow(tx,2)*pow(ty,2) + 24*tx*pow(nx,2)*pow(ny,2)*pow(ty,3) + 
    4*ny*pow(nx,3)*pow(ty,4),6*pow(ny,4)*pow(tx,2)*pow(ty,2) + 16*nx*tx*pow(ny,3)*pow(ty,3) + 6*pow(nx,2)*pow(ny,2)*pow(ty,4),
   4*tx*pow(ny,4)*pow(ty,3) + 4*nx*pow(ny,3)*pow(ty,4),pow(ny,4)*pow(ty,4),pow(nx,3)*pow(tx,5),5*ty*pow(nx,3)*pow(tx,4) + 3*ny*pow(nx,2)*pow(tx,5),
   15*ny*ty*pow(nx,2)*pow(tx,4) + 3*nx*pow(ny,2)*pow(tx,5) + 10*pow(nx,3)*pow(tx,3)*pow(ty,2),
   15*nx*ty*pow(ny,2)*pow(tx,4) + pow(ny,3)*pow(tx,5) + 30*ny*pow(nx,2)*pow(tx,3)*pow(ty,2) + 10*pow(nx,3)*pow(tx,2)*pow(ty,3),
   5*ty*pow(ny,3)*pow(tx,4) + 30*nx*pow(ny,2)*pow(tx,3)*pow(ty,2) + 30*ny*pow(nx,2)*pow(tx,2)*pow(ty,3) + 5*tx*pow(nx,3)*pow(ty,4),
   10*pow(ny,3)*pow(tx,3)*pow(ty,2) + 30*nx*pow(ny,2)*pow(tx,2)*pow(ty,3) + 15*ny*tx*pow(nx,2)*pow(ty,4) + pow(nx,3)*pow(ty,5),
   10*pow(ny,3)*pow(tx,2)*pow(ty,3) + 15*nx*tx*pow(ny,2)*pow(ty,4) + 3*ny*pow(nx,2)*pow(ty,5),5*tx*pow(ny,3)*pow(ty,4) + 3*nx*pow(ny,2)*pow(ty,5),
   pow(ny,3)*pow(ty,5),pow(nx,2)*pow(tx,6),6*ty*pow(nx,2)*pow(tx,5) + 2*nx*ny*pow(tx,6),
   12*nx*ny*ty*pow(tx,5) + pow(ny,2)*pow(tx,6) + 15*pow(nx,2)*pow(tx,4)*pow(ty,2),
   6*ty*pow(ny,2)*pow(tx,5) + 30*nx*ny*pow(tx,4)*pow(ty,2) + 20*pow(nx,2)*pow(tx,3)*pow(ty,3),
   15*pow(ny,2)*pow(tx,4)*pow(ty,2) + 40*nx*ny*pow(tx,3)*pow(ty,3) + 15*pow(nx,2)*pow(tx,2)*pow(ty,4),
   20*pow(ny,2)*pow(tx,3)*pow(ty,3) + 30*nx*ny*pow(tx,2)*pow(ty,4) + 6*tx*pow(nx,2)*pow(ty,5),
   15*pow(ny,2)*pow(tx,2)*pow(ty,4) + 12*nx*ny*tx*pow(ty,5) + pow(nx,2)*pow(ty,6),6*tx*pow(ny,2)*pow(ty,5) + 2*nx*ny*pow(ty,6),pow(ny,2)*pow(ty,6),
   nx*pow(tx,7),7*nx*ty*pow(tx,6) + ny*pow(tx,7),7*ny*ty*pow(tx,6) + 21*nx*pow(tx,5)*pow(ty,2),21*ny*pow(tx,5)*pow(ty,2) + 35*nx*pow(tx,4)*pow(ty,3),
   35*ny*pow(tx,4)*pow(ty,3) + 35*nx*pow(tx,3)*pow(ty,4),35*ny*pow(tx,3)*pow(ty,4) + 21*nx*pow(tx,2)*pow(ty,5),
   21*ny*pow(tx,2)*pow(ty,5) + 7*nx*tx*pow(ty,6),7*ny*tx*pow(ty,6) + nx*pow(ty,7),ny*pow(ty,7),pow(tx,8),8*ty*pow(tx,7),28*pow(tx,6)*pow(ty,2),
   56*pow(tx,5)*pow(ty,3),70*pow(tx,4)*pow(ty,4),56*pow(tx,3)*pow(ty,5),28*pow(tx,2)*pow(ty,6),8*tx*pow(ty,7),pow(ty,8);

    degree = 9;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,9);

    tensor_project[degree].P << pow(nx,9),9*ny*pow(nx,8),36*pow(nx,7)*pow(ny,2),84*pow(nx,6)*pow(ny,3),126*pow(nx,5)*pow(ny,4),126*pow(nx,4)*pow(ny,5),84*pow(nx,3)*pow(ny,6),
   36*pow(nx,2)*pow(ny,7),9*nx*pow(ny,8),pow(ny,9),tx*pow(nx,8),8*ny*tx*pow(nx,7) + ty*pow(nx,8),8*ny*ty*pow(nx,7) + 28*tx*pow(nx,6)*pow(ny,2),
   28*ty*pow(nx,6)*pow(ny,2) + 56*tx*pow(nx,5)*pow(ny,3),56*ty*pow(nx,5)*pow(ny,3) + 70*tx*pow(nx,4)*pow(ny,4),
   70*ty*pow(nx,4)*pow(ny,4) + 56*tx*pow(nx,3)*pow(ny,5),56*ty*pow(nx,3)*pow(ny,5) + 28*tx*pow(nx,2)*pow(ny,6),
   28*ty*pow(nx,2)*pow(ny,6) + 8*nx*tx*pow(ny,7),8*nx*ty*pow(ny,7) + tx*pow(ny,8),ty*pow(ny,8),pow(nx,7)*pow(tx,2),
   2*tx*ty*pow(nx,7) + 7*ny*pow(nx,6)*pow(tx,2),14*ny*tx*ty*pow(nx,6) + 21*pow(nx,5)*pow(ny,2)*pow(tx,2) + pow(nx,7)*pow(ty,2),
   42*tx*ty*pow(nx,5)*pow(ny,2) + 35*pow(nx,4)*pow(ny,3)*pow(tx,2) + 7*ny*pow(nx,6)*pow(ty,2),
   70*tx*ty*pow(nx,4)*pow(ny,3) + 35*pow(nx,3)*pow(ny,4)*pow(tx,2) + 21*pow(nx,5)*pow(ny,2)*pow(ty,2),
   70*tx*ty*pow(nx,3)*pow(ny,4) + 21*pow(nx,2)*pow(ny,5)*pow(tx,2) + 35*pow(nx,4)*pow(ny,3)*pow(ty,2),
   42*tx*ty*pow(nx,2)*pow(ny,5) + 7*nx*pow(ny,6)*pow(tx,2) + 35*pow(nx,3)*pow(ny,4)*pow(ty,2),
   14*nx*tx*ty*pow(ny,6) + pow(ny,7)*pow(tx,2) + 21*pow(nx,2)*pow(ny,5)*pow(ty,2),2*tx*ty*pow(ny,7) + 7*nx*pow(ny,6)*pow(ty,2),pow(ny,7)*pow(ty,2),
   pow(nx,6)*pow(tx,3),3*ty*pow(nx,6)*pow(tx,2) + 6*ny*pow(nx,5)*pow(tx,3),
   18*ny*ty*pow(nx,5)*pow(tx,2) + 15*pow(nx,4)*pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,6)*pow(ty,2),
   45*ty*pow(nx,4)*pow(ny,2)*pow(tx,2) + 20*pow(nx,3)*pow(ny,3)*pow(tx,3) + 18*ny*tx*pow(nx,5)*pow(ty,2) + pow(nx,6)*pow(ty,3),
   60*ty*pow(nx,3)*pow(ny,3)*pow(tx,2) + 15*pow(nx,2)*pow(ny,4)*pow(tx,3) + 45*tx*pow(nx,4)*pow(ny,2)*pow(ty,2) + 6*ny*pow(nx,5)*pow(ty,3),
   45*ty*pow(nx,2)*pow(ny,4)*pow(tx,2) + 6*nx*pow(ny,5)*pow(tx,3) + 60*tx*pow(nx,3)*pow(ny,3)*pow(ty,2) + 15*pow(nx,4)*pow(ny,2)*pow(ty,3),
   18*nx*ty*pow(ny,5)*pow(tx,2) + pow(ny,6)*pow(tx,3) + 45*tx*pow(nx,2)*pow(ny,4)*pow(ty,2) + 20*pow(nx,3)*pow(ny,3)*pow(ty,3),
   3*ty*pow(ny,6)*pow(tx,2) + 18*nx*tx*pow(ny,5)*pow(ty,2) + 15*pow(nx,2)*pow(ny,4)*pow(ty,3),3*tx*pow(ny,6)*pow(ty,2) + 6*nx*pow(ny,5)*pow(ty,3),
   pow(ny,6)*pow(ty,3),pow(nx,5)*pow(tx,4),4*ty*pow(nx,5)*pow(tx,3) + 5*ny*pow(nx,4)*pow(tx,4),
   20*ny*ty*pow(nx,4)*pow(tx,3) + 10*pow(nx,3)*pow(ny,2)*pow(tx,4) + 6*pow(nx,5)*pow(tx,2)*pow(ty,2),
   40*ty*pow(nx,3)*pow(ny,2)*pow(tx,3) + 10*pow(nx,2)*pow(ny,3)*pow(tx,4) + 30*ny*pow(nx,4)*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,5)*pow(ty,3),
   40*ty*pow(nx,2)*pow(ny,3)*pow(tx,3) + 5*nx*pow(ny,4)*pow(tx,4) + 60*pow(nx,3)*pow(ny,2)*pow(tx,2)*pow(ty,2) + 20*ny*tx*pow(nx,4)*pow(ty,3) + 
    pow(nx,5)*pow(ty,4),20*nx*ty*pow(ny,4)*pow(tx,3) + pow(ny,5)*pow(tx,4) + 60*pow(nx,2)*pow(ny,3)*pow(tx,2)*pow(ty,2) + 
    40*tx*pow(nx,3)*pow(ny,2)*pow(ty,3) + 5*ny*pow(nx,4)*pow(ty,4),
   4*ty*pow(ny,5)*pow(tx,3) + 30*nx*pow(ny,4)*pow(tx,2)*pow(ty,2) + 40*tx*pow(nx,2)*pow(ny,3)*pow(ty,3) + 10*pow(nx,3)*pow(ny,2)*pow(ty,4),
   6*pow(ny,5)*pow(tx,2)*pow(ty,2) + 20*nx*tx*pow(ny,4)*pow(ty,3) + 10*pow(nx,2)*pow(ny,3)*pow(ty,4),
   4*tx*pow(ny,5)*pow(ty,3) + 5*nx*pow(ny,4)*pow(ty,4),pow(ny,5)*pow(ty,4),pow(nx,4)*pow(tx,5),5*ty*pow(nx,4)*pow(tx,4) + 4*ny*pow(nx,3)*pow(tx,5),
   20*ny*ty*pow(nx,3)*pow(tx,4) + 6*pow(nx,2)*pow(ny,2)*pow(tx,5) + 10*pow(nx,4)*pow(tx,3)*pow(ty,2),
   30*ty*pow(nx,2)*pow(ny,2)*pow(tx,4) + 4*nx*pow(ny,3)*pow(tx,5) + 40*ny*pow(nx,3)*pow(tx,3)*pow(ty,2) + 10*pow(nx,4)*pow(tx,2)*pow(ty,3),
   20*nx*ty*pow(ny,3)*pow(tx,4) + pow(ny,4)*pow(tx,5) + 60*pow(nx,2)*pow(ny,2)*pow(tx,3)*pow(ty,2) + 40*ny*pow(nx,3)*pow(tx,2)*pow(ty,3) + 
    5*tx*pow(nx,4)*pow(ty,4),5*ty*pow(ny,4)*pow(tx,4) + 40*nx*pow(ny,3)*pow(tx,3)*pow(ty,2) + 60*pow(nx,2)*pow(ny,2)*pow(tx,2)*pow(ty,3) + 
    20*ny*tx*pow(nx,3)*pow(ty,4) + pow(nx,4)*pow(ty,5),10*pow(ny,4)*pow(tx,3)*pow(ty,2) + 40*nx*pow(ny,3)*pow(tx,2)*pow(ty,3) + 
    30*tx*pow(nx,2)*pow(ny,2)*pow(ty,4) + 4*ny*pow(nx,3)*pow(ty,5),
   10*pow(ny,4)*pow(tx,2)*pow(ty,3) + 20*nx*tx*pow(ny,3)*pow(ty,4) + 6*pow(nx,2)*pow(ny,2)*pow(ty,5),
   5*tx*pow(ny,4)*pow(ty,4) + 4*nx*pow(ny,3)*pow(ty,5),pow(ny,4)*pow(ty,5),pow(nx,3)*pow(tx,6),6*ty*pow(nx,3)*pow(tx,5) + 3*ny*pow(nx,2)*pow(tx,6),
   18*ny*ty*pow(nx,2)*pow(tx,5) + 3*nx*pow(ny,2)*pow(tx,6) + 15*pow(nx,3)*pow(tx,4)*pow(ty,2),
   18*nx*ty*pow(ny,2)*pow(tx,5) + pow(ny,3)*pow(tx,6) + 45*ny*pow(nx,2)*pow(tx,4)*pow(ty,2) + 20*pow(nx,3)*pow(tx,3)*pow(ty,3),
   6*ty*pow(ny,3)*pow(tx,5) + 45*nx*pow(ny,2)*pow(tx,4)*pow(ty,2) + 60*ny*pow(nx,2)*pow(tx,3)*pow(ty,3) + 15*pow(nx,3)*pow(tx,2)*pow(ty,4),
   15*pow(ny,3)*pow(tx,4)*pow(ty,2) + 60*nx*pow(ny,2)*pow(tx,3)*pow(ty,3) + 45*ny*pow(nx,2)*pow(tx,2)*pow(ty,4) + 6*tx*pow(nx,3)*pow(ty,5),
   20*pow(ny,3)*pow(tx,3)*pow(ty,3) + 45*nx*pow(ny,2)*pow(tx,2)*pow(ty,4) + 18*ny*tx*pow(nx,2)*pow(ty,5) + pow(nx,3)*pow(ty,6),
   15*pow(ny,3)*pow(tx,2)*pow(ty,4) + 18*nx*tx*pow(ny,2)*pow(ty,5) + 3*ny*pow(nx,2)*pow(ty,6),6*tx*pow(ny,3)*pow(ty,5) + 3*nx*pow(ny,2)*pow(ty,6),
   pow(ny,3)*pow(ty,6),pow(nx,2)*pow(tx,7),7*ty*pow(nx,2)*pow(tx,6) + 2*nx*ny*pow(tx,7),
   14*nx*ny*ty*pow(tx,6) + pow(ny,2)*pow(tx,7) + 21*pow(nx,2)*pow(tx,5)*pow(ty,2),
   7*ty*pow(ny,2)*pow(tx,6) + 42*nx*ny*pow(tx,5)*pow(ty,2) + 35*pow(nx,2)*pow(tx,4)*pow(ty,3),
   21*pow(ny,2)*pow(tx,5)*pow(ty,2) + 70*nx*ny*pow(tx,4)*pow(ty,3) + 35*pow(nx,2)*pow(tx,3)*pow(ty,4),
   35*pow(ny,2)*pow(tx,4)*pow(ty,3) + 70*nx*ny*pow(tx,3)*pow(ty,4) + 21*pow(nx,2)*pow(tx,2)*pow(ty,5),
   35*pow(ny,2)*pow(tx,3)*pow(ty,4) + 42*nx*ny*pow(tx,2)*pow(ty,5) + 7*tx*pow(nx,2)*pow(ty,6),
   21*pow(ny,2)*pow(tx,2)*pow(ty,5) + 14*nx*ny*tx*pow(ty,6) + pow(nx,2)*pow(ty,7),7*tx*pow(ny,2)*pow(ty,6) + 2*nx*ny*pow(ty,7),pow(ny,2)*pow(ty,7),
   nx*pow(tx,8),8*nx*ty*pow(tx,7) + ny*pow(tx,8),8*ny*ty*pow(tx,7) + 28*nx*pow(tx,6)*pow(ty,2),28*ny*pow(tx,6)*pow(ty,2) + 56*nx*pow(tx,5)*pow(ty,3),
   56*ny*pow(tx,5)*pow(ty,3) + 70*nx*pow(tx,4)*pow(ty,4),70*ny*pow(tx,4)*pow(ty,4) + 56*nx*pow(tx,3)*pow(ty,5),
   56*ny*pow(tx,3)*pow(ty,5) + 28*nx*pow(tx,2)*pow(ty,6),28*ny*pow(tx,2)*pow(ty,6) + 8*nx*tx*pow(ty,7),8*ny*tx*pow(ty,7) + nx*pow(ty,8),ny*pow(ty,8),
   pow(tx,9),9*ty*pow(tx,8),36*pow(tx,7)*pow(ty,2),84*pow(tx,6)*pow(ty,3),126*pow(tx,5)*pow(ty,4),126*pow(tx,4)*pow(ty,5),84*pow(tx,3)*pow(ty,6),
   36*pow(tx,2)*pow(ty,7),9*tx*pow(ty,8),pow(ty,9);

    degree = 10;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,10);

    tensor_project[degree].P <<pow(nx,10),10*ny*pow(nx,9),45*pow(nx,8)*pow(ny,2),120*pow(nx,7)*pow(ny,3),210*pow(nx,6)*pow(ny,4),252*pow(nx,5)*pow(ny,5),
   210*pow(nx,4)*pow(ny,6),120*pow(nx,3)*pow(ny,7),45*pow(nx,2)*pow(ny,8),10*nx*pow(ny,9),pow(ny,10),tx*pow(nx,9),9*ny*tx*pow(nx,8) + ty*pow(nx,9),
   9*ny*ty*pow(nx,8) + 36*tx*pow(nx,7)*pow(ny,2),36*ty*pow(nx,7)*pow(ny,2) + 84*tx*pow(nx,6)*pow(ny,3),
   84*ty*pow(nx,6)*pow(ny,3) + 126*tx*pow(nx,5)*pow(ny,4),126*ty*pow(nx,5)*pow(ny,4) + 126*tx*pow(nx,4)*pow(ny,5),
   126*ty*pow(nx,4)*pow(ny,5) + 84*tx*pow(nx,3)*pow(ny,6),84*ty*pow(nx,3)*pow(ny,6) + 36*tx*pow(nx,2)*pow(ny,7),
   36*ty*pow(nx,2)*pow(ny,7) + 9*nx*tx*pow(ny,8),9*nx*ty*pow(ny,8) + tx*pow(ny,9),ty*pow(ny,9),pow(nx,8)*pow(tx,2),
   2*tx*ty*pow(nx,8) + 8*ny*pow(nx,7)*pow(tx,2),16*ny*tx*ty*pow(nx,7) + 28*pow(nx,6)*pow(ny,2)*pow(tx,2) + pow(nx,8)*pow(ty,2),
   56*tx*ty*pow(nx,6)*pow(ny,2) + 56*pow(nx,5)*pow(ny,3)*pow(tx,2) + 8*ny*pow(nx,7)*pow(ty,2),
   112*tx*ty*pow(nx,5)*pow(ny,3) + 70*pow(nx,4)*pow(ny,4)*pow(tx,2) + 28*pow(nx,6)*pow(ny,2)*pow(ty,2),
   140*tx*ty*pow(nx,4)*pow(ny,4) + 56*pow(nx,3)*pow(ny,5)*pow(tx,2) + 56*pow(nx,5)*pow(ny,3)*pow(ty,2),
   112*tx*ty*pow(nx,3)*pow(ny,5) + 28*pow(nx,2)*pow(ny,6)*pow(tx,2) + 70*pow(nx,4)*pow(ny,4)*pow(ty,2),
   56*tx*ty*pow(nx,2)*pow(ny,6) + 8*nx*pow(ny,7)*pow(tx,2) + 56*pow(nx,3)*pow(ny,5)*pow(ty,2),
   16*nx*tx*ty*pow(ny,7) + pow(ny,8)*pow(tx,2) + 28*pow(nx,2)*pow(ny,6)*pow(ty,2),2*tx*ty*pow(ny,8) + 8*nx*pow(ny,7)*pow(ty,2),pow(ny,8)*pow(ty,2),
   pow(nx,7)*pow(tx,3),3*ty*pow(nx,7)*pow(tx,2) + 7*ny*pow(nx,6)*pow(tx,3),
   21*ny*ty*pow(nx,6)*pow(tx,2) + 21*pow(nx,5)*pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,7)*pow(ty,2),
   63*ty*pow(nx,5)*pow(ny,2)*pow(tx,2) + 35*pow(nx,4)*pow(ny,3)*pow(tx,3) + 21*ny*tx*pow(nx,6)*pow(ty,2) + pow(nx,7)*pow(ty,3),
   105*ty*pow(nx,4)*pow(ny,3)*pow(tx,2) + 35*pow(nx,3)*pow(ny,4)*pow(tx,3) + 63*tx*pow(nx,5)*pow(ny,2)*pow(ty,2) + 7*ny*pow(nx,6)*pow(ty,3),
   105*ty*pow(nx,3)*pow(ny,4)*pow(tx,2) + 21*pow(nx,2)*pow(ny,5)*pow(tx,3) + 105*tx*pow(nx,4)*pow(ny,3)*pow(ty,2) + 21*pow(nx,5)*pow(ny,2)*pow(ty,3),
   63*ty*pow(nx,2)*pow(ny,5)*pow(tx,2) + 7*nx*pow(ny,6)*pow(tx,3) + 105*tx*pow(nx,3)*pow(ny,4)*pow(ty,2) + 35*pow(nx,4)*pow(ny,3)*pow(ty,3),
   21*nx*ty*pow(ny,6)*pow(tx,2) + pow(ny,7)*pow(tx,3) + 63*tx*pow(nx,2)*pow(ny,5)*pow(ty,2) + 35*pow(nx,3)*pow(ny,4)*pow(ty,3),
   3*ty*pow(ny,7)*pow(tx,2) + 21*nx*tx*pow(ny,6)*pow(ty,2) + 21*pow(nx,2)*pow(ny,5)*pow(ty,3),3*tx*pow(ny,7)*pow(ty,2) + 7*nx*pow(ny,6)*pow(ty,3),
   pow(ny,7)*pow(ty,3),pow(nx,6)*pow(tx,4),4*ty*pow(nx,6)*pow(tx,3) + 6*ny*pow(nx,5)*pow(tx,4),
   24*ny*ty*pow(nx,5)*pow(tx,3) + 15*pow(nx,4)*pow(ny,2)*pow(tx,4) + 6*pow(nx,6)*pow(tx,2)*pow(ty,2),
   60*ty*pow(nx,4)*pow(ny,2)*pow(tx,3) + 20*pow(nx,3)*pow(ny,3)*pow(tx,4) + 36*ny*pow(nx,5)*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,6)*pow(ty,3),
   80*ty*pow(nx,3)*pow(ny,3)*pow(tx,3) + 15*pow(nx,2)*pow(ny,4)*pow(tx,4) + 90*pow(nx,4)*pow(ny,2)*pow(tx,2)*pow(ty,2) + 
    24*ny*tx*pow(nx,5)*pow(ty,3) + pow(nx,6)*pow(ty,4),60*ty*pow(nx,2)*pow(ny,4)*pow(tx,3) + 6*nx*pow(ny,5)*pow(tx,4) + 
    120*pow(nx,3)*pow(ny,3)*pow(tx,2)*pow(ty,2) + 60*tx*pow(nx,4)*pow(ny,2)*pow(ty,3) + 6*ny*pow(nx,5)*pow(ty,4),
   24*nx*ty*pow(ny,5)*pow(tx,3) + pow(ny,6)*pow(tx,4) + 90*pow(nx,2)*pow(ny,4)*pow(tx,2)*pow(ty,2) + 80*tx*pow(nx,3)*pow(ny,3)*pow(ty,3) + 
    15*pow(nx,4)*pow(ny,2)*pow(ty,4),4*ty*pow(ny,6)*pow(tx,3) + 36*nx*pow(ny,5)*pow(tx,2)*pow(ty,2) + 60*tx*pow(nx,2)*pow(ny,4)*pow(ty,3) + 
    20*pow(nx,3)*pow(ny,3)*pow(ty,4),6*pow(ny,6)*pow(tx,2)*pow(ty,2) + 24*nx*tx*pow(ny,5)*pow(ty,3) + 15*pow(nx,2)*pow(ny,4)*pow(ty,4),
   4*tx*pow(ny,6)*pow(ty,3) + 6*nx*pow(ny,5)*pow(ty,4),pow(ny,6)*pow(ty,4),pow(nx,5)*pow(tx,5),5*ty*pow(nx,5)*pow(tx,4) + 5*ny*pow(nx,4)*pow(tx,5),
   25*ny*ty*pow(nx,4)*pow(tx,4) + 10*pow(nx,3)*pow(ny,2)*pow(tx,5) + 10*pow(nx,5)*pow(tx,3)*pow(ty,2),
   50*ty*pow(nx,3)*pow(ny,2)*pow(tx,4) + 10*pow(nx,2)*pow(ny,3)*pow(tx,5) + 50*ny*pow(nx,4)*pow(tx,3)*pow(ty,2) + 10*pow(nx,5)*pow(tx,2)*pow(ty,3),
   50*ty*pow(nx,2)*pow(ny,3)*pow(tx,4) + 5*nx*pow(ny,4)*pow(tx,5) + 100*pow(nx,3)*pow(ny,2)*pow(tx,3)*pow(ty,2) + 
    50*ny*pow(nx,4)*pow(tx,2)*pow(ty,3) + 5*tx*pow(nx,5)*pow(ty,4),
   25*nx*ty*pow(ny,4)*pow(tx,4) + pow(ny,5)*pow(tx,5) + 100*pow(nx,2)*pow(ny,3)*pow(tx,3)*pow(ty,2) + 100*pow(nx,3)*pow(ny,2)*pow(tx,2)*pow(ty,3) + 
    25*ny*tx*pow(nx,4)*pow(ty,4) + pow(nx,5)*pow(ty,5),5*ty*pow(ny,5)*pow(tx,4) + 50*nx*pow(ny,4)*pow(tx,3)*pow(ty,2) + 
    100*pow(nx,2)*pow(ny,3)*pow(tx,2)*pow(ty,3) + 50*tx*pow(nx,3)*pow(ny,2)*pow(ty,4) + 5*ny*pow(nx,4)*pow(ty,5),
   10*pow(ny,5)*pow(tx,3)*pow(ty,2) + 50*nx*pow(ny,4)*pow(tx,2)*pow(ty,3) + 50*tx*pow(nx,2)*pow(ny,3)*pow(ty,4) + 10*pow(nx,3)*pow(ny,2)*pow(ty,5),
   10*pow(ny,5)*pow(tx,2)*pow(ty,3) + 25*nx*tx*pow(ny,4)*pow(ty,4) + 10*pow(nx,2)*pow(ny,3)*pow(ty,5),
   5*tx*pow(ny,5)*pow(ty,4) + 5*nx*pow(ny,4)*pow(ty,5),pow(ny,5)*pow(ty,5),pow(nx,4)*pow(tx,6),6*ty*pow(nx,4)*pow(tx,5) + 4*ny*pow(nx,3)*pow(tx,6),
   24*ny*ty*pow(nx,3)*pow(tx,5) + 6*pow(nx,2)*pow(ny,2)*pow(tx,6) + 15*pow(nx,4)*pow(tx,4)*pow(ty,2),
   36*ty*pow(nx,2)*pow(ny,2)*pow(tx,5) + 4*nx*pow(ny,3)*pow(tx,6) + 60*ny*pow(nx,3)*pow(tx,4)*pow(ty,2) + 20*pow(nx,4)*pow(tx,3)*pow(ty,3),
   24*nx*ty*pow(ny,3)*pow(tx,5) + pow(ny,4)*pow(tx,6) + 90*pow(nx,2)*pow(ny,2)*pow(tx,4)*pow(ty,2) + 80*ny*pow(nx,3)*pow(tx,3)*pow(ty,3) + 
    15*pow(nx,4)*pow(tx,2)*pow(ty,4),6*ty*pow(ny,4)*pow(tx,5) + 60*nx*pow(ny,3)*pow(tx,4)*pow(ty,2) + 120*pow(nx,2)*pow(ny,2)*pow(tx,3)*pow(ty,3) + 
    60*ny*pow(nx,3)*pow(tx,2)*pow(ty,4) + 6*tx*pow(nx,4)*pow(ty,5),
   15*pow(ny,4)*pow(tx,4)*pow(ty,2) + 80*nx*pow(ny,3)*pow(tx,3)*pow(ty,3) + 90*pow(nx,2)*pow(ny,2)*pow(tx,2)*pow(ty,4) + 
    24*ny*tx*pow(nx,3)*pow(ty,5) + pow(nx,4)*pow(ty,6),20*pow(ny,4)*pow(tx,3)*pow(ty,3) + 60*nx*pow(ny,3)*pow(tx,2)*pow(ty,4) + 
    36*tx*pow(nx,2)*pow(ny,2)*pow(ty,5) + 4*ny*pow(nx,3)*pow(ty,6),
   15*pow(ny,4)*pow(tx,2)*pow(ty,4) + 24*nx*tx*pow(ny,3)*pow(ty,5) + 6*pow(nx,2)*pow(ny,2)*pow(ty,6),
   6*tx*pow(ny,4)*pow(ty,5) + 4*nx*pow(ny,3)*pow(ty,6),pow(ny,4)*pow(ty,6),pow(nx,3)*pow(tx,7),7*ty*pow(nx,3)*pow(tx,6) + 3*ny*pow(nx,2)*pow(tx,7),
   21*ny*ty*pow(nx,2)*pow(tx,6) + 3*nx*pow(ny,2)*pow(tx,7) + 21*pow(nx,3)*pow(tx,5)*pow(ty,2),
   21*nx*ty*pow(ny,2)*pow(tx,6) + pow(ny,3)*pow(tx,7) + 63*ny*pow(nx,2)*pow(tx,5)*pow(ty,2) + 35*pow(nx,3)*pow(tx,4)*pow(ty,3),
   7*ty*pow(ny,3)*pow(tx,6) + 63*nx*pow(ny,2)*pow(tx,5)*pow(ty,2) + 105*ny*pow(nx,2)*pow(tx,4)*pow(ty,3) + 35*pow(nx,3)*pow(tx,3)*pow(ty,4),
   21*pow(ny,3)*pow(tx,5)*pow(ty,2) + 105*nx*pow(ny,2)*pow(tx,4)*pow(ty,3) + 105*ny*pow(nx,2)*pow(tx,3)*pow(ty,4) + 21*pow(nx,3)*pow(tx,2)*pow(ty,5),
   35*pow(ny,3)*pow(tx,4)*pow(ty,3) + 105*nx*pow(ny,2)*pow(tx,3)*pow(ty,4) + 63*ny*pow(nx,2)*pow(tx,2)*pow(ty,5) + 7*tx*pow(nx,3)*pow(ty,6),
   35*pow(ny,3)*pow(tx,3)*pow(ty,4) + 63*nx*pow(ny,2)*pow(tx,2)*pow(ty,5) + 21*ny*tx*pow(nx,2)*pow(ty,6) + pow(nx,3)*pow(ty,7),
   21*pow(ny,3)*pow(tx,2)*pow(ty,5) + 21*nx*tx*pow(ny,2)*pow(ty,6) + 3*ny*pow(nx,2)*pow(ty,7),7*tx*pow(ny,3)*pow(ty,6) + 3*nx*pow(ny,2)*pow(ty,7),
   pow(ny,3)*pow(ty,7),pow(nx,2)*pow(tx,8),8*ty*pow(nx,2)*pow(tx,7) + 2*nx*ny*pow(tx,8),
   16*nx*ny*ty*pow(tx,7) + pow(ny,2)*pow(tx,8) + 28*pow(nx,2)*pow(tx,6)*pow(ty,2),
   8*ty*pow(ny,2)*pow(tx,7) + 56*nx*ny*pow(tx,6)*pow(ty,2) + 56*pow(nx,2)*pow(tx,5)*pow(ty,3),
   28*pow(ny,2)*pow(tx,6)*pow(ty,2) + 112*nx*ny*pow(tx,5)*pow(ty,3) + 70*pow(nx,2)*pow(tx,4)*pow(ty,4),
   56*pow(ny,2)*pow(tx,5)*pow(ty,3) + 140*nx*ny*pow(tx,4)*pow(ty,4) + 56*pow(nx,2)*pow(tx,3)*pow(ty,5),
   70*pow(ny,2)*pow(tx,4)*pow(ty,4) + 112*nx*ny*pow(tx,3)*pow(ty,5) + 28*pow(nx,2)*pow(tx,2)*pow(ty,6),
   56*pow(ny,2)*pow(tx,3)*pow(ty,5) + 56*nx*ny*pow(tx,2)*pow(ty,6) + 8*tx*pow(nx,2)*pow(ty,7),
   28*pow(ny,2)*pow(tx,2)*pow(ty,6) + 16*nx*ny*tx*pow(ty,7) + pow(nx,2)*pow(ty,8),8*tx*pow(ny,2)*pow(ty,7) + 2*nx*ny*pow(ty,8),pow(ny,2)*pow(ty,8),
   nx*pow(tx,9),9*nx*ty*pow(tx,8) + ny*pow(tx,9),9*ny*ty*pow(tx,8) + 36*nx*pow(tx,7)*pow(ty,2),36*ny*pow(tx,7)*pow(ty,2) + 84*nx*pow(tx,6)*pow(ty,3),
   84*ny*pow(tx,6)*pow(ty,3) + 126*nx*pow(tx,5)*pow(ty,4),126*ny*pow(tx,5)*pow(ty,4) + 126*nx*pow(tx,4)*pow(ty,5),
   126*ny*pow(tx,4)*pow(ty,5) + 84*nx*pow(tx,3)*pow(ty,6),84*ny*pow(tx,3)*pow(ty,6) + 36*nx*pow(tx,2)*pow(ty,7),
   36*ny*pow(tx,2)*pow(ty,7) + 9*nx*tx*pow(ty,8),9*ny*tx*pow(ty,8) + nx*pow(ty,9),ny*pow(ty,9),pow(tx,10),10*ty*pow(tx,9),45*pow(tx,8)*pow(ty,2),
   120*pow(tx,7)*pow(ty,3),210*pow(tx,6)*pow(ty,4),252*pow(tx,5)*pow(ty,5),210*pow(tx,4)*pow(ty,6),120*pow(tx,3)*pow(ty,7),45*pow(tx,2)*pow(ty,8),
   10*tx*pow(ty,9),pow(ty,10);


    degree = 11;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,11);


    tensor_project[degree].P <<pow(nx,11),11*ny*pow(nx,10),55*pow(nx,9)*pow(ny,2),165*pow(nx,8)*pow(ny,3),330*pow(nx,7)*pow(ny,4),462*pow(nx,6)*pow(ny,5),
   462*pow(nx,5)*pow(ny,6),330*pow(nx,4)*pow(ny,7),165*pow(nx,3)*pow(ny,8),55*pow(nx,2)*pow(ny,9),11*nx*pow(ny,10),pow(ny,11),tx*pow(nx,10),
   10*ny*tx*pow(nx,9) + ty*pow(nx,10),10*ny*ty*pow(nx,9) + 45*tx*pow(nx,8)*pow(ny,2),45*ty*pow(nx,8)*pow(ny,2) + 120*tx*pow(nx,7)*pow(ny,3),
   120*ty*pow(nx,7)*pow(ny,3) + 210*tx*pow(nx,6)*pow(ny,4),210*ty*pow(nx,6)*pow(ny,4) + 252*tx*pow(nx,5)*pow(ny,5),
   252*ty*pow(nx,5)*pow(ny,5) + 210*tx*pow(nx,4)*pow(ny,6),210*ty*pow(nx,4)*pow(ny,6) + 120*tx*pow(nx,3)*pow(ny,7),
   120*ty*pow(nx,3)*pow(ny,7) + 45*tx*pow(nx,2)*pow(ny,8),45*ty*pow(nx,2)*pow(ny,8) + 10*nx*tx*pow(ny,9),10*nx*ty*pow(ny,9) + tx*pow(ny,10),
   ty*pow(ny,10),pow(nx,9)*pow(tx,2),2*tx*ty*pow(nx,9) + 9*ny*pow(nx,8)*pow(tx,2),
   18*ny*tx*ty*pow(nx,8) + 36*pow(nx,7)*pow(ny,2)*pow(tx,2) + pow(nx,9)*pow(ty,2),
   72*tx*ty*pow(nx,7)*pow(ny,2) + 84*pow(nx,6)*pow(ny,3)*pow(tx,2) + 9*ny*pow(nx,8)*pow(ty,2),
   168*tx*ty*pow(nx,6)*pow(ny,3) + 126*pow(nx,5)*pow(ny,4)*pow(tx,2) + 36*pow(nx,7)*pow(ny,2)*pow(ty,2),
   252*tx*ty*pow(nx,5)*pow(ny,4) + 126*pow(nx,4)*pow(ny,5)*pow(tx,2) + 84*pow(nx,6)*pow(ny,3)*pow(ty,2),
   252*tx*ty*pow(nx,4)*pow(ny,5) + 84*pow(nx,3)*pow(ny,6)*pow(tx,2) + 126*pow(nx,5)*pow(ny,4)*pow(ty,2),
   168*tx*ty*pow(nx,3)*pow(ny,6) + 36*pow(nx,2)*pow(ny,7)*pow(tx,2) + 126*pow(nx,4)*pow(ny,5)*pow(ty,2),
   72*tx*ty*pow(nx,2)*pow(ny,7) + 9*nx*pow(ny,8)*pow(tx,2) + 84*pow(nx,3)*pow(ny,6)*pow(ty,2),
   18*nx*tx*ty*pow(ny,8) + pow(ny,9)*pow(tx,2) + 36*pow(nx,2)*pow(ny,7)*pow(ty,2),2*tx*ty*pow(ny,9) + 9*nx*pow(ny,8)*pow(ty,2),pow(ny,9)*pow(ty,2),
   pow(nx,8)*pow(tx,3),3*ty*pow(nx,8)*pow(tx,2) + 8*ny*pow(nx,7)*pow(tx,3),
   24*ny*ty*pow(nx,7)*pow(tx,2) + 28*pow(nx,6)*pow(ny,2)*pow(tx,3) + 3*tx*pow(nx,8)*pow(ty,2),
   84*ty*pow(nx,6)*pow(ny,2)*pow(tx,2) + 56*pow(nx,5)*pow(ny,3)*pow(tx,3) + 24*ny*tx*pow(nx,7)*pow(ty,2) + pow(nx,8)*pow(ty,3),
   168*ty*pow(nx,5)*pow(ny,3)*pow(tx,2) + 70*pow(nx,4)*pow(ny,4)*pow(tx,3) + 84*tx*pow(nx,6)*pow(ny,2)*pow(ty,2) + 8*ny*pow(nx,7)*pow(ty,3),
   210*ty*pow(nx,4)*pow(ny,4)*pow(tx,2) + 56*pow(nx,3)*pow(ny,5)*pow(tx,3) + 168*tx*pow(nx,5)*pow(ny,3)*pow(ty,2) + 28*pow(nx,6)*pow(ny,2)*pow(ty,3),
   168*ty*pow(nx,3)*pow(ny,5)*pow(tx,2) + 28*pow(nx,2)*pow(ny,6)*pow(tx,3) + 210*tx*pow(nx,4)*pow(ny,4)*pow(ty,2) + 56*pow(nx,5)*pow(ny,3)*pow(ty,3),
   84*ty*pow(nx,2)*pow(ny,6)*pow(tx,2) + 8*nx*pow(ny,7)*pow(tx,3) + 168*tx*pow(nx,3)*pow(ny,5)*pow(ty,2) + 70*pow(nx,4)*pow(ny,4)*pow(ty,3),
   24*nx*ty*pow(ny,7)*pow(tx,2) + pow(ny,8)*pow(tx,3) + 84*tx*pow(nx,2)*pow(ny,6)*pow(ty,2) + 56*pow(nx,3)*pow(ny,5)*pow(ty,3),
   3*ty*pow(ny,8)*pow(tx,2) + 24*nx*tx*pow(ny,7)*pow(ty,2) + 28*pow(nx,2)*pow(ny,6)*pow(ty,3),3*tx*pow(ny,8)*pow(ty,2) + 8*nx*pow(ny,7)*pow(ty,3),
   pow(ny,8)*pow(ty,3),pow(nx,7)*pow(tx,4),4*ty*pow(nx,7)*pow(tx,3) + 7*ny*pow(nx,6)*pow(tx,4),
   28*ny*ty*pow(nx,6)*pow(tx,3) + 21*pow(nx,5)*pow(ny,2)*pow(tx,4) + 6*pow(nx,7)*pow(tx,2)*pow(ty,2),
   84*ty*pow(nx,5)*pow(ny,2)*pow(tx,3) + 35*pow(nx,4)*pow(ny,3)*pow(tx,4) + 42*ny*pow(nx,6)*pow(tx,2)*pow(ty,2) + 4*tx*pow(nx,7)*pow(ty,3),
   140*ty*pow(nx,4)*pow(ny,3)*pow(tx,3) + 35*pow(nx,3)*pow(ny,4)*pow(tx,4) + 126*pow(nx,5)*pow(ny,2)*pow(tx,2)*pow(ty,2) + 
    28*ny*tx*pow(nx,6)*pow(ty,3) + pow(nx,7)*pow(ty,4),140*ty*pow(nx,3)*pow(ny,4)*pow(tx,3) + 21*pow(nx,2)*pow(ny,5)*pow(tx,4) + 
    210*pow(nx,4)*pow(ny,3)*pow(tx,2)*pow(ty,2) + 84*tx*pow(nx,5)*pow(ny,2)*pow(ty,3) + 7*ny*pow(nx,6)*pow(ty,4),
   84*ty*pow(nx,2)*pow(ny,5)*pow(tx,3) + 7*nx*pow(ny,6)*pow(tx,4) + 210*pow(nx,3)*pow(ny,4)*pow(tx,2)*pow(ty,2) + 
    140*tx*pow(nx,4)*pow(ny,3)*pow(ty,3) + 21*pow(nx,5)*pow(ny,2)*pow(ty,4),
   28*nx*ty*pow(ny,6)*pow(tx,3) + pow(ny,7)*pow(tx,4) + 126*pow(nx,2)*pow(ny,5)*pow(tx,2)*pow(ty,2) + 140*tx*pow(nx,3)*pow(ny,4)*pow(ty,3) + 
    35*pow(nx,4)*pow(ny,3)*pow(ty,4),4*ty*pow(ny,7)*pow(tx,3) + 42*nx*pow(ny,6)*pow(tx,2)*pow(ty,2) + 84*tx*pow(nx,2)*pow(ny,5)*pow(ty,3) + 
    35*pow(nx,3)*pow(ny,4)*pow(ty,4),6*pow(ny,7)*pow(tx,2)*pow(ty,2) + 28*nx*tx*pow(ny,6)*pow(ty,3) + 21*pow(nx,2)*pow(ny,5)*pow(ty,4),
   4*tx*pow(ny,7)*pow(ty,3) + 7*nx*pow(ny,6)*pow(ty,4),pow(ny,7)*pow(ty,4),pow(nx,6)*pow(tx,5),5*ty*pow(nx,6)*pow(tx,4) + 6*ny*pow(nx,5)*pow(tx,5),
   30*ny*ty*pow(nx,5)*pow(tx,4) + 15*pow(nx,4)*pow(ny,2)*pow(tx,5) + 10*pow(nx,6)*pow(tx,3)*pow(ty,2),
   75*ty*pow(nx,4)*pow(ny,2)*pow(tx,4) + 20*pow(nx,3)*pow(ny,3)*pow(tx,5) + 60*ny*pow(nx,5)*pow(tx,3)*pow(ty,2) + 10*pow(nx,6)*pow(tx,2)*pow(ty,3),
   100*ty*pow(nx,3)*pow(ny,3)*pow(tx,4) + 15*pow(nx,2)*pow(ny,4)*pow(tx,5) + 150*pow(nx,4)*pow(ny,2)*pow(tx,3)*pow(ty,2) + 
    60*ny*pow(nx,5)*pow(tx,2)*pow(ty,3) + 5*tx*pow(nx,6)*pow(ty,4),
   75*ty*pow(nx,2)*pow(ny,4)*pow(tx,4) + 6*nx*pow(ny,5)*pow(tx,5) + 200*pow(nx,3)*pow(ny,3)*pow(tx,3)*pow(ty,2) + 
    150*pow(nx,4)*pow(ny,2)*pow(tx,2)*pow(ty,3) + 30*ny*tx*pow(nx,5)*pow(ty,4) + pow(nx,6)*pow(ty,5),
   30*nx*ty*pow(ny,5)*pow(tx,4) + pow(ny,6)*pow(tx,5) + 150*pow(nx,2)*pow(ny,4)*pow(tx,3)*pow(ty,2) + 200*pow(nx,3)*pow(ny,3)*pow(tx,2)*pow(ty,3) + 
    75*tx*pow(nx,4)*pow(ny,2)*pow(ty,4) + 6*ny*pow(nx,5)*pow(ty,5),
   5*ty*pow(ny,6)*pow(tx,4) + 60*nx*pow(ny,5)*pow(tx,3)*pow(ty,2) + 150*pow(nx,2)*pow(ny,4)*pow(tx,2)*pow(ty,3) + 
    100*tx*pow(nx,3)*pow(ny,3)*pow(ty,4) + 15*pow(nx,4)*pow(ny,2)*pow(ty,5),
   10*pow(ny,6)*pow(tx,3)*pow(ty,2) + 60*nx*pow(ny,5)*pow(tx,2)*pow(ty,3) + 75*tx*pow(nx,2)*pow(ny,4)*pow(ty,4) + 20*pow(nx,3)*pow(ny,3)*pow(ty,5),
   10*pow(ny,6)*pow(tx,2)*pow(ty,3) + 30*nx*tx*pow(ny,5)*pow(ty,4) + 15*pow(nx,2)*pow(ny,4)*pow(ty,5),
   5*tx*pow(ny,6)*pow(ty,4) + 6*nx*pow(ny,5)*pow(ty,5),pow(ny,6)*pow(ty,5),pow(nx,5)*pow(tx,6),6*ty*pow(nx,5)*pow(tx,5) + 5*ny*pow(nx,4)*pow(tx,6),
   30*ny*ty*pow(nx,4)*pow(tx,5) + 10*pow(nx,3)*pow(ny,2)*pow(tx,6) + 15*pow(nx,5)*pow(tx,4)*pow(ty,2),
   60*ty*pow(nx,3)*pow(ny,2)*pow(tx,5) + 10*pow(nx,2)*pow(ny,3)*pow(tx,6) + 75*ny*pow(nx,4)*pow(tx,4)*pow(ty,2) + 20*pow(nx,5)*pow(tx,3)*pow(ty,3),
   60*ty*pow(nx,2)*pow(ny,3)*pow(tx,5) + 5*nx*pow(ny,4)*pow(tx,6) + 150*pow(nx,3)*pow(ny,2)*pow(tx,4)*pow(ty,2) + 
    100*ny*pow(nx,4)*pow(tx,3)*pow(ty,3) + 15*pow(nx,5)*pow(tx,2)*pow(ty,4),
   30*nx*ty*pow(ny,4)*pow(tx,5) + pow(ny,5)*pow(tx,6) + 150*pow(nx,2)*pow(ny,3)*pow(tx,4)*pow(ty,2) + 200*pow(nx,3)*pow(ny,2)*pow(tx,3)*pow(ty,3) + 
    75*ny*pow(nx,4)*pow(tx,2)*pow(ty,4) + 6*tx*pow(nx,5)*pow(ty,5),
   6*ty*pow(ny,5)*pow(tx,5) + 75*nx*pow(ny,4)*pow(tx,4)*pow(ty,2) + 200*pow(nx,2)*pow(ny,3)*pow(tx,3)*pow(ty,3) + 
    150*pow(nx,3)*pow(ny,2)*pow(tx,2)*pow(ty,4) + 30*ny*tx*pow(nx,4)*pow(ty,5) + pow(nx,5)*pow(ty,6),
   15*pow(ny,5)*pow(tx,4)*pow(ty,2) + 100*nx*pow(ny,4)*pow(tx,3)*pow(ty,3) + 150*pow(nx,2)*pow(ny,3)*pow(tx,2)*pow(ty,4) + 
    60*tx*pow(nx,3)*pow(ny,2)*pow(ty,5) + 5*ny*pow(nx,4)*pow(ty,6),
   20*pow(ny,5)*pow(tx,3)*pow(ty,3) + 75*nx*pow(ny,4)*pow(tx,2)*pow(ty,4) + 60*tx*pow(nx,2)*pow(ny,3)*pow(ty,5) + 10*pow(nx,3)*pow(ny,2)*pow(ty,6),
   15*pow(ny,5)*pow(tx,2)*pow(ty,4) + 30*nx*tx*pow(ny,4)*pow(ty,5) + 10*pow(nx,2)*pow(ny,3)*pow(ty,6),
   6*tx*pow(ny,5)*pow(ty,5) + 5*nx*pow(ny,4)*pow(ty,6),pow(ny,5)*pow(ty,6),pow(nx,4)*pow(tx,7),7*ty*pow(nx,4)*pow(tx,6) + 4*ny*pow(nx,3)*pow(tx,7),
   28*ny*ty*pow(nx,3)*pow(tx,6) + 6*pow(nx,2)*pow(ny,2)*pow(tx,7) + 21*pow(nx,4)*pow(tx,5)*pow(ty,2),
   42*ty*pow(nx,2)*pow(ny,2)*pow(tx,6) + 4*nx*pow(ny,3)*pow(tx,7) + 84*ny*pow(nx,3)*pow(tx,5)*pow(ty,2) + 35*pow(nx,4)*pow(tx,4)*pow(ty,3),
   28*nx*ty*pow(ny,3)*pow(tx,6) + pow(ny,4)*pow(tx,7) + 126*pow(nx,2)*pow(ny,2)*pow(tx,5)*pow(ty,2) + 140*ny*pow(nx,3)*pow(tx,4)*pow(ty,3) + 
    35*pow(nx,4)*pow(tx,3)*pow(ty,4),7*ty*pow(ny,4)*pow(tx,6) + 84*nx*pow(ny,3)*pow(tx,5)*pow(ty,2) + 210*pow(nx,2)*pow(ny,2)*pow(tx,4)*pow(ty,3) + 
    140*ny*pow(nx,3)*pow(tx,3)*pow(ty,4) + 21*pow(nx,4)*pow(tx,2)*pow(ty,5),
   21*pow(ny,4)*pow(tx,5)*pow(ty,2) + 140*nx*pow(ny,3)*pow(tx,4)*pow(ty,3) + 210*pow(nx,2)*pow(ny,2)*pow(tx,3)*pow(ty,4) + 
    84*ny*pow(nx,3)*pow(tx,2)*pow(ty,5) + 7*tx*pow(nx,4)*pow(ty,6),
   35*pow(ny,4)*pow(tx,4)*pow(ty,3) + 140*nx*pow(ny,3)*pow(tx,3)*pow(ty,4) + 126*pow(nx,2)*pow(ny,2)*pow(tx,2)*pow(ty,5) + 
    28*ny*tx*pow(nx,3)*pow(ty,6) + pow(nx,4)*pow(ty,7),35*pow(ny,4)*pow(tx,3)*pow(ty,4) + 84*nx*pow(ny,3)*pow(tx,2)*pow(ty,5) + 
    42*tx*pow(nx,2)*pow(ny,2)*pow(ty,6) + 4*ny*pow(nx,3)*pow(ty,7),
   21*pow(ny,4)*pow(tx,2)*pow(ty,5) + 28*nx*tx*pow(ny,3)*pow(ty,6) + 6*pow(nx,2)*pow(ny,2)*pow(ty,7),
   7*tx*pow(ny,4)*pow(ty,6) + 4*nx*pow(ny,3)*pow(ty,7),pow(ny,4)*pow(ty,7),pow(nx,3)*pow(tx,8),8*ty*pow(nx,3)*pow(tx,7) + 3*ny*pow(nx,2)*pow(tx,8),
   24*ny*ty*pow(nx,2)*pow(tx,7) + 3*nx*pow(ny,2)*pow(tx,8) + 28*pow(nx,3)*pow(tx,6)*pow(ty,2),
   24*nx*ty*pow(ny,2)*pow(tx,7) + pow(ny,3)*pow(tx,8) + 84*ny*pow(nx,2)*pow(tx,6)*pow(ty,2) + 56*pow(nx,3)*pow(tx,5)*pow(ty,3),
   8*ty*pow(ny,3)*pow(tx,7) + 84*nx*pow(ny,2)*pow(tx,6)*pow(ty,2) + 168*ny*pow(nx,2)*pow(tx,5)*pow(ty,3) + 70*pow(nx,3)*pow(tx,4)*pow(ty,4),
   28*pow(ny,3)*pow(tx,6)*pow(ty,2) + 168*nx*pow(ny,2)*pow(tx,5)*pow(ty,3) + 210*ny*pow(nx,2)*pow(tx,4)*pow(ty,4) + 56*pow(nx,3)*pow(tx,3)*pow(ty,5),
   56*pow(ny,3)*pow(tx,5)*pow(ty,3) + 210*nx*pow(ny,2)*pow(tx,4)*pow(ty,4) + 168*ny*pow(nx,2)*pow(tx,3)*pow(ty,5) + 28*pow(nx,3)*pow(tx,2)*pow(ty,6),
   70*pow(ny,3)*pow(tx,4)*pow(ty,4) + 168*nx*pow(ny,2)*pow(tx,3)*pow(ty,5) + 84*ny*pow(nx,2)*pow(tx,2)*pow(ty,6) + 8*tx*pow(nx,3)*pow(ty,7),
   56*pow(ny,3)*pow(tx,3)*pow(ty,5) + 84*nx*pow(ny,2)*pow(tx,2)*pow(ty,6) + 24*ny*tx*pow(nx,2)*pow(ty,7) + pow(nx,3)*pow(ty,8),
   28*pow(ny,3)*pow(tx,2)*pow(ty,6) + 24*nx*tx*pow(ny,2)*pow(ty,7) + 3*ny*pow(nx,2)*pow(ty,8),8*tx*pow(ny,3)*pow(ty,7) + 3*nx*pow(ny,2)*pow(ty,8),
   pow(ny,3)*pow(ty,8),pow(nx,2)*pow(tx,9),9*ty*pow(nx,2)*pow(tx,8) + 2*nx*ny*pow(tx,9),
   18*nx*ny*ty*pow(tx,8) + pow(ny,2)*pow(tx,9) + 36*pow(nx,2)*pow(tx,7)*pow(ty,2),
   9*ty*pow(ny,2)*pow(tx,8) + 72*nx*ny*pow(tx,7)*pow(ty,2) + 84*pow(nx,2)*pow(tx,6)*pow(ty,3),
   36*pow(ny,2)*pow(tx,7)*pow(ty,2) + 168*nx*ny*pow(tx,6)*pow(ty,3) + 126*pow(nx,2)*pow(tx,5)*pow(ty,4),
   84*pow(ny,2)*pow(tx,6)*pow(ty,3) + 252*nx*ny*pow(tx,5)*pow(ty,4) + 126*pow(nx,2)*pow(tx,4)*pow(ty,5),
   126*pow(ny,2)*pow(tx,5)*pow(ty,4) + 252*nx*ny*pow(tx,4)*pow(ty,5) + 84*pow(nx,2)*pow(tx,3)*pow(ty,6),
   126*pow(ny,2)*pow(tx,4)*pow(ty,5) + 168*nx*ny*pow(tx,3)*pow(ty,6) + 36*pow(nx,2)*pow(tx,2)*pow(ty,7),
   84*pow(ny,2)*pow(tx,3)*pow(ty,6) + 72*nx*ny*pow(tx,2)*pow(ty,7) + 9*tx*pow(nx,2)*pow(ty,8),
   36*pow(ny,2)*pow(tx,2)*pow(ty,7) + 18*nx*ny*tx*pow(ty,8) + pow(nx,2)*pow(ty,9),9*tx*pow(ny,2)*pow(ty,8) + 2*nx*ny*pow(ty,9),pow(ny,2)*pow(ty,9),
   nx*pow(tx,10),10*nx*ty*pow(tx,9) + ny*pow(tx,10),10*ny*ty*pow(tx,9) + 45*nx*pow(tx,8)*pow(ty,2),
   45*ny*pow(tx,8)*pow(ty,2) + 120*nx*pow(tx,7)*pow(ty,3),120*ny*pow(tx,7)*pow(ty,3) + 210*nx*pow(tx,6)*pow(ty,4),
   210*ny*pow(tx,6)*pow(ty,4) + 252*nx*pow(tx,5)*pow(ty,5),252*ny*pow(tx,5)*pow(ty,5) + 210*nx*pow(tx,4)*pow(ty,6),
   210*ny*pow(tx,4)*pow(ty,6) + 120*nx*pow(tx,3)*pow(ty,7),120*ny*pow(tx,3)*pow(ty,7) + 45*nx*pow(tx,2)*pow(ty,8),
   45*ny*pow(tx,2)*pow(ty,8) + 10*nx*tx*pow(ty,9),10*ny*tx*pow(ty,9) + nx*pow(ty,10),ny*pow(ty,10),pow(tx,11),11*ty*pow(tx,10),
   55*pow(tx,9)*pow(ty,2),165*pow(tx,8)*pow(ty,3),330*pow(tx,7)*pow(ty,4),462*pow(tx,6)*pow(ty,5),462*pow(tx,5)*pow(ty,6),330*pow(tx,4)*pow(ty,7),
   165*pow(tx,3)*pow(ty,8),55*pow(tx,2)*pow(ty,9),11*tx*pow(ty,10),pow(ty,11);    
 }


  template<>
  void 
  Base_TensorInfo<3>
  ::reinit_local(const double nx,const double ny,const double nz,
                 const double tx,const double ty,const double tz,
                 const double rx,const double ry,const double rz,
                std::vector<projector_data> &tensor_project)
  {
    // first we allocate the required memory
    tensor_project.resize(max_tensorial_degree + 1);

    // now we allocate memory for every individual projector
    allocate_tensor_memory(tensor_project);

    AssertDimension(tensor_project.size(),max_tensorial_degree + 1);

    // we first check whether P has been properly allocated or not
    for (unsigned int i = 0 ; i < max_tensorial_degree ; i++)
    {
      AssertDimension(tensor_project[i].P.rows(),components(i));
      AssertDimension(tensor_project[i].P.cols(),components(i));
    }

    unsigned int degree = 0;

    tensor_project[0].P << 1;

    degree = 1;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P << nx,ny,nz,tx,ty,tz,rx,ry,rz;  
    

    degree = 2;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P << pow(nx,2) - pow(nz,2),2*nx*ny,2*nx*nz,pow(ny,2) - pow(nz,2),2*ny*nz,nx*tx - nz*tz,
   ny*tx + nx*ty,nz*tx + nx*tz,ny*ty - nz*tz,nz*ty + ny*tz,nx*rx - nz*rz,ny*rx + nx*ry,
   nz*rx + nx*rz,ny*ry - nz*rz,nz*ry + ny*rz,pow(tx,2) - pow(tz,2),2*tx*ty,2*tx*tz,
   pow(ty,2) - pow(tz,2),2*ty*tz,rx*tx - rz*tz,ry*tx + rx*ty,rz*tx + rx*tz,ry*ty - rz*tz,
   rz*ty + ry*tz; // end of the third row


    degree = 3;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P << pow(nx,3) - 3*nx*pow(nz,2),3*ny*(pow(nx,2) - pow(nz,2)),3*nz*pow(nx,2) - pow(nz,3),
   3*nx*(pow(ny,2) - pow(nz,2)),6*nx*ny*nz,pow(ny,3) - 3*ny*pow(nz,2),
   3*nz*pow(ny,2) - pow(nz,3),-2*nx*nz*tz + tx*pow(nx,2) - tx*pow(nz,2),
   2*nx*ny*tx - nz*(nz*ty + 2*ny*tz) + ty*pow(nx,2),
   2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2),
   2*nx*ny*ty - nz*(nz*tx + 2*nx*tz) + tx*pow(ny,2),2*(ny*nz*tx + nx*nz*ty + nx*ny*tz),
   -2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2),2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2),
   -2*nx*nz*rz + rx*pow(nx,2) - rx*pow(nz,2),
   2*nx*ny*rx - nz*(nz*ry + 2*ny*rz) + ry*pow(nx,2),
   2*nx*nz*rx + rz*pow(nx,2) - rz*pow(nz,2),
   2*nx*ny*ry - nz*(nz*rx + 2*nx*rz) + rx*pow(ny,2),2*(ny*nz*rx + nx*nz*ry + nx*ny*rz),
   -2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2),2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2),
   -2*nz*tx*tz + nx*(pow(tx,2) - pow(tz,2)),
   2*nx*tx*ty - 2*nz*ty*tz + ny*pow(tx,2) - ny*pow(tz,2),
   2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2)),
   2*ny*tx*ty - 2*nz*tx*tz + nx*(pow(ty,2) - pow(tz,2)),
   2*(nz*tx*ty + ny*tx*tz + nx*ty*tz),-2*nz*ty*tz + ny*(pow(ty,2) - pow(tz,2)),
   2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2)),-(nz*(rz*tx + rx*tz)) + nx*(rx*tx - rz*tz),
   ny*rx*tx + nx*ry*tx + nx*rx*ty - nz*rz*ty - nz*ry*tz - ny*rz*tz,
   nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz,
   ny*ry*tx - nz*rz*tx + ny*rx*ty + nx*ry*ty - nz*rx*tz - nx*rz*tz,
   nz*ry*tx + ny*rz*tx + nz*rx*ty + nx*rz*ty + ny*rx*tz + nx*ry*tz,
   -(nz*(rz*ty + ry*tz)) + ny*(ry*ty - rz*tz),nz*ry*ty + ny*rz*ty + ny*ry*tz - nz*rz*tz,
   pow(tx,3) - 3*tx*pow(tz,2),3*ty*(pow(tx,2) - pow(tz,2)),3*tz*pow(tx,2) - pow(tz,3),
   3*tx*(pow(ty,2) - pow(tz,2)),6*tx*ty*tz,pow(ty,3) - 3*ty*pow(tz,2),
   3*tz*pow(ty,2) - pow(tz,3),-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)),
   2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2),
   2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2)),
   2*ry*tx*ty - 2*rz*tx*tz + rx*(pow(ty,2) - pow(tz,2)),
   2*(rz*tx*ty + ry*tx*tz + rx*ty*tz),-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)),
   2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2)); // fourth row



    degree = 4;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
      tensor_project[degree].P << pow(nx,4) - 6*pow(nx,2)*pow(nz,2) + pow(nz,4),4*nx*ny*(pow(nx,2) - 3*pow(nz,2)),4*nx*nz*(pow(nx,2) - pow(nz,2)),
   2*(3*pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 3*pow(ny,2)*pow(nz,2) + pow(nz,4)),-4*ny*nz*(-3*pow(nx,2) + pow(nz,2)),4*nx*ny*(pow(ny,2) - 3*pow(nz,2)),
   -4*nx*nz*(-3*pow(ny,2) + pow(nz,2)),pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4),4*ny*nz*(pow(ny,2) - pow(nz,2)),
   -3*nz*tz*pow(nx,2) + tx*pow(nx,3) - 3*nx*tx*pow(nz,2) + tz*pow(nz,3),
   -3*nx*nz*(nz*ty + 2*ny*tz) + 3*ny*tx*pow(nx,2) + ty*pow(nx,3) - 3*ny*tx*pow(nz,2),
   3*nz*tx*pow(nx,2) + tz*pow(nx,3) - 3*nx*tz*pow(nz,2) - tx*pow(nz,3),
   3*(ny*ty - nz*tz)*pow(nx,2) + 3*nx*tx*(pow(ny,2) - pow(nz,2)) + nz*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)),
   6*nx*ny*nz*tx + 3*(nz*ty + ny*tz)*pow(nx,2) - (nz*ty + 3*ny*tz)*pow(nz,2),
   -3*ny*nz*(nz*tx + 2*nx*tz) + 3*nx*ty*pow(ny,2) + tx*pow(ny,3) - 3*nx*ty*pow(nz,2),
   6*nx*ny*nz*ty + 3*(nz*tx + nx*tz)*pow(ny,2) - (nz*tx + 3*nx*tz)*pow(nz,2),-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3),
   3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3),-3*nz*rz*pow(nx,2) + rx*pow(nx,3) - 3*nx*rx*pow(nz,2) + rz*pow(nz,3),
   -3*nx*nz*(nz*ry + 2*ny*rz) + 3*ny*rx*pow(nx,2) + ry*pow(nx,3) - 3*ny*rx*pow(nz,2),
   3*nz*rx*pow(nx,2) + rz*pow(nx,3) - 3*nx*rz*pow(nz,2) - rx*pow(nz,3),
   3*(ny*ry - nz*rz)*pow(nx,2) + 3*nx*rx*(pow(ny,2) - pow(nz,2)) + nz*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2)),
   6*nx*ny*nz*rx + 3*(nz*ry + ny*rz)*pow(nx,2) - (nz*ry + 3*ny*rz)*pow(nz,2),
   -3*ny*nz*(nz*rx + 2*nx*rz) + 3*nx*ry*pow(ny,2) + rx*pow(ny,3) - 3*nx*ry*pow(nz,2),
   6*nx*ny*nz*ry + 3*(nz*rx + nx*rz)*pow(ny,2) - (nz*rx + 3*nx*rz)*pow(nz,2),-3*nz*rz*pow(ny,2) + ry*pow(ny,3) - 3*ny*ry*pow(nz,2) + rz*pow(nz,3),
   3*nz*ry*pow(ny,2) + rz*pow(ny,3) - 3*ny*rz*pow(nz,2) - ry*pow(nz,3),
   -4*nx*nz*tx*tz + pow(nx,2)*(pow(tx,2) - pow(tz,2)) + pow(nz,2)*(-pow(tx,2) + pow(tz,2)),
   2*(-(nz*tx*(nz*ty + 2*ny*tz)) + tx*ty*pow(nx,2) + nx*(-2*nz*ty*tz + ny*(pow(tx,2) - pow(tz,2)))),
   2*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))),
   -4*nx*nz*tx*tz + 4*ny*ty*(nx*tx - nz*tz) - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
    pow(nx,2)*(pow(ty,2) - pow(tz,2)),2*(ty*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + ny*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2)))),
   2*(-(nz*ty*(nz*tx + 2*nx*tz)) + tx*ty*pow(ny,2) + ny*(-2*nz*tx*tz + nx*(pow(ty,2) - pow(tz,2)))),
   2*(2*ny*ty*(nz*tx + nx*tz) + tx*tz*pow(ny,2) - nz*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   -4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2)),
   2*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2))),
   -2*nx*nz*(rz*tx + rx*tz) + (rx*tx - rz*tz)*pow(nx,2) + (-(rx*tx) + rz*tz)*pow(nz,2),
   -(nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz)) - 2*nx*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz) + (ry*tx + rx*ty)*pow(nx,2),
   2*nx*nz*(rx*tx - rz*tz) + (rz*tx + rx*tz)*pow(nx,2) - (rz*tx + rx*tz)*pow(nz,2),
   -2*nx*nz*(rz*tx + rx*tz) + 2*ny*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz)) + (ry*ty - rz*tz)*pow(nx,2) + (rx*tx - rz*tz)*pow(ny,2) - 
    (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2),2*nx*nz*(ry*tx + rx*ty) + 2*ny*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz) + (rz*ty + ry*tz)*pow(nx,2) - 
    (rz*ty + ry*tz)*pow(nz,2),-(nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz)) - 2*ny*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz) + 
    (ry*tx + rx*ty)*pow(ny,2),2*ny*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz) - nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz) + 
    (rz*tx + rx*tz)*pow(ny,2),-2*ny*nz*(rz*ty + ry*tz) + (ry*ty - rz*tz)*pow(ny,2) + (-(ry*ty) + rz*tz)*pow(nz,2),
   2*ny*nz*(ry*ty - rz*tz) + (rz*ty + ry*tz)*pow(ny,2) - (rz*ty + ry*tz)*pow(nz,2),
   nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2)),
   ny*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))),
   3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3),
   3*ny*ty*(pow(tx,2) - pow(tz,2)) + 3*nx*tx*(pow(ty,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)),
   3*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*nx*tx*ty + 3*ny*pow(tx,2) - ny*pow(tz,2)),
   ty*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 3*ny*tx*(pow(ty,2) - pow(tz,2)),
   3*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ny*tx*ty + 3*nx*pow(ty,2) - nx*pow(tz,2)),
   nz*tz*(-3*pow(ty,2) + pow(tz,2)) + ny*(pow(ty,3) - 3*ty*pow(tz,2)),3*ny*tz*pow(ty,2) + nz*pow(ty,3) - 3*nz*ty*pow(tz,2) - ny*pow(tz,3),
   nx*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + nz*(-2*rx*tx*tz + rz*(-pow(tx,2) + pow(tz,2))),
   -2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + 
    nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)),
   nz*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + nx*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))),
   -(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + nx*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    ny*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)),
   2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)),-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
    ny*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2)),
   2*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)),ny*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + nz*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))),
   nz*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + ny*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))),pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4),
   4*tx*ty*(pow(tx,2) - 3*pow(tz,2)),4*tx*tz*(pow(tx,2) - pow(tz,2)),2*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   -4*ty*tz*(-3*pow(tx,2) + pow(tz,2)),4*tx*ty*(pow(ty,2) - 3*pow(tz,2)),-4*tx*tz*(-3*pow(ty,2) + pow(tz,2)),
   pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4),4*ty*tz*(pow(ty,2) - pow(tz,2)),rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)),
   ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3),
   3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)),
   3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2)),
   ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2)),
   3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)),
   rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2)),3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3);

   

    degree = 5;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
      tensor_project[degree].P << pow(nx,5) - 10*pow(nx,3)*pow(nz,2) + 5*nx*pow(nz,4),5*ny*(pow(nx,4) - 6*pow(nx,2)*pow(nz,2) + pow(nz,4)),
   5*nz*pow(nx,4) - 10*pow(nx,2)*pow(nz,3) + pow(nz,5),10*nx*(pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 3*pow(ny,2)*pow(nz,2) + pow(nz,4)),
   20*nx*ny*nz*(pow(nx,2) - pow(nz,2)),10*ny*(pow(nx,2)*(pow(ny,2) - 3*pow(nz,2)) - pow(ny,2)*pow(nz,2) + pow(nz,4)),
   2*(5*pow(nx,2)*(3*nz*pow(ny,2) - pow(nz,3)) - 5*pow(ny,2)*pow(nz,3) + pow(nz,5)),5*nx*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)),
   20*nx*ny*nz*(pow(ny,2) - pow(nz,2)),pow(ny,5) - 10*pow(ny,3)*pow(nz,2) + 5*ny*pow(nz,4),5*nz*pow(ny,4) - 10*pow(ny,2)*pow(nz,3) + pow(nz,5),
   -4*nz*tz*pow(nx,3) + tx*pow(nx,4) - 6*tx*pow(nx,2)*pow(nz,2) + 4*nx*tz*pow(nz,3) + tx*pow(nz,4),
   -6*nz*(nz*ty + 2*ny*tz)*pow(nx,2) + 4*ny*tx*pow(nx,3) + ty*pow(nx,4) - 12*nx*ny*tx*pow(nz,2) + (nz*ty + 4*ny*tz)*pow(nz,3),
   4*nz*tx*pow(nx,3) + tz*pow(nx,4) - 6*tz*pow(nx,2)*pow(nz,2) - 4*nx*tx*pow(nz,3) + tz*pow(nz,4),
   2*(2*(ny*ty - nz*tz)*pow(nx,3) + 3*tx*pow(nx,2)*(pow(ny,2) - pow(nz,2)) + tx*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
      2*nx*nz*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2))),
   4*(3*ny*nz*tx*pow(nx,2) + (nz*ty + ny*tz)*pow(nx,3) - nx*(nz*ty + 3*ny*tz)*pow(nz,2) - ny*tx*pow(nz,3)),
   2*(2*nx*ny*tx*(pow(ny,2) - 3*pow(nz,2)) + 3*pow(nx,2)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
      nz*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3))),
   2*(3*pow(nx,2)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 
      nx*(6*nz*tx*pow(ny,2) - 2*tx*pow(nz,3))),-6*nz*(nz*tx + 2*nx*tz)*pow(ny,2) + 4*nx*ty*pow(ny,3) + tx*pow(ny,4) - 12*nx*ny*ty*pow(nz,2) + 
    (nz*tx + 4*nx*tz)*pow(nz,3),4*(3*nx*nz*ty*pow(ny,2) + (nz*tx + nx*tz)*pow(ny,3) - ny*(nz*tx + 3*nx*tz)*pow(nz,2) - nx*ty*pow(nz,3)),
   -4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4),
   4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4),
   -4*nz*rz*pow(nx,3) + rx*pow(nx,4) - 6*rx*pow(nx,2)*pow(nz,2) + 4*nx*rz*pow(nz,3) + rx*pow(nz,4),
   -6*nz*(nz*ry + 2*ny*rz)*pow(nx,2) + 4*ny*rx*pow(nx,3) + ry*pow(nx,4) - 12*nx*ny*rx*pow(nz,2) + (nz*ry + 4*ny*rz)*pow(nz,3),
   4*nz*rx*pow(nx,3) + rz*pow(nx,4) - 6*rz*pow(nx,2)*pow(nz,2) - 4*nx*rx*pow(nz,3) + rz*pow(nz,4),
   2*(2*(ny*ry - nz*rz)*pow(nx,3) + 3*rx*pow(nx,2)*(pow(ny,2) - pow(nz,2)) + rx*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
      2*nx*nz*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2))),
   4*(3*ny*nz*rx*pow(nx,2) + (nz*ry + ny*rz)*pow(nx,3) - nx*(nz*ry + 3*ny*rz)*pow(nz,2) - ny*rx*pow(nz,3)),
   2*(2*nx*ny*rx*(pow(ny,2) - 3*pow(nz,2)) + 3*pow(nx,2)*(-2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2)) + 
      nz*(-3*nz*ry*pow(ny,2) - 2*rz*pow(ny,3) + 4*ny*rz*pow(nz,2) + ry*pow(nz,3))),
   2*(3*pow(nx,2)*(2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2)) + pow(nz,2)*(-2*ny*nz*ry - 3*rz*pow(ny,2) + rz*pow(nz,2)) + 
      nx*(6*nz*rx*pow(ny,2) - 2*rx*pow(nz,3))),-6*nz*(nz*rx + 2*nx*rz)*pow(ny,2) + 4*nx*ry*pow(ny,3) + rx*pow(ny,4) - 12*nx*ny*ry*pow(nz,2) + 
    (nz*rx + 4*nx*rz)*pow(nz,3),4*(3*nx*nz*ry*pow(ny,2) + (nz*rx + nx*rz)*pow(ny,3) - ny*(nz*rx + 3*nx*rz)*pow(nz,2) - nx*ry*pow(nz,3)),
   -4*nz*rz*pow(ny,3) + ry*pow(ny,4) - 6*ry*pow(ny,2)*pow(nz,2) + 4*ny*rz*pow(nz,3) + ry*pow(nz,4),
   4*nz*ry*pow(ny,3) + rz*pow(ny,4) - 6*rz*pow(ny,2)*pow(nz,2) - 4*ny*ry*pow(nz,3) + rz*pow(nz,4),
   -6*nz*tx*tz*pow(nx,2) + 2*tx*tz*pow(nz,3) + pow(nx,3)*(pow(tx,2) - pow(tz,2)) + 3*nx*pow(nz,2)*(-pow(tx,2) + pow(tz,2)),
   -6*nx*nz*tx*(nz*ty + 2*ny*tz) + 2*tx*ty*pow(nx,3) + pow(nz,2)*(2*nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 
    3*pow(nx,2)*(2*nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2))),2*tx*tz*pow(nx,3) - 6*nx*tx*tz*pow(nz,2) + 3*nz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 
    pow(nz,3)*(-pow(tx,2) + pow(tz,2)),6*tx*(ny*ty - nz*tz)*pow(nx,2) + 2*nz*tx*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 
    3*nx*(-4*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2))) + pow(nx,3)*(pow(ty,2) - pow(tz,2)),
   6*tx*(nz*ty + ny*tz)*pow(nx,2) + 2*ty*tz*pow(nx,3) - 2*tx*(nz*ty + 3*ny*tz)*pow(nz,2) - 6*nx*nz*(nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2))),
   6*ty*(nx*tx - nz*tz)*pow(ny,2) + 2*nz*ty*(-3*nx*nz*tx - 3*tz*pow(nx,2) + 2*tz*pow(nz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2)) - 
    3*ny*(4*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   6*ny*ty*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + 3*pow(ny,2)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
    nz*(6*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 3*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   -6*ny*nz*ty*(nz*tx + 2*nx*tz) + 2*tx*ty*pow(ny,3) + pow(nz,2)*(2*nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 
    3*pow(ny,2)*(2*nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2))),6*ty*(nz*tx + nx*tz)*pow(ny,2) + 2*tx*tz*pow(ny,3) - 2*ty*(nz*tx + 3*nx*tz)*pow(nz,2) - 
    6*ny*nz*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2))),-6*nz*ty*tz*pow(ny,2) + 2*ty*tz*pow(nz,3) + pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 
    3*ny*pow(nz,2)*(-pow(ty,2) + pow(tz,2)),2*ty*tz*pow(ny,3) - 6*ny*ty*tz*pow(nz,2) + 3*nz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + 
    pow(nz,3)*(-pow(ty,2) + pow(tz,2)),-3*nz*(rz*tx + rx*tz)*pow(nx,2) + (rx*tx - rz*tz)*pow(nx,3) + 3*nx*(-(rx*tx) + rz*tz)*pow(nz,2) + 
    (rz*tx + rx*tz)*pow(nz,3),-3*nx*nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz) - 
    3*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,2) + (ry*tx + rx*ty)*pow(nx,3) + 
    (-3*ny*rx*tx + nz*rz*ty + nz*ry*tz + 3*ny*rz*tz)*pow(nz,2),
   3*nz*(rx*tx - rz*tz)*pow(nx,2) + (rz*tx + rx*tz)*pow(nx,3) - 3*nx*(rz*tx + rx*tz)*pow(nz,2) + (-(rx*tx) + rz*tz)*pow(nz,3),
   3*(ny*(ry*tx + rx*ty) - nz*(rz*tx + rx*tz))*pow(nx,2) + (ry*ty - rz*tz)*pow(nx,3) + 
    nz*(-3*ny*nz*(ry*tx + rx*ty) - 3*(rz*tx + rx*tz)*pow(ny,2) + 2*(rz*tx + rx*tz)*pow(nz,2)) + 
    3*nx*(-2*ny*nz*(rz*ty + ry*tz) + (rx*tx - rz*tz)*pow(ny,2) - (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)),
   -3*nx*nz*(-2*ny*rx*tx + nz*rz*ty + nz*ry*tz + 2*ny*rz*tz) + 3*(nz*ry*tx + ny*rz*tx + nz*rx*ty + ny*rx*tz)*pow(nx,2) + (rz*ty + ry*tz)*pow(nx,3) - 
    (nz*ry*tx + 3*ny*rz*tx + nz*rx*ty + 3*ny*rx*tz)*pow(nz,2),
   3*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz))*pow(ny,2) + (rx*tx - rz*tz)*pow(ny,3) + 
    nz*(-3*nx*nz*(ry*tx + rx*ty) - 3*(rz*ty + ry*tz)*pow(nx,2) + 2*(rz*ty + ry*tz)*pow(nz,2)) - 
    3*ny*(2*nx*nz*(rz*tx + rx*tz) + (-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)),
   3*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz)*pow(ny,2) + 
    3*ny*(2*nx*nz*(ry*tx + rx*ty) + (rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    nz*(3*nx*nz*(rz*tx + rx*tz) + 3*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)),
   -3*ny*nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz) - 3*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz)*pow(ny,2) + (ry*tx + rx*ty)*pow(ny,3) + 
    (nz*rz*tx - 3*nx*ry*ty + nz*rx*tz + 3*nx*rz*tz)*pow(nz,2),
   -3*ny*nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz) + 3*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz)*pow(ny,2) + (rz*tx + rx*tz)*pow(ny,3) - 
    (nz*ry*tx + nz*rx*ty + 3*nx*rz*ty + 3*nx*ry*tz)*pow(nz,2),
   -3*nz*(rz*ty + ry*tz)*pow(ny,2) + (ry*ty - rz*tz)*pow(ny,3) + 3*ny*(-(ry*ty) + rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3),
   3*nz*(ry*ty - rz*tz)*pow(ny,2) + (rz*ty + ry*tz)*pow(ny,3) - 3*ny*(rz*ty + ry*tz)*pow(nz,2) + (-(ry*ty) + rz*tz)*pow(nz,3),
   -(tx*pow(nz,2)*(pow(tx,2) - 3*pow(tz,2))) + 2*nx*nz*tz*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(pow(tx,3) - 3*tx*pow(tz,2)),
   2*nx*tx*(-6*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 3*ty*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 
    nz*(2*ny*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*nz*ty*(-pow(tx,2) + pow(tz,2))),
   2*nx*nz*tx*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(3*tz*pow(tx,2) - pow(tz,3)),
   -(tx*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
    2*nx*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,2)*(pow(tx,3) - 3*tx*pow(tz,2)) - 
    6*ny*ty*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))),6*ty*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))) + 
    2*ny*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)),
   -(ty*(12*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)))) + 
    3*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + ny*(6*nx*tx*(pow(ty,2) - pow(tz,2)) + 2*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))),
   6*ny*ty*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) + 6*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + 
    tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3)) + pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3)),
   2*ny*ty*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 3*tx*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + 
    nz*(2*nx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*nz*tx*(-pow(ty,2) + pow(tz,2))),
   6*tx*ty*tz*pow(ny,2) + 2*nz*ty*(-3*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
    ny*(6*nx*tz*pow(ty,2) + 6*nz*tx*(pow(ty,2) - pow(tz,2)) - 2*nx*pow(tz,3)),
   -(ty*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2))) + 2*ny*nz*tz*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,2)*(pow(ty,3) - 3*ty*pow(tz,2)),
   2*ny*nz*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,2)*(3*tz*pow(ty,2) - pow(tz,3)),
   pow(nx,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - 2*nx*nz*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    pow(nz,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),pow(nx,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    2*nx*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    nz*(nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) - 2*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(nx,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 2*nx*nz*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    pow(nz,2)*(-2*rx*tx*tz + rz*(-pow(tx,2) + pow(tz,2))),-(pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
    2*nx*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    2*ny*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
      nx*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
      ny*(nz*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + nx*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))),
   -4*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) - pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
    pow(nx,2)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) - 
    2*ny*(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + nx*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))),
   -(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + pow(ny,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 2*nx*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    2*ny*(2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   pow(ny,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) - 
    2*ny*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    nz*(nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) - 2*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) - nz*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
      ny*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))),
   pow(ny,2)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 2*ny*nz*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    pow(nz,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),pow(ny,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    2*ny*nz*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + pow(nz,2)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))),
   4*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + nx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   4*ty*(nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2))) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   2*(2*ny*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
      nx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(nz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(nx*ty*(3*pow(tx,2) - pow(tz,2)) + ny*(pow(tx,3) - tx*pow(tz,2)))),
   2*(2*ty*(nx*tx*(pow(ty,2) - 3*pow(tz,2)) + nz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
      ny*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   2*(6*ny*ty*tz*pow(tx,2) + 6*nx*tx*tz*pow(ty,2) - 2*nx*tx*pow(tz,3) - 2*ny*ty*pow(tz,3) + 
      nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*ny*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*(tz*(nx*ty*(pow(ty,2) - pow(tz,2)) + ny*tx*(3*pow(ty,2) - pow(tz,2))) + nz*tx*ty*(pow(ty,2) - 3*pow(tz,2))),
   4*nz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ny*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*ny*ty*tz*(pow(ty,2) - pow(tz,2)) + nz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   nx*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    nz*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    nx*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3)),
   nz*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    nx*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)),
   -(nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)))) + 
    nx*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    ny*(ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))),
   nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    nx*(6*rx*tx*ty*tz + 3*rz*ty*pow(tx,2) + 3*ry*tz*pow(tx,2) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   -(nz*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
    nx*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    ny*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))),
   nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    nx*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    ny*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))),
   ny*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
    nx*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    nz*(-6*ry*tx*ty*tz - 3*rx*tz*pow(ty,2) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)) + rx*pow(tz,3)),
   nz*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
    ny*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   ny*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    nz*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   nz*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    ny*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4),
   5*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5),
   10*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)),20*tx*ty*tz*(pow(tx,2) - pow(tz,2)),
   10*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)),
   2*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)),5*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   20*tx*ty*tz*(pow(ty,2) - pow(tz,2)),pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4),5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5),
   4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   4*ty*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   2*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
      rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))),
   2*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
      ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   2*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
      rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))),
   4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4));




    degree = 6;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P << pow(nx,6) - 15*pow(nx,4)*pow(nz,2) + 15*pow(nx,2)*pow(nz,4) - pow(nz,6),6*nx*ny*(pow(nx,4) - 10*pow(nx,2)*pow(nz,2) + 5*pow(nz,4)),
   6*nz*pow(nx,5) - 20*pow(nx,3)*pow(nz,3) + 6*nx*pow(nz,5),
   3*(5*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 10*pow(nx,2)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 5*pow(ny,2)*pow(nz,4) - pow(nz,6)),
   6*ny*nz*(5*pow(nx,4) - 10*pow(nx,2)*pow(nz,2) + pow(nz,4)),20*nx*ny*(pow(nx,2)*(pow(ny,2) - 3*pow(nz,2)) - 3*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)),
   4*nx*nz*(5*pow(nx,2)*(3*pow(ny,2) - pow(nz,2)) + 3*pow(nz,2)*(-5*pow(ny,2) + pow(nz,2))),
   3*(-(pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4))) + 5*pow(nx,2)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4))),
   4*ny*nz*(15*pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 5*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)),6*nx*ny*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)),
   6*nx*nz*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)),pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6),
   6*nz*pow(ny,5) - 20*pow(ny,3)*pow(nz,3) + 6*ny*pow(nz,5),
   -5*nz*tz*pow(nx,4) + tx*pow(nx,5) - 10*tx*pow(nx,3)*pow(nz,2) + 10*tz*pow(nx,2)*pow(nz,3) + 5*nx*tx*pow(nz,4) - tz*pow(nz,5),
   -10*nz*(nz*ty + 2*ny*tz)*pow(nx,3) + 5*ny*tx*pow(nx,4) + ty*pow(nx,5) - 30*ny*tx*pow(nx,2)*pow(nz,2) + 5*nx*(nz*ty + 4*ny*tz)*pow(nz,3) + 
    5*ny*tx*pow(nz,4),5*nz*tx*pow(nx,4) + tz*pow(nx,5) - 10*tz*pow(nx,3)*pow(nz,2) - 10*tx*pow(nx,2)*pow(nz,3) + 5*nx*tz*pow(nz,4) + tx*pow(nz,5),
   5*(ny*ty - nz*tz)*pow(nx,4) + 10*tx*pow(nx,3)*(pow(ny,2) - pow(nz,2)) + 10*nx*tx*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    10*nz*pow(nx,2)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + (5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3),
   20*ny*nz*tx*pow(nx,3) + 5*(nz*ty + ny*tz)*pow(nx,4) - 10*(nz*ty + 3*ny*tz)*pow(nx,2)*pow(nz,2) - 20*nx*ny*tx*pow(nz,3) + 
    (nz*ty + 5*ny*tz)*pow(nz,4),10*(ny*tx*pow(nx,2)*(pow(ny,2) - 3*pow(nz,2)) + ny*tx*pow(nz,2)*(-pow(ny,2) + pow(nz,2)) + 
      pow(nx,3)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + nx*nz*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3))),
   2*(-5*nz*tx*pow(nx,2)*(-3*pow(ny,2) + pow(nz,2)) + 5*pow(nx,3)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 
      5*nx*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + tx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3)),
   10*pow(nx,2)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) + 5*nx*tx*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    nz*(-10*nz*ty*pow(ny,3) - 5*tz*pow(ny,4) + 20*tz*pow(ny,2)*pow(nz,2) + 10*ny*ty*pow(nz,3) - 3*tz*pow(nz,4)),
   2*(10*nx*ny*nz*tx*(pow(ny,2) - pow(nz,2)) + 5*pow(nx,2)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3))),
   -10*nz*(nz*tx + 2*nx*tz)*pow(ny,3) + 5*nx*ty*pow(ny,4) + tx*pow(ny,5) - 30*nx*ty*pow(ny,2)*pow(nz,2) + 5*ny*(nz*tx + 4*nx*tz)*pow(nz,3) + 
    5*nx*ty*pow(nz,4),20*nx*nz*ty*pow(ny,3) + 5*(nz*tx + nx*tz)*pow(ny,4) - 10*(nz*tx + 3*nx*tz)*pow(ny,2)*pow(nz,2) - 20*nx*ny*ty*pow(nz,3) + 
    (nz*tx + 5*nx*tz)*pow(nz,4),-5*nz*tz*pow(ny,4) + ty*pow(ny,5) - 10*ty*pow(ny,3)*pow(nz,2) + 10*tz*pow(ny,2)*pow(nz,3) + 5*ny*ty*pow(nz,4) - 
    tz*pow(nz,5),5*nz*ty*pow(ny,4) + tz*pow(ny,5) - 10*tz*pow(ny,3)*pow(nz,2) - 10*ty*pow(ny,2)*pow(nz,3) + 5*ny*tz*pow(nz,4) + ty*pow(nz,5),
   -5*nz*rz*pow(nx,4) + rx*pow(nx,5) - 10*rx*pow(nx,3)*pow(nz,2) + 10*rz*pow(nx,2)*pow(nz,3) + 5*nx*rx*pow(nz,4) - rz*pow(nz,5),
   -10*nz*(nz*ry + 2*ny*rz)*pow(nx,3) + 5*ny*rx*pow(nx,4) + ry*pow(nx,5) - 30*ny*rx*pow(nx,2)*pow(nz,2) + 5*nx*(nz*ry + 4*ny*rz)*pow(nz,3) + 
    5*ny*rx*pow(nz,4),5*nz*rx*pow(nx,4) + rz*pow(nx,5) - 10*rz*pow(nx,3)*pow(nz,2) - 10*rx*pow(nx,2)*pow(nz,3) + 5*nx*rz*pow(nz,4) + rx*pow(nz,5),
   5*(ny*ry - nz*rz)*pow(nx,4) + 10*rx*pow(nx,3)*(pow(ny,2) - pow(nz,2)) + 10*nx*rx*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    10*nz*pow(nx,2)*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2)) + (5*ny*nz*ry + 10*rz*pow(ny,2) - 3*rz*pow(nz,2))*pow(nz,3),
   20*ny*nz*rx*pow(nx,3) + 5*(nz*ry + ny*rz)*pow(nx,4) - 10*(nz*ry + 3*ny*rz)*pow(nx,2)*pow(nz,2) - 20*nx*ny*rx*pow(nz,3) + 
    (nz*ry + 5*ny*rz)*pow(nz,4),10*(ny*rx*pow(nx,2)*(pow(ny,2) - 3*pow(nz,2)) + ny*rx*pow(nz,2)*(-pow(ny,2) + pow(nz,2)) + 
      pow(nx,3)*(-2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2)) + nx*nz*(-3*nz*ry*pow(ny,2) - 2*rz*pow(ny,3) + 4*ny*rz*pow(nz,2) + ry*pow(nz,3))),
   2*(-5*nz*rx*pow(nx,2)*(-3*pow(ny,2) + pow(nz,2)) + 5*pow(nx,3)*(2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2)) + 
      5*nx*pow(nz,2)*(-2*ny*nz*ry - 3*rz*pow(ny,2) + rz*pow(nz,2)) + rx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3)),
   10*pow(nx,2)*(-3*nz*rz*pow(ny,2) + ry*pow(ny,3) - 3*ny*ry*pow(nz,2) + rz*pow(nz,3)) + 5*nx*rx*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    nz*(-10*nz*ry*pow(ny,3) - 5*rz*pow(ny,4) + 20*rz*pow(ny,2)*pow(nz,2) + 10*ny*ry*pow(nz,3) - 3*rz*pow(nz,4)),
   2*(10*nx*ny*nz*rx*(pow(ny,2) - pow(nz,2)) + 5*pow(nx,2)*(3*nz*ry*pow(ny,2) + rz*pow(ny,3) - 3*ny*rz*pow(nz,2) - ry*pow(nz,3)) + 
      pow(nz,2)*(-5*nz*ry*pow(ny,2) - 5*rz*pow(ny,3) + 5*ny*rz*pow(nz,2) + ry*pow(nz,3))),
   -10*nz*(nz*rx + 2*nx*rz)*pow(ny,3) + 5*nx*ry*pow(ny,4) + rx*pow(ny,5) - 30*nx*ry*pow(ny,2)*pow(nz,2) + 5*ny*(nz*rx + 4*nx*rz)*pow(nz,3) + 
    5*nx*ry*pow(nz,4),20*nx*nz*ry*pow(ny,3) + 5*(nz*rx + nx*rz)*pow(ny,4) - 10*(nz*rx + 3*nx*rz)*pow(ny,2)*pow(nz,2) - 20*nx*ny*ry*pow(nz,3) + 
    (nz*rx + 5*nx*rz)*pow(nz,4),-5*nz*rz*pow(ny,4) + ry*pow(ny,5) - 10*ry*pow(ny,3)*pow(nz,2) + 10*rz*pow(ny,2)*pow(nz,3) + 5*ny*ry*pow(nz,4) - 
    rz*pow(nz,5),5*nz*ry*pow(ny,4) + rz*pow(ny,5) - 10*rz*pow(ny,3)*pow(nz,2) - 10*ry*pow(ny,2)*pow(nz,3) + 5*ny*rz*pow(nz,4) + ry*pow(nz,5),
   -8*nz*tx*tz*pow(nx,3) + 8*nx*tx*tz*pow(nz,3) + pow(nx,4)*(pow(tx,2) - pow(tz,2)) + pow(nz,4)*(pow(tx,2) - pow(tz,2)) + 
    6*pow(nx,2)*pow(nz,2)*(-pow(tx,2) + pow(tz,2)),2*(-6*nz*tx*(nz*ty + 2*ny*tz)*pow(nx,2) + tx*ty*pow(nx,4) + tx*(nz*ty + 4*ny*tz)*pow(nz,3) + 
      2*nx*pow(nz,2)*(2*nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 2*pow(nx,3)*(2*nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   2*(tx*tz*pow(nx,4) - 6*tx*tz*pow(nx,2)*pow(nz,2) + tx*tz*pow(nz,4) + 2*nz*pow(nx,3)*(pow(tx,2) - pow(tz,2)) + 
      2*nx*pow(nz,3)*(-pow(tx,2) + pow(tz,2))),8*tx*(ny*ty - nz*tz)*pow(nx,3) + 8*nx*nz*tx*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 
    pow(nz,2)*(8*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    6*pow(nx,2)*(-4*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,4)*(pow(ty,2) - pow(tz,2)),2*(4*tx*(nz*ty + ny*tz)*pow(nx,3) + ty*tz*pow(nx,4) - 4*nx*tx*(nz*ty + 3*ny*tz)*pow(nz,2) + 
      pow(nz,3)*(nz*ty*tz - 2*ny*pow(tx,2) + 2*ny*pow(tz,2)) - 6*nz*pow(nx,2)*(nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   4*(3*tx*pow(nx,2)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + nz*tx*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      nx*(-6*nz*ty*tz*pow(ny,2) + 4*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 
      pow(nx,3)*(-2*nz*ty*tz + ny*(pow(ty,2) - pow(tz,2)))),
   4*(3*tx*pow(nx,2)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + tx*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) - 
      nx*nz*(6*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
      pow(nx,3)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2)))),8*ty*(nx*tx - nz*tz)*pow(ny,3) + 
    8*ny*nz*ty*(-3*nx*nz*tx - 3*tz*pow(nx,2) + 2*tz*pow(nz,2)) + pow(ny,4)*(pow(tx,2) - pow(tz,2)) - 
    6*pow(ny,2)*(4*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    pow(nz,2)*(8*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 6*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   4*(3*ty*pow(ny,2)*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + ty*pow(nz,2)*(-2*nx*nz*tx - 3*tz*pow(nx,2) + tz*pow(nz,2)) + 
      pow(ny,3)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
      ny*nz*(6*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 3*pow(nx,2)*(-pow(ty,2) + pow(tz,2)))),
   2*(-6*nz*ty*(nz*tx + 2*nx*tz)*pow(ny,2) + tx*ty*pow(ny,4) + ty*(nz*tx + 4*nx*tz)*pow(nz,3) + 
      2*ny*pow(nz,2)*(2*nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 2*pow(ny,3)*(2*nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   2*(4*ty*(nz*tx + nx*tz)*pow(ny,3) + tx*tz*pow(ny,4) - 4*ny*ty*(nz*tx + 3*nx*tz)*pow(nz,2) + 
      pow(nz,3)*(nz*tx*tz - 2*nx*pow(ty,2) + 2*nx*pow(tz,2)) - 6*nz*pow(ny,2)*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   -8*nz*ty*tz*pow(ny,3) + 8*ny*ty*tz*pow(nz,3) + pow(ny,4)*(pow(ty,2) - pow(tz,2)) + pow(nz,4)*(pow(ty,2) - pow(tz,2)) + 
    6*pow(ny,2)*pow(nz,2)*(-pow(ty,2) + pow(tz,2)),2*(ty*tz*pow(ny,4) - 6*ty*tz*pow(ny,2)*pow(nz,2) + ty*tz*pow(nz,4) + 
      2*nz*pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 2*ny*pow(nz,3)*(-pow(ty,2) + pow(tz,2))),
   -4*nz*(rz*tx + rx*tz)*pow(nx,3) + (rx*tx - rz*tz)*pow(nx,4) + 6*(-(rx*tx) + rz*tz)*pow(nx,2)*pow(nz,2) + 4*nx*(rz*tx + rx*tz)*pow(nz,3) + 
    (rx*tx - rz*tz)*pow(nz,4),-6*nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz)*pow(nx,2) - 
    4*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,3) + (ry*tx + rx*ty)*pow(nx,4) + 
    4*nx*(-3*ny*rx*tx + nz*rz*ty + nz*ry*tz + 3*ny*rz*tz)*pow(nz,2) + (nz*ry*tx + 4*ny*rz*tx + nz*rx*ty + 4*ny*rx*tz)*pow(nz,3),
   4*nz*(rx*tx - rz*tz)*pow(nx,3) + (rz*tx + rx*tz)*pow(nx,4) - 6*(rz*tx + rx*tz)*pow(nx,2)*pow(nz,2) + 4*nx*(-(rx*tx) + rz*tz)*pow(nz,3) + 
    (rz*tx + rx*tz)*pow(nz,4),4*(ny*(ry*tx + rx*ty) - nz*(rz*tx + rx*tz))*pow(nx,3) + (ry*ty - rz*tz)*pow(nx,4) - 
    4*nx*nz*(3*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(nz,2)) + 
    pow(nz,2)*(4*ny*nz*(rz*ty + ry*tz) + (-6*rx*tx + 6*rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2)) + 
    6*pow(nx,2)*(-2*ny*nz*(rz*ty + ry*tz) + (rx*tx - rz*tz)*pow(ny,2) - (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)),
   -6*nz*(-2*ny*rx*tx + nz*rz*ty + nz*ry*tz + 2*ny*rz*tz)*pow(nx,2) + 4*(nz*ry*tx + ny*rz*tx + nz*rx*ty + ny*rx*tz)*pow(nx,3) + 
    (rz*ty + ry*tz)*pow(nx,4) - 4*nx*(nz*ry*tx + 3*ny*rz*tx + nz*rx*ty + 3*ny*rx*tz)*pow(nz,2) + 
    (-4*ny*rx*tx + nz*rz*ty + nz*ry*tz + 4*ny*rz*tz)*pow(nz,3),
   2*(-2*(-(ny*ry*ty) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,3) + 
      3*pow(nx,2)*(-2*ny*nz*(rz*tx + rx*tz) + (ry*tx + rx*ty)*pow(ny,2) - (ry*tx + rx*ty)*pow(nz,2)) + 
      nz*(-3*nz*(ry*tx + rx*ty)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(ny,3) + 4*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
      2*nx*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (rx*tx - rz*tz)*pow(ny,3) - 3*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + 2*(rz*ty + ry*tz)*pow(nz,3))),
   2*(2*(nz*ry*ty + ny*rz*ty + ny*ry*tz - nz*rz*tz)*pow(nx,3) + 
      3*pow(nx,2)*(2*ny*nz*(ry*tx + rx*ty) + (rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) + 
      pow(nz,2)*(-2*ny*nz*(ry*tx + rx*ty) - 3*(rz*tx + rx*tz)*pow(ny,2) + (rz*tx + rx*tz)*pow(nz,2)) - 
      2*nx*nz*(3*ny*nz*(rz*ty + ry*tz) + (-3*rx*tx + 3*rz*tz)*pow(ny,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2))),
   4*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz))*pow(ny,3) + (rx*tx - rz*tz)*pow(ny,4) - 
    4*ny*nz*(3*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - 2*(rz*ty + ry*tz)*pow(nz,2)) + 
    pow(nz,2)*(4*nx*nz*(rz*tx + rx*tz) + 6*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2)) - 
    6*pow(ny,2)*(2*nx*nz*(rz*tx + rx*tz) + (-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)),
   2*(2*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz)*pow(ny,3) + 
      3*pow(ny,2)*(2*nx*nz*(ry*tx + rx*ty) + (rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) + 
      pow(nz,2)*(-2*nx*nz*(ry*tx + rx*ty) - 3*(rz*ty + ry*tz)*pow(nx,2) + (rz*ty + ry*tz)*pow(nz,2)) - 
      2*ny*nz*(3*nx*nz*(rz*tx + rx*tz) + 3*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2))),
   -6*nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz)*pow(ny,2) - 4*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz)*pow(ny,3) + 
    (ry*tx + rx*ty)*pow(ny,4) + 4*ny*(nz*rz*tx - 3*nx*ry*ty + nz*rx*tz + 3*nx*rz*tz)*pow(nz,2) + 
    (nz*ry*tx + nz*rx*ty + 4*nx*rz*ty + 4*nx*ry*tz)*pow(nz,3),
   -6*nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz)*pow(ny,2) + 4*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz)*pow(ny,3) + 
    (rz*tx + rx*tz)*pow(ny,4) - 4*ny*(nz*ry*tx + nz*rx*ty + 3*nx*rz*ty + 3*nx*ry*tz)*pow(nz,2) + 
    (nz*rz*tx - 4*nx*ry*ty + nz*rx*tz + 4*nx*rz*tz)*pow(nz,3),
   -4*nz*(rz*ty + ry*tz)*pow(ny,3) + (ry*ty - rz*tz)*pow(ny,4) + 6*(-(ry*ty) + rz*tz)*pow(ny,2)*pow(nz,2) + 4*ny*(rz*ty + ry*tz)*pow(nz,3) + 
    (ry*ty - rz*tz)*pow(nz,4),4*nz*(ry*ty - rz*tz)*pow(ny,3) + (rz*ty + ry*tz)*pow(ny,4) - 6*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,2) + 
    4*ny*(-(ry*ty) + rz*tz)*pow(nz,3) + (rz*ty + ry*tz)*pow(nz,4),
   -3*nx*tx*pow(nz,2)*(pow(tx,2) - 3*pow(tz,2)) + 3*nz*tz*pow(nx,2)*(-3*pow(tx,2) + pow(tz,2)) - tz*pow(nz,3)*(-3*pow(tx,2) + pow(tz,2)) + 
    pow(nx,3)*(pow(tx,3) - 3*tx*pow(tz,2)),3*(tx*pow(nx,2)*(-6*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + ty*pow(nx,3)*(pow(tx,2) - pow(tz,2)) + 
      tx*pow(nz,2)*(2*nz*ty*tz - ny*pow(tx,2) + 3*ny*pow(tz,2)) + nx*nz*(2*ny*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*nz*ty*(-pow(tx,2) + pow(tz,2)))),
   3*nz*tx*pow(nx,2)*(pow(tx,2) - 3*pow(tz,2)) - tx*pow(nz,3)*(pow(tx,2) - 3*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + 
    pow(nx,3)*(3*tz*pow(tx,2) - pow(tz,3)),3*(nx*tx*(-12*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
         pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + tx*pow(nx,3)*(pow(ty,2) - pow(tz,2)) + 
      pow(nx,2)*(3*ny*ty*(pow(tx,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      nz*(tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 3*ny*nz*ty*(-pow(tx,2) + pow(tz,2)) + pow(ny,2)*(-3*tz*pow(tx,2) + pow(tz,3)))),
   3*(2*tx*ty*tz*pow(nx,3) + 2*nx*nz*tx*(-3*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 
      pow(nx,2)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2))) + 
      pow(nz,2)*(ny*tz*(-3*pow(tx,2) + pow(tz,2)) + nz*ty*(-pow(tx,2) + pow(tz,2)))),
   ty*(-18*nz*tx*tz*pow(nx,2) + 12*tx*tz*pow(nz,3) - 3*nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(ny,3)*(pow(tx,3) - 3*tx*pow(tz,2)) - 9*ty*pow(ny,2)*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))) - 
    3*ny*(tx*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 2*nx*nz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       3*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),18*ny*ty*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))) - 
    tx*pow(nz,3)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 9*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) - tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 
    3*nx*tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + 3*pow(ny,2)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)),
   3*(-(ny*ty*(12*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)))) + 
      ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + pow(ny,2)*(3*nx*tx*(pow(ty,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      nz*(tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 3*nx*nz*tx*(-pow(ty,2) + pow(tz,2)) + pow(nx,2)*(-3*tz*pow(ty,2) + pow(tz,3)))),
   -(nz*ty*(18*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 3*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)))) + 
    9*ty*pow(ny,2)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) + 
    3*ny*(6*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + 
    pow(ny,3)*(3*tz*pow(tx,2) - pow(tz,3)),3*(ty*pow(ny,2)*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + tx*pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 
      ty*pow(nz,2)*(2*nz*tx*tz - nx*pow(ty,2) + 3*nx*pow(tz,2)) + ny*nz*(2*nx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*nz*tx*(-pow(ty,2) + pow(tz,2)))),
   3*(2*tx*ty*tz*pow(ny,3) + 2*ny*nz*ty*(-3*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
      pow(ny,2)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2))) + 
      pow(nz,2)*(nx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*tx*(-pow(ty,2) + pow(tz,2)))),
   -3*ny*ty*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*nz*tz*pow(ny,2)*(-3*pow(ty,2) + pow(tz,2)) - tz*pow(nz,3)*(-3*pow(ty,2) + pow(tz,2)) + 
    pow(ny,3)*(pow(ty,3) - 3*ty*pow(tz,2)),3*nz*ty*pow(ny,2)*(pow(ty,2) - 3*pow(tz,2)) - ty*pow(nz,3)*(pow(ty,2) - 3*pow(tz,2)) + 
    3*ny*tz*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,3)*(3*tz*pow(ty,2) - pow(tz,3)),
   pow(nx,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - 3*nz*pow(nx,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    pow(nz,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 3*nx*pow(nz,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   pow(nx,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    3*pow(nx,2)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) - 
    3*nx*nz*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 2*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(nx,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 3*nx*pow(nz,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
    3*nz*pow(nx,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + pow(nz,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   3*nx*(-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) - pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + pow(nx,3)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    3*pow(nx,2)*(-(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       ny*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    nz*(3*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
    ,2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,3) - 6*nx*nz*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*rz*tx*tz - rx*pow(tx,2) + rx*pow(tz,2))) + 
    pow(nz,2)*(nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) - 3*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
    3*pow(nx,2)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   -6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
    3*nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 3*ny*
     (pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))
       ) + 3*pow(ny,2)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   -(pow(nz,3)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
    3*nx*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    3*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    6*ny*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
       nx*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    3*pow(ny,2)*(nz*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + nx*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) - 
    3*pow(ny,2)*(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       nx*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    nz*(3*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 3*ny*(4*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 3*pow(ny,2)*
     (2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    nz*(6*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))) - 
    3*ny*(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,2)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,3)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) - 
    3*pow(ny,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 3*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
    3*ny*nz*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 2*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) - 6*ny*nz*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    pow(nz,2)*(nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) - 3*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
    3*pow(ny,2)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   pow(ny,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 3*nz*pow(ny,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    pow(nz,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 3*ny*pow(nz,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),
   pow(ny,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 3*ny*pow(nz,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    3*nz*pow(ny,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + pow(nz,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),
   8*nx*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + pow(nx,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    pow(nz,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   2*(-2*nz*tx*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*ny*tz*(pow(tx,2) - pow(tz,2))) + 2*tx*ty*pow(nx,2)*(pow(tx,2) - 3*pow(tz,2)) + 
      nx*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   2*(2*tx*tz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 2*tx*tz*pow(nz,2)*(-pow(tx,2) + pow(tz,2)) + 
      nx*nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   -8*nx*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 8*ny*ty*(nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    2*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),
   2*(2*ty*(2*nx*nz*tx*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
      ny*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*(-(ty*(tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
           2*nx*nz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + tx*ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) + 
      ny*(-2*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + nx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   4*(-(tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + tx*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) - 
      tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 2*ny*ty*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)) + 
      nx*nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 8*ny*ty*(nx*tx*(pow(ty,2) - 3*pow(tz,2)) + nz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
    pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    2*pow(ny,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),
   4*(-(ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2))) + ty*(2*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
      ny*(-2*nx*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(-2*nz*ty*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 2*nx*tz*(pow(ty,2) - pow(tz,2))) + 2*tx*ty*pow(ny,2)*(pow(ty,2) - 3*pow(tz,2)) + 
      ny*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(4*ny*ty*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + nx*tz*(pow(ty,2) - pow(tz,2))) + tx*pow(ny,2)*(6*tz*pow(ty,2) - 2*pow(tz,3)) + 
      nz*(2*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   8*ny*nz*ty*tz*(-pow(ty,2) + pow(tz,2)) + pow(ny,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    pow(nz,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   2*(2*ty*tz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + 2*ty*tz*pow(nz,2)*(-pow(ty,2) + pow(tz,2)) + 
      ny*nz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   pow(nx,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    pow(nz,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    2*nx*nz*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   pow(nx,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    nz*(nz*(6*rz*tx*ty*tz - 3*rx*ty*pow(tx,2) - ry*pow(tx,3) + 3*ry*tx*pow(tz,2) + 3*rx*ty*pow(tz,2)) + 
       2*ny*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3))) + 
    2*nx*(ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   2*nx*nz*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nx,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    pow(nz,2)*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   -2*nx*nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    pow(nz,2)*(-(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
       3*ry*ty*(-pow(tx,2) + pow(tz,2))) + pow(nx,2)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
       rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    2*ny*(nx*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(nx,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    2*nx*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nz,2)*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3)) + 
    2*ny*(nx*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
       nz*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3))),
   -2*nx*nz*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) - 
    pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    pow(nx,2)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    pow(ny,2)*(ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))) - 
    2*ny*(nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       nx*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2)))),
   -(pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)))) + 
    2*nx*nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    2*ny*(nx*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)),
   pow(nz,2)*(-(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
    pow(ny,2)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    2*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    2*ny*(nx*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))),
   -(pow(nz,2)*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
    2*nx*nz*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    2*ny*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       nx*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)))) + 
    pow(nx,2)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   pow(ny,2)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
    nz*(nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       2*nx*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))) - 
    2*ny*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    nz*(2*nx*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       nz*(-6*ry*tx*ty*tz - 3*rx*tz*pow(ty,2) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)) + rx*pow(tz,3))) + 
    2*ny*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))),
   pow(ny,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    pow(nz,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    2*ny*nz*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   2*ny*nz*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(ny,2)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    pow(nz,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   -(nz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + nx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   ny*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
    5*ty*(4*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + nx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
    5*ny*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*nx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    tz*(20*nx*tx*ty*(pow(tx,2) - pow(tz,2)) + ny*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   10*(ny*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      ty*(-2*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + nx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(5*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      tz*(10*ny*tx*ty*(pow(tx,2) - pow(tz,2)) + nx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   nz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
    5*nx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ny*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)),
   2*(5*nz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      tz*(10*nx*tx*ty*(pow(ty,2) - pow(tz,2)) + ny*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   5*ny*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    ty*(20*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    tz*(20*ny*tx*ty*(pow(ty,2) - pow(tz,2)) + nx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   -(nz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ny*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   ny*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   nx*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    nz*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   -4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
    ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    nx*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   nz*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    nx*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   ny*(4*ty*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*nx*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
       rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   4*nx*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
    nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   2*(-2*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
         rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ny*
       (2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nx*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(2*ny*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))) + 
      nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nx*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
         rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   nx*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*ny*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   2*(2*nx*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
      nz*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
         rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   -4*nz*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    ny*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*ny*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   ny*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    nz*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   nz*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    ny*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6),6*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)),
   6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5),
   3*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   6*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)),20*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),
   4*tx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))),
   3*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*ty*tz*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),6*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)),
   6*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)),pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6),
   6*tz*pow(ty,5) - 20*pow(ty,3)*pow(tz,3) + 6*ty*pow(tz,5),
   -(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
    5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
    5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   10*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
    5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)),
   2*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   -(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4));

    

    degree = 7;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    tensor_project[degree].P <<pow(nx,7) - 21*pow(nx,5)*pow(nz,2) + 35*pow(nx,3)*pow(nz,4) - 7*nx*pow(nz,6),
   7*ny*(pow(nx,6) - 15*pow(nx,4)*pow(nz,2) + 15*pow(nx,2)*pow(nz,4) - pow(nz,6)),
   7*nz*pow(nx,6) - 35*pow(nx,4)*pow(nz,3) + 21*pow(nx,2)*pow(nz,5) - pow(nz,7),
   7*nx*(3*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 10*pow(nx,2)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 15*pow(ny,2)*pow(nz,4) - 3*pow(nz,6)),
   14*nx*ny*nz*(3*pow(nx,4) - 10*pow(nx,2)*pow(nz,2) + 3*pow(nz,4)),
   7*ny*(5*pow(nx,4)*(pow(ny,2) - 3*pow(nz,2)) + 30*pow(nx,2)*pow(nz,2)*(-pow(ny,2) + pow(nz,2)) + 5*pow(ny,2)*pow(nz,4) - 3*pow(nz,6)),
   35*pow(nx,4)*(3*nz*pow(ny,2) - pow(nz,3)) - 42*pow(nx,2)*(5*pow(ny,2)*pow(nz,3) - pow(nz,5)) + 21*pow(ny,2)*pow(nz,5) - 3*pow(nz,7),
   7*nx*(5*pow(nx,2)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 3*(5*pow(ny,4)*pow(nz,2) - 10*pow(ny,2)*pow(nz,4) + pow(nz,6))),
   28*nx*ny*nz*(5*pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 5*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)),
   7*ny*(-3*pow(ny,4)*pow(nz,2) + 10*pow(ny,2)*pow(nz,4) + 3*pow(nx,2)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) - 3*pow(nz,6)),
   -35*pow(ny,4)*pow(nz,3) + 42*pow(ny,2)*pow(nz,5) + 21*pow(nx,2)*(5*nz*pow(ny,4) - 10*pow(ny,2)*pow(nz,3) + pow(nz,5)) - 3*pow(nz,7),
   7*nx*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)),14*nx*ny*nz*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)),
   pow(ny,7) - 21*pow(ny,5)*pow(nz,2) + 35*pow(ny,3)*pow(nz,4) - 7*ny*pow(nz,6),
   7*nz*pow(ny,6) - 35*pow(ny,4)*pow(nz,3) + 21*pow(ny,2)*pow(nz,5) - pow(nz,7),
   -6*nz*tz*pow(nx,5) + tx*pow(nx,6) - 15*tx*pow(nx,4)*pow(nz,2) + 20*tz*pow(nx,3)*pow(nz,3) + 15*tx*pow(nx,2)*pow(nz,4) - 6*nx*tz*pow(nz,5) - 
    tx*pow(nz,6),-15*nz*(nz*ty + 2*ny*tz)*pow(nx,4) + 6*ny*tx*pow(nx,5) + ty*pow(nx,6) - 60*ny*tx*pow(nx,3)*pow(nz,2) + 
    15*(nz*ty + 4*ny*tz)*pow(nx,2)*pow(nz,3) + 30*nx*ny*tx*pow(nz,4) - (nz*ty + 6*ny*tz)*pow(nz,5),
   6*nz*tx*pow(nx,5) + tz*pow(nx,6) - 15*tz*pow(nx,4)*pow(nz,2) - 20*tx*pow(nx,3)*pow(nz,3) + 15*tz*pow(nx,2)*pow(nz,4) + 6*nx*tx*pow(nz,5) - 
    tz*pow(nz,6),6*(ny*ty - nz*tz)*pow(nx,5) + 15*tx*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 30*tx*pow(nx,2)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    20*nz*pow(nx,3)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 6*nx*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) - 
    3*tx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4),30*ny*nz*tx*pow(nx,4) + 6*(nz*ty + ny*tz)*pow(nx,5) - 20*(nz*ty + 3*ny*tz)*pow(nx,3)*pow(nz,2) - 
    60*ny*tx*pow(nx,2)*pow(nz,3) + 6*nx*(nz*ty + 5*ny*tz)*pow(nz,4) + 6*ny*tx*pow(nz,5),
   20*ny*tx*pow(nx,3)*(pow(ny,2) - 3*pow(nz,2)) - 60*nx*ny*tx*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
    15*pow(nx,4)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
    pow(nz,3)*(15*nz*ty*pow(ny,2) + 20*tz*pow(ny,3) - 18*ny*tz*pow(nz,2) - 3*ty*pow(nz,3)) + 
    30*nz*pow(nx,2)*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)),
   -20*nz*tx*pow(nx,3)*(-3*pow(ny,2) + pow(nz,2)) + 15*pow(nx,4)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 
    30*pow(nx,2)*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 12*nx*tx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
    3*(2*ny*nz*ty + 5*tz*pow(ny,2) - tz*pow(nz,2))*pow(nz,4),
   20*pow(nx,3)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
    3*tx*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 15*tx*pow(nx,2)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    6*nx*nz*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)),
   4*(15*ny*nz*tx*pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 5*tx*pow(ny,3)*pow(nz,3) + 
      5*pow(nx,3)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      3*nx*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 3*ny*tx*pow(nz,5)),
   6*nx*ny*tx*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
    15*pow(nx,2)*(-4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4)) - 
    nz*(15*nz*ty*pow(ny,4) + 6*tz*pow(ny,5) - 40*tz*pow(ny,3)*pow(nz,2) - 30*ty*pow(ny,2)*pow(nz,3) + 18*ny*tz*pow(nz,4) + 3*ty*pow(nz,5)),
   6*nx*nz*tx*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    pow(nz,2)*(-20*nz*ty*pow(ny,3) - 15*tz*pow(ny,4) + 30*tz*pow(ny,2)*pow(nz,2) + 12*ny*ty*pow(nz,3) - 3*tz*pow(nz,4)) + 
    15*pow(nx,2)*(4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4)),
   -15*nz*(nz*tx + 2*nx*tz)*pow(ny,4) + 6*nx*ty*pow(ny,5) + tx*pow(ny,6) - 60*nx*ty*pow(ny,3)*pow(nz,2) + 15*(nz*tx + 4*nx*tz)*pow(ny,2)*pow(nz,3) + 
    30*nx*ny*ty*pow(nz,4) - (nz*tx + 6*nx*tz)*pow(nz,5),30*nx*nz*ty*pow(ny,4) + 6*(nz*tx + nx*tz)*pow(ny,5) - 
    20*(nz*tx + 3*nx*tz)*pow(ny,3)*pow(nz,2) - 60*nx*ty*pow(ny,2)*pow(nz,3) + 6*ny*(nz*tx + 5*nx*tz)*pow(nz,4) + 6*nx*ty*pow(nz,5),
   -6*nz*tz*pow(ny,5) + ty*pow(ny,6) - 15*ty*pow(ny,4)*pow(nz,2) + 20*tz*pow(ny,3)*pow(nz,3) + 15*ty*pow(ny,2)*pow(nz,4) - 6*ny*tz*pow(nz,5) - 
    ty*pow(nz,6),6*nz*ty*pow(ny,5) + tz*pow(ny,6) - 15*tz*pow(ny,4)*pow(nz,2) - 20*ty*pow(ny,3)*pow(nz,3) + 15*tz*pow(ny,2)*pow(nz,4) + 
    6*ny*ty*pow(nz,5) - tz*pow(nz,6),-6*nz*rz*pow(nx,5) + rx*pow(nx,6) - 15*rx*pow(nx,4)*pow(nz,2) + 20*rz*pow(nx,3)*pow(nz,3) + 
    15*rx*pow(nx,2)*pow(nz,4) - 6*nx*rz*pow(nz,5) - rx*pow(nz,6),
   -15*nz*(nz*ry + 2*ny*rz)*pow(nx,4) + 6*ny*rx*pow(nx,5) + ry*pow(nx,6) - 60*ny*rx*pow(nx,3)*pow(nz,2) + 15*(nz*ry + 4*ny*rz)*pow(nx,2)*pow(nz,3) + 
    30*nx*ny*rx*pow(nz,4) - (nz*ry + 6*ny*rz)*pow(nz,5),6*nz*rx*pow(nx,5) + rz*pow(nx,6) - 15*rz*pow(nx,4)*pow(nz,2) - 20*rx*pow(nx,3)*pow(nz,3) + 
    15*rz*pow(nx,2)*pow(nz,4) + 6*nx*rx*pow(nz,5) - rz*pow(nz,6),
   6*(ny*ry - nz*rz)*pow(nx,5) + 15*rx*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 30*rx*pow(nx,2)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    20*nz*pow(nx,3)*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2)) + 6*nx*(5*ny*nz*ry + 10*rz*pow(ny,2) - 3*rz*pow(nz,2))*pow(nz,3) - 
    3*rx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4),30*ny*nz*rx*pow(nx,4) + 6*(nz*ry + ny*rz)*pow(nx,5) - 20*(nz*ry + 3*ny*rz)*pow(nx,3)*pow(nz,2) - 
    60*ny*rx*pow(nx,2)*pow(nz,3) + 6*nx*(nz*ry + 5*ny*rz)*pow(nz,4) + 6*ny*rx*pow(nz,5),
   20*ny*rx*pow(nx,3)*(pow(ny,2) - 3*pow(nz,2)) - 60*nx*ny*rx*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
    15*pow(nx,4)*(-2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2)) + 
    pow(nz,3)*(15*nz*ry*pow(ny,2) + 20*rz*pow(ny,3) - 18*ny*rz*pow(nz,2) - 3*ry*pow(nz,3)) + 
    30*nz*pow(nx,2)*(-3*nz*ry*pow(ny,2) - 2*rz*pow(ny,3) + 4*ny*rz*pow(nz,2) + ry*pow(nz,3)),
   -20*nz*rx*pow(nx,3)*(-3*pow(ny,2) + pow(nz,2)) + 15*pow(nx,4)*(2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2)) + 
    30*pow(nx,2)*pow(nz,2)*(-2*ny*nz*ry - 3*rz*pow(ny,2) + rz*pow(nz,2)) + 12*nx*rx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
    3*(2*ny*nz*ry + 5*rz*pow(ny,2) - rz*pow(nz,2))*pow(nz,4),
   20*pow(nx,3)*(-3*nz*rz*pow(ny,2) + ry*pow(ny,3) - 3*ny*ry*pow(nz,2) + rz*pow(nz,3)) - 
    3*rx*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 15*rx*pow(nx,2)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    6*nx*nz*(10*nz*ry*pow(ny,3) + 5*rz*pow(ny,4) - 20*rz*pow(ny,2)*pow(nz,2) - 10*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)),
   4*(15*ny*nz*rx*pow(nx,2)*(pow(ny,2) - pow(nz,2)) - 5*rx*pow(ny,3)*pow(nz,3) + 
      5*pow(nx,3)*(3*nz*ry*pow(ny,2) + rz*pow(ny,3) - 3*ny*rz*pow(nz,2) - ry*pow(nz,3)) + 
      3*nx*pow(nz,2)*(-5*nz*ry*pow(ny,2) - 5*rz*pow(ny,3) + 5*ny*rz*pow(nz,2) + ry*pow(nz,3)) + 3*ny*rx*pow(nz,5)),
   6*nx*ny*rx*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
    15*pow(nx,2)*(-4*nz*rz*pow(ny,3) + ry*pow(ny,4) - 6*ry*pow(ny,2)*pow(nz,2) + 4*ny*rz*pow(nz,3) + ry*pow(nz,4)) - 
    nz*(15*nz*ry*pow(ny,4) + 6*rz*pow(ny,5) - 40*rz*pow(ny,3)*pow(nz,2) - 30*ry*pow(ny,2)*pow(nz,3) + 18*ny*rz*pow(nz,4) + 3*ry*pow(nz,5)),
   6*nx*nz*rx*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    pow(nz,2)*(-20*nz*ry*pow(ny,3) - 15*rz*pow(ny,4) + 30*rz*pow(ny,2)*pow(nz,2) + 12*ny*ry*pow(nz,3) - 3*rz*pow(nz,4)) + 
    15*pow(nx,2)*(4*nz*ry*pow(ny,3) + rz*pow(ny,4) - 6*rz*pow(ny,2)*pow(nz,2) - 4*ny*ry*pow(nz,3) + rz*pow(nz,4)),
   -15*nz*(nz*rx + 2*nx*rz)*pow(ny,4) + 6*nx*ry*pow(ny,5) + rx*pow(ny,6) - 60*nx*ry*pow(ny,3)*pow(nz,2) + 15*(nz*rx + 4*nx*rz)*pow(ny,2)*pow(nz,3) + 
    30*nx*ny*ry*pow(nz,4) - (nz*rx + 6*nx*rz)*pow(nz,5),30*nx*nz*ry*pow(ny,4) + 6*(nz*rx + nx*rz)*pow(ny,5) - 
    20*(nz*rx + 3*nx*rz)*pow(ny,3)*pow(nz,2) - 60*nx*ry*pow(ny,2)*pow(nz,3) + 6*ny*(nz*rx + 5*nx*rz)*pow(nz,4) + 6*nx*ry*pow(nz,5),
   -6*nz*rz*pow(ny,5) + ry*pow(ny,6) - 15*ry*pow(ny,4)*pow(nz,2) + 20*rz*pow(ny,3)*pow(nz,3) + 15*ry*pow(ny,2)*pow(nz,4) - 6*ny*rz*pow(nz,5) - 
    ry*pow(nz,6),6*nz*ry*pow(ny,5) + rz*pow(ny,6) - 15*rz*pow(ny,4)*pow(nz,2) - 20*ry*pow(ny,3)*pow(nz,3) + 15*rz*pow(ny,2)*pow(nz,4) + 
    6*ny*ry*pow(nz,5) - rz*pow(nz,6),-10*nz*tx*tz*pow(nx,4) + 20*tx*tz*pow(nx,2)*pow(nz,3) - 2*tx*tz*pow(nz,5) + pow(nx,5)*(pow(tx,2) - pow(tz,2)) - 
    10*pow(nx,3)*pow(nz,2)*(pow(tx,2) - pow(tz,2)) + 5*nx*pow(nz,4)*(pow(tx,2) - pow(tz,2)),
   -20*nz*tx*(nz*ty + 2*ny*tz)*pow(nx,3) + 2*tx*ty*pow(nx,5) + 10*nx*tx*(nz*ty + 4*ny*tz)*pow(nz,3) + 
    pow(nz,4)*(-2*nz*ty*tz + 5*ny*(pow(tx,2) - pow(tz,2))) + 10*pow(nx,2)*pow(nz,2)*(2*nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 
    5*pow(nx,4)*(2*nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2))),2*tx*tz*pow(nx,5) - 20*tx*tz*pow(nx,3)*pow(nz,2) + 10*nx*tx*tz*pow(nz,4) + 
    5*nz*pow(nx,4)*(pow(tx,2) - pow(tz,2)) - 10*pow(nx,2)*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + pow(nz,5)*(pow(tx,2) - pow(tz,2)),
   10*tx*(ny*ty - nz*tz)*pow(nx,4) + 20*nz*tx*pow(nx,2)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 
    2*tx*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) + 
    5*nx*pow(nz,2)*(8*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    10*pow(nx,3)*(-4*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,5)*(pow(ty,2) - pow(tz,2)),2*(5*tx*(nz*ty + ny*tz)*pow(nx,4) + ty*tz*pow(nx,5) - 10*tx*(nz*ty + 3*ny*tz)*pow(nx,2)*pow(nz,2) + 
      tx*(nz*ty + 5*ny*tz)*pow(nz,4) + 5*nx*pow(nz,3)*(nz*ty*tz - 2*ny*pow(tx,2) + 2*ny*pow(tz,2)) - 
      10*nz*pow(nx,3)*(nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   20*tx*pow(nx,3)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
    20*nx*nz*tx*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
    pow(nz,2)*(20*nz*ty*tz*pow(ny,2) - 6*ty*tz*pow(nz,3) + 5*ny*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 
       10*pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 10*pow(nx,2)*
     (-6*nz*ty*tz*pow(ny,2) + 4*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2))) - 
    5*pow(nx,4)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2))),20*tx*pow(nx,3)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 
    20*nx*tx*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 
    pow(nz,3)*(10*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 10*pow(ny,2)*(pow(tx,2) - pow(tz,2))) - 
    10*nz*pow(nx,2)*(6*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    5*pow(nx,4)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2))),20*tx*pow(nx,2)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
    2*nz*tx*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
    5*nx*(-8*nz*ty*tz*pow(ny,3) + 16*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       6*pow(ny,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    10*pow(nx,3)*(-4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2))),
   4*(5*tx*pow(nx,2)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      tx*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      5*nx*nz*(-3*nz*ty*tz*pow(ny,2) + ty*tz*pow(nz,3) - ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 
      5*pow(nx,3)*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2)))),
   10*ty*(nx*tx - nz*tz)*pow(ny,4) + 20*nz*ty*pow(ny,2)*(-3*nx*nz*tx - 3*tz*pow(nx,2) + 2*tz*pow(nz,2)) + 
    2*ty*(5*nx*nz*tx + 10*tz*pow(nx,2) - 3*tz*pow(nz,2))*pow(nz,3) + pow(ny,5)*(pow(tx,2) - pow(tz,2)) - 
    10*pow(ny,3)*(4*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    5*ny*pow(nz,2)*(8*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 6*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   20*ty*pow(ny,3)*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + 20*ny*ty*pow(nz,2)*(-2*nx*nz*tx - 3*tz*pow(nx,2) + tz*pow(nz,2)) + 
    5*pow(ny,4)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
    10*nz*pow(ny,2)*(6*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 3*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    pow(nz,3)*(10*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 10*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   -20*nz*ty*(nz*tx + 2*nx*tz)*pow(ny,3) + 2*tx*ty*pow(ny,5) + 10*ny*ty*(nz*tx + 4*nx*tz)*pow(nz,3) + 
    pow(nz,4)*(-2*nz*tx*tz + 5*nx*(pow(ty,2) - pow(tz,2))) + 10*pow(ny,2)*pow(nz,2)*(2*nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 
    5*pow(ny,4)*(2*nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2))),2*
    (5*ty*(nz*tx + nx*tz)*pow(ny,4) + tx*tz*pow(ny,5) - 10*ty*(nz*tx + 3*nx*tz)*pow(ny,2)*pow(nz,2) + ty*(nz*tx + 5*nx*tz)*pow(nz,4) + 
      5*ny*pow(nz,3)*(nz*tx*tz - 2*nx*pow(ty,2) + 2*nx*pow(tz,2)) - 10*nz*pow(ny,3)*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   -10*nz*ty*tz*pow(ny,4) + 20*ty*tz*pow(ny,2)*pow(nz,3) - 2*ty*tz*pow(nz,5) + pow(ny,5)*(pow(ty,2) - pow(tz,2)) - 
    10*pow(ny,3)*pow(nz,2)*(pow(ty,2) - pow(tz,2)) + 5*ny*pow(nz,4)*(pow(ty,2) - pow(tz,2)),
   2*ty*tz*pow(ny,5) - 20*ty*tz*pow(ny,3)*pow(nz,2) + 10*ny*ty*tz*pow(nz,4) + 5*nz*pow(ny,4)*(pow(ty,2) - pow(tz,2)) - 
    10*pow(ny,2)*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + pow(nz,5)*(pow(ty,2) - pow(tz,2)),
   -5*nz*(rz*tx + rx*tz)*pow(nx,4) + (rx*tx - rz*tz)*pow(nx,5) - 10*(rx*tx - rz*tz)*pow(nx,3)*pow(nz,2) + 10*(rz*tx + rx*tz)*pow(nx,2)*pow(nz,3) + 
    5*nx*(rx*tx - rz*tz)*pow(nz,4) - (rz*tx + rx*tz)*pow(nz,5),
   -10*nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz)*pow(nx,3) - 5*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,4) + 
    (ry*tx + rx*ty)*pow(nx,5) + 10*(-3*ny*rx*tx + nz*rz*ty + nz*ry*tz + 3*ny*rz*tz)*pow(nx,2)*pow(nz,2) + 
    5*nx*(nz*ry*tx + 4*ny*rz*tx + nz*rx*ty + 4*ny*rx*tz)*pow(nz,3) - (-5*ny*rx*tx + nz*rz*ty + nz*ry*tz + 5*ny*rz*tz)*pow(nz,4),
   5*nz*(rx*tx - rz*tz)*pow(nx,4) + (rz*tx + rx*tz)*pow(nx,5) - 10*(rz*tx + rx*tz)*pow(nx,3)*pow(nz,2) - 10*(rx*tx - rz*tz)*pow(nx,2)*pow(nz,3) + 
    5*nx*(rz*tx + rx*tz)*pow(nz,4) + (rx*tx - rz*tz)*pow(nz,5),
   5*(ny*(ry*tx + rx*ty) - nz*(rz*tx + rx*tz))*pow(nx,4) + (ry*ty - rz*tz)*pow(nx,5) - 
    10*nz*pow(nx,2)*(3*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(nz,2)) + 
    5*nx*pow(nz,2)*(4*ny*nz*(rz*ty + ry*tz) + (-6*rx*tx + 6*rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2)) + 
    10*pow(nx,3)*(-2*ny*nz*(rz*ty + ry*tz) + (rx*tx - rz*tz)*pow(ny,2) - (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    (5*ny*nz*(ry*tx + rx*ty) + 10*(rz*tx + rx*tz)*pow(ny,2) - 3*(rz*tx + rx*tz)*pow(nz,2))*pow(nz,3),
   -10*nz*(-2*ny*rx*tx + nz*rz*ty + nz*ry*tz + 2*ny*rz*tz)*pow(nx,3) + 5*(nz*ry*tx + ny*rz*tx + nz*rx*ty + ny*rx*tz)*pow(nx,4) + 
    (rz*ty + ry*tz)*pow(nx,5) - 10*(nz*ry*tx + 3*ny*rz*tx + nz*rx*ty + 3*ny*rx*tz)*pow(nx,2)*pow(nz,2) + 
    5*nx*(-4*ny*rx*tx + nz*rz*ty + nz*ry*tz + 4*ny*rz*tz)*pow(nz,3) + (nz*ry*tx + 5*ny*rz*tx + nz*rx*ty + 5*ny*rx*tz)*pow(nz,4),
   -5*(-(ny*ry*ty) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,4) + 
    10*pow(nx,3)*(-2*ny*nz*(rz*tx + rx*tz) + (ry*tx + rx*ty)*pow(ny,2) - (ry*tx + rx*ty)*pow(nz,2)) + 
    10*nx*nz*(-3*nz*(ry*tx + rx*ty)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(ny,3) + 4*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
    pow(nz,2)*(10*nz*(rz*ty + ry*tz)*pow(ny,2) - 10*(rx*tx - rz*tz)*pow(ny,3) + 5*ny*(2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2) - 
       3*(rz*ty + ry*tz)*pow(nz,3)) + 10*pow(nx,2)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (rx*tx - rz*tz)*pow(ny,3) - 
       3*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + 2*(rz*ty + ry*tz)*pow(nz,3)),
   5*(nz*ry*ty + ny*rz*ty + ny*ry*tz - nz*rz*tz)*pow(nx,4) + 
    10*pow(nx,3)*(2*ny*nz*(ry*tx + rx*ty) + (rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    10*nx*pow(nz,2)*(2*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    10*nz*pow(nx,2)*(3*ny*nz*(rz*ty + ry*tz) + (-3*rx*tx + 3*rz*tz)*pow(ny,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    (5*ny*nz*(rz*ty + ry*tz) - 10*(rx*tx - rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3),
   10*pow(nx,3)*(-2*ny*nz*(rz*ty + ry*tz) + (ry*ty - rz*tz)*pow(ny,2) + (-(ry*ty) + rz*tz)*pow(nz,2)) + 
    10*pow(nx,2)*(-3*nz*(rz*tx + rx*tz)*pow(ny,2) + (ry*tx + rx*ty)*pow(ny,3) - 3*ny*(ry*tx + rx*ty)*pow(nz,2) + (rz*tx + rx*tz)*pow(nz,3)) - 
    nz*(10*nz*(ry*tx + rx*ty)*pow(ny,3) + 5*(rz*tx + rx*tz)*pow(ny,4) - 20*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 10*ny*(ry*tx + rx*ty)*pow(nz,3) + 
       3*(rz*tx + rx*tz)*pow(nz,4)) + 5*nx*(-4*nz*(rz*ty + ry*tz)*pow(ny,3) + (rx*tx - rz*tz)*pow(ny,4) - 
       6*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,2)*pow(nz,2) + 8*ny*(rz*ty + ry*tz)*pow(nz,3) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4)),
   2*(5*pow(nx,3)*(2*ny*nz*(ry*ty - rz*tz) + (rz*ty + ry*tz)*pow(ny,2) - (rz*ty + ry*tz)*pow(nz,2)) + 
      5*pow(nx,2)*(3*nz*(ry*tx + rx*ty)*pow(ny,2) + (rz*tx + rx*tz)*pow(ny,3) - 3*ny*(rz*tx + rx*tz)*pow(nz,2) - (ry*tx + rx*ty)*pow(nz,3)) + 
      pow(nz,2)*(-5*nz*(ry*tx + rx*ty)*pow(ny,2) - 5*(rz*tx + rx*tz)*pow(ny,3) + 5*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
      5*nx*nz*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + 2*(rx*tx - rz*tz)*pow(ny,3) - 2*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3)))
    ,5*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz))*pow(ny,4) + (rx*tx - rz*tz)*pow(ny,5) - 
    10*nz*pow(ny,2)*(3*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - 2*(rz*ty + ry*tz)*pow(nz,2)) + 
    5*ny*pow(nz,2)*(4*nx*nz*(rz*tx + rx*tz) + 6*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2)) - 
    10*pow(ny,3)*(2*nx*nz*(rz*tx + rx*tz) + (-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    (5*nx*nz*(ry*tx + rx*ty) + 10*(rz*ty + ry*tz)*pow(nx,2) - 3*(rz*ty + ry*tz)*pow(nz,2))*pow(nz,3),
   5*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz)*pow(ny,4) + 
    10*pow(ny,3)*(2*nx*nz*(ry*tx + rx*ty) + (rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    10*ny*pow(nz,2)*(2*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    10*nz*pow(ny,2)*(3*nx*nz*(rz*tx + rx*tz) + 3*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    (5*nx*nz*(rz*tx + rx*tz) + 10*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3),
   -10*nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz)*pow(ny,3) - 5*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz)*pow(ny,4) + 
    (ry*tx + rx*ty)*pow(ny,5) + 10*(nz*rz*tx - 3*nx*ry*ty + nz*rx*tz + 3*nx*rz*tz)*pow(ny,2)*pow(nz,2) + 
    5*ny*(nz*ry*tx + nz*rx*ty + 4*nx*rz*ty + 4*nx*ry*tz)*pow(nz,3) - (nz*rz*tx - 5*nx*ry*ty + nz*rx*tz + 5*nx*rz*tz)*pow(nz,4),
   -10*nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz)*pow(ny,3) + 5*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz)*pow(ny,4) + 
    (rz*tx + rx*tz)*pow(ny,5) - 10*(nz*ry*tx + nz*rx*ty + 3*nx*rz*ty + 3*nx*ry*tz)*pow(ny,2)*pow(nz,2) + 
    5*ny*(nz*rz*tx - 4*nx*ry*ty + nz*rx*tz + 4*nx*rz*tz)*pow(nz,3) + (nz*ry*tx + nz*rx*ty + 5*nx*rz*ty + 5*nx*ry*tz)*pow(nz,4),
   -5*nz*(rz*ty + ry*tz)*pow(ny,4) + (ry*ty - rz*tz)*pow(ny,5) - 10*(ry*ty - rz*tz)*pow(ny,3)*pow(nz,2) + 10*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,3) + 
    5*ny*(ry*ty - rz*tz)*pow(nz,4) - (rz*ty + ry*tz)*pow(nz,5),
   5*nz*(ry*ty - rz*tz)*pow(ny,4) + (rz*ty + ry*tz)*pow(ny,5) - 10*(rz*ty + ry*tz)*pow(ny,3)*pow(nz,2) - 10*(ry*ty - rz*tz)*pow(ny,2)*pow(nz,3) + 
    5*ny*(rz*ty + ry*tz)*pow(nz,4) + (ry*ty - rz*tz)*pow(nz,5),
   -6*tx*pow(nx,2)*pow(nz,2)*(pow(tx,2) - 3*pow(tz,2)) + tx*pow(nz,4)*(pow(tx,2) - 3*pow(tz,2)) + 4*nz*tz*pow(nx,3)*(-3*pow(tx,2) + pow(tz,2)) - 
    4*nx*tz*pow(nz,3)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,4)*(pow(tx,3) - 3*tx*pow(tz,2)),
   4*tx*pow(nx,3)*(-6*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 3*ty*pow(nx,4)*(pow(tx,2) - pow(tz,2)) + 
    12*nx*tx*pow(nz,2)*(2*nz*ty*tz - ny*pow(tx,2) + 3*ny*pow(tz,2)) + 
    pow(nz,3)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - 4*ny*tz*(-3*pow(tx,2) + pow(tz,2))) - 
    6*nz*pow(nx,2)*(6*ny*tz*pow(tx,2) + 3*nz*ty*(pow(tx,2) - pow(tz,2)) - 2*ny*pow(tz,3)),
   4*nz*tx*pow(nx,3)*(pow(tx,2) - 3*pow(tz,2)) - 4*nx*tx*pow(nz,3)*(pow(tx,2) - 3*pow(tz,2)) + 6*tz*pow(nx,2)*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) - 
    tz*pow(nz,4)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,4)*(3*tz*pow(tx,2) - pow(tz,3)),
   tx*pow(nz,2)*(24*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    6*tx*pow(nx,2)*(-12*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    3*tx*pow(nx,4)*(pow(ty,2) - pow(tz,2)) + 4*pow(nx,3)*(3*ny*ty*(pow(tx,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    12*nx*nz*(tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 3*ny*nz*ty*(-pow(tx,2) + pow(tz,2)) + pow(ny,2)*(-3*tz*pow(tx,2) + pow(tz,3))),
   2*(3*tx*ty*tz*pow(nx,4) + 6*nz*tx*pow(nx,2)*(-3*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 
      tx*pow(nz,3)*(3*nz*ty*tz - 2*ny*pow(tx,2) + 6*ny*pow(tz,2)) - 
      6*nx*pow(nz,2)*(nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2))) + 
      pow(nx,3)*(6*ny*tz*pow(tx,2) + 6*nz*ty*(pow(tx,2) - pow(tz,2)) - 2*ny*pow(tz,3))),
   4*nx*tx*(-18*nz*ty*tz*pow(ny,2) + 12*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
       pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + pow(nx,4)*(pow(ty,3) - 3*ty*pow(tz,2)) + 
    nz*(ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 18*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       12*ny*tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 4*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) - 
    12*tx*pow(nx,3)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2))) + 
    6*pow(nx,2)*(-(ty*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       2*ny*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))),
   -4*nx*nz*tx*(18*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    12*tx*pow(nx,3)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2))) + 
    6*pow(nx,2)*(6*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + 
       pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + pow(nx,4)*(3*tz*pow(ty,2) - pow(tz,3)) + 
    3*pow(nz,2)*(tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 4*ny*nz*ty*(-pow(tx,2) + pow(tz,2)) + 
       pow(ny,2)*(-6*tz*pow(tx,2) + 2*pow(tz,3))),4*ny*ty*(-18*nz*tx*tz*pow(nx,2) + 12*tx*tz*pow(nz,3) - 
       3*nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + pow(ny,4)*(pow(tx,3) - 3*tx*pow(tz,2)) - 
    12*ty*pow(ny,3)*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))) + 
    nz*(tx*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 18*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       12*nx*tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 4*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
    6*pow(ny,2)*(tx*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 2*nx*nz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       3*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),4*(nz*ty*(-9*nz*tx*tz*pow(nx,2) + 3*tx*tz*pow(nz,3) - 
         nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 
      9*ty*pow(ny,2)*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))) - 
      ny*(tx*pow(nz,3)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 9*nz*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
      pow(ny,3)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3))),
   ty*pow(nz,2)*(24*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
    6*ty*pow(ny,2)*(12*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
    3*ty*pow(ny,4)*(pow(tx,2) - pow(tz,2)) + 4*pow(ny,3)*(3*nx*tx*(pow(ty,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    12*ny*nz*(tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 3*nx*nz*tx*(-pow(ty,2) + pow(tz,2)) + pow(nx,2)*(-3*tz*pow(ty,2) + pow(tz,3))),
   -4*ny*nz*ty*(18*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 3*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
    12*ty*pow(ny,3)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) + 
    6*pow(ny,2)*(6*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + 
       pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + pow(ny,4)*(3*tz*pow(tx,2) - pow(tz,3)) + 
    3*pow(nz,2)*(tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 4*nx*nz*tx*(-pow(ty,2) + pow(tz,2)) + 
       pow(nx,2)*(-6*tz*pow(ty,2) + 2*pow(tz,3))),4*ty*pow(ny,3)*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
    3*tx*pow(ny,4)*(pow(ty,2) - pow(tz,2)) + 12*ny*ty*pow(nz,2)*(2*nz*tx*tz - nx*pow(ty,2) + 3*nx*pow(tz,2)) + 
    pow(nz,3)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - 4*nx*tz*(-3*pow(ty,2) + pow(tz,2))) - 
    6*nz*pow(ny,2)*(6*nx*tz*pow(ty,2) + 3*nz*tx*(pow(ty,2) - pow(tz,2)) - 2*nx*pow(tz,3)),
   2*(3*tx*ty*tz*pow(ny,4) + 6*nz*ty*pow(ny,2)*(-3*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
      ty*pow(nz,3)*(3*nz*tx*tz - 2*nx*pow(ty,2) + 6*nx*pow(tz,2)) - 
      6*ny*pow(nz,2)*(nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2))) + 
      pow(ny,3)*(6*nx*tz*pow(ty,2) + 6*nz*tx*(pow(ty,2) - pow(tz,2)) - 2*nx*pow(tz,3))),
   -6*ty*pow(ny,2)*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2)) + ty*pow(nz,4)*(pow(ty,2) - 3*pow(tz,2)) + 4*nz*tz*pow(ny,3)*(-3*pow(ty,2) + pow(tz,2)) - 
    4*ny*tz*pow(nz,3)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,4)*(pow(ty,3) - 3*ty*pow(tz,2)),
   4*nz*ty*pow(ny,3)*(pow(ty,2) - 3*pow(tz,2)) - 4*ny*ty*pow(nz,3)*(pow(ty,2) - 3*pow(tz,2)) + 6*tz*pow(ny,2)*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) - 
    tz*pow(nz,4)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,4)*(3*tz*pow(ty,2) - pow(tz,3)),
   pow(nx,4)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + pow(nz,4)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - 
    4*nz*pow(nx,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 4*nx*pow(nz,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    6*pow(nx,2)*pow(nz,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   pow(nx,4)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    4*pow(nx,3)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    4*nx*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) - 
    6*nz*pow(nx,2)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 2*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
    pow(nz,3)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 4*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(nx,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*pow(nx,2)*pow(nz,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    pow(nz,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 4*nz*pow(nx,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    4*nx*pow(nz,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   6*pow(nx,2)*(-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) - pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + pow(nx,4)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    pow(nz,2)*(8*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ry*tx*ty - 6*rz*tx*tz + 2*rx*pow(tx,2) + rx*pow(ty,2) - 3*rx*pow(tz,2)) + 
       pow(ny,2)*(12*rz*tx*tz - 6*rx*pow(tx,2) + 6*rx*pow(tz,2))) + 
    4*pow(nx,3)*(-(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       ny*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    4*nx*nz*(3*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
    ,2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,4) - 6*nz*pow(nx,2)*
       (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*rz*tx*tz - rx*pow(tx,2) + rx*pow(tz,2))) + 
      pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(4*rz*tx*tz - 2*rx*pow(tx,2) + 2*rx*pow(tz,2))) + 
      2*pow(nx,3)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
      2*nx*pow(nz,2)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 3*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))),
   4*nx*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*ny*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    pow(nx,4)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + 4*pow(nx,3)*
     (-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    6*pow(nx,2)*(-(pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
       2*ny*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)))
      + nz*(pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       4*ny*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   pow(nx,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 4*pow(nx,3)*
     (2*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    4*nx*nz*(6*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) + 
    6*pow(nx,2)*(-(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       pow(ny,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 2*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,2)*(4*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
    ,4*ny*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
    pow(ny,4)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + nz*
     (pow(nz,3)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
       4*nx*pow(nz,2)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 6*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    6*pow(ny,2)*(pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))
       ) + 4*pow(ny,3)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   4*(nz*(-3*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
         nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
      3*pow(ny,2)*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
         nx*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
      pow(ny,3)*(nz*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + nx*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
      ny*(pow(nz,3)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
         3*nx*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
         3*nz*pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,3)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))))),
   pow(ny,4)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) - 
    4*pow(ny,3)*(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       nx*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    pow(nz,2)*(8*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 
       6*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    4*ny*nz*(3*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 6*pow(ny,2)*(4*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 4*pow(ny,3)*
     (2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,2)*(4*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 4*ny*nz*(6*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))) - 
    6*pow(ny,2)*(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,2)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,4)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) - 
    4*pow(ny,3)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    4*ny*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 3*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
    6*nz*pow(ny,2)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 2*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
    pow(nz,3)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 4*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,4) - 6*nz*pow(ny,2)*
       (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
      pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 2*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
      2*pow(ny,3)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) - 
      2*ny*pow(nz,2)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 3*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))),
   pow(ny,4)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + pow(nz,4)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 
    4*nz*pow(ny,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 4*ny*pow(nz,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    6*pow(ny,2)*pow(nz,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),
   pow(ny,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 6*pow(ny,2)*pow(nz,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    pow(nz,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 4*nz*pow(ny,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + 
    4*ny*pow(nz,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),
   -12*nz*tx*tz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 4*tx*tz*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + 
    pow(nx,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 3*nx*pow(nz,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   -12*nx*nz*tx*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*ny*tz*(pow(tx,2) - pow(tz,2))) + 4*tx*ty*pow(nx,3)*(pow(tx,2) - 3*pow(tz,2)) + 
    pow(nz,2)*(-4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) - 3*ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    3*pow(nx,2)*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   4*tx*tz*pow(nx,3)*(pow(tx,2) - pow(tz,2)) - 12*nx*tx*tz*pow(nz,2)*(pow(tx,2) - pow(tz,2)) + 
    3*nz*pow(nx,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - pow(nz,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   12*tx*pow(nx,2)*(ny*ty*(pow(tx,2) - 3*pow(tz,2)) + nz*tz*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    4*nz*tx*(3*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       tz*pow(nz,2)*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
    2*pow(nx,3)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    3*nx*(8*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
       pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   2*(6*tx*pow(nx,2)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + ny*tz*(pow(tx,2) - pow(tz,2))) - 
      2*tx*pow(nz,2)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*ny*tz*(pow(tx,2) - pow(tz,2))) + pow(nx,3)*(6*ty*tz*pow(tx,2) - 2*ty*pow(tz,3)) + 
      3*nx*nz*(2*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*ty*(-3*nx*tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
       tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 3*nz*tz*pow(nx,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
    12*ty*pow(ny,2)*(nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    3*ny*(8*nx*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -12*nx*tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 4*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 
    12*ny*ty*(2*nx*nz*tx*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
    6*nz*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    pow(nz,3)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    3*pow(ny,2)*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   -12*ny*ty*(tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
       2*nx*nz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 4*tx*ty*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2)) + 
    4*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 12*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
    pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*nx*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    6*pow(ny,2)*(-2*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + nx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(tx*tz*pow(ny,3)*(pow(tx,2) - pow(tz,2)) - ty*(tx*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 
         3*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 3*ty*pow(ny,2)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)) + 
      3*ny*(-(tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) - tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         nx*nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   12*ty*pow(ny,2)*(nx*tx*(pow(ty,2) - 3*pow(tz,2)) + nz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) - 
    4*nz*ty*(3*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       tz*pow(nz,2)*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
    2*pow(ny,3)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*ny*(-8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -4*ty*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2)) + 12*ny*ty*
     (2*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
    6*pow(ny,2)*(-2*nx*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nz*(12*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -12*ny*nz*ty*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 2*nx*tz*(pow(ty,2) - pow(tz,2))) + 4*tx*ty*pow(ny,3)*(pow(ty,2) - 3*pow(tz,2)) + 
    pow(nz,2)*(-4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 3*nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    3*pow(ny,2)*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   2*(6*ty*pow(ny,2)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + nx*tz*(pow(ty,2) - pow(tz,2))) - 
      2*ty*pow(nz,2)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*nx*tz*(pow(ty,2) - pow(tz,2))) + tx*pow(ny,3)*(6*tz*pow(ty,2) - 2*pow(tz,3)) + 
      3*ny*nz*(2*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   -12*nz*ty*tz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + 4*ty*tz*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + 
    pow(ny,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*ny*pow(nz,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*ty*tz*pow(ny,3)*(pow(ty,2) - pow(tz,2)) - 12*ny*ty*tz*pow(nz,2)*(pow(ty,2) - pow(tz,2)) + 
    3*nz*pow(ny,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - pow(nz,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   pow(nx,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    3*nx*pow(nz,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nz,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    3*nz*pow(nx,2)*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   pow(nx,3)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    3*nx*nz*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    3*pow(nx,2)*(ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) + 
    pow(nz,2)*(nz*(6*rx*tx*ty*tz + 3*rz*ty*pow(tx,2) + 3*ry*tz*pow(tx,2) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) - 
       3*ny*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3))),
   3*nz*pow(nx,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    pow(nz,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nx,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    3*nx*pow(nz,2)*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   pow(nx,3)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    3*nx*(pow(nz,2)*(-(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
          3*ry*ty*(-pow(tx,2) + pow(tz,2))) - 2*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    3*pow(nx,2)*(-(nz*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2)))) + 
       ny*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    nz*(pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       3*ny*nz*(6*rz*tx*ty*tz - 3*rx*ty*pow(tx,2) - ry*pow(tx,3) + 3*ry*tx*pow(tz,2) + 3*rx*ty*pow(tz,2)) + 
       pow(ny,2)*(-9*rx*tz*pow(tx,2) - 3*rz*pow(tx,3) + 9*rz*tx*pow(tz,2) + 3*rx*pow(tz,3))),
   pow(nx,3)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    3*pow(nx,2)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    pow(nz,2)*(nz*(6*rz*tx*ty*tz - 3*rx*ty*pow(tx,2) - ry*pow(tx,3) + 3*ry*tx*pow(tz,2) + 3*rx*ty*pow(tz,2)) + 
       3*ny*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3))) + 
    3*nx*nz*(2*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
    3*nz*pow(nx,2)*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) - 
    3*nx*pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    pow(nx,3)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    pow(ny,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    3*ny*(2*nx*nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + 3*pow(ny,2)*(nx*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   -3*nx*pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    pow(nz,3)*(-(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
       3*ry*ty*(-pow(tx,2) + pow(tz,2))) + 3*nz*pow(nx,2)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
       rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + pow(nx,3)*
     (3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    3*pow(ny,2)*(nz*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nx*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    3*ny*(pow(nx,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       2*nx*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
    3*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(ny,3)*(ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))) - 
    3*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
       3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))) - 
    3*pow(ny,2)*(nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       nx*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2)))) - 
    3*ny*(pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       2*nx*nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))),
   -3*nx*pow(nz,2)*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) - 
    pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    3*nz*pow(nx,2)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) - 
    3*ny*(pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       2*nx*nz*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nx,2)*(-6*ry*tx*ty*tz + rx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)))) + 
    3*pow(ny,2)*(nx*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    pow(ny,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    pow(nx,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   pow(ny,3)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    3*pow(ny,2)*(nx*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) - 
    3*ny*(2*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     + nz*(pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       3*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,3)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    3*pow(ny,2)*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       nx*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)))) - 
    nz*(3*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       3*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     - 3*ny*(2*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,3)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
    3*ny*nz*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    pow(nz,2)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       3*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))) - 
    3*pow(ny,2)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,3)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    3*pow(ny,2)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    pow(nz,2)*(nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*nx*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))) - 
    3*ny*nz*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       2*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    3*ny*pow(nz,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(nz,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    3*nz*pow(ny,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   3*nz*pow(ny,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    pow(nz,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(ny,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    3*ny*pow(nz,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   -2*nx*nz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - tx*pow(nz,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(nx,2)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   5*ty*pow(nx,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    nz*(2*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*nx*tx*(20*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))),
   -(tz*pow(nz,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 2*nx*nz*tx*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(nx,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)),
   -(tx*pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)))) + 
    10*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    2*nx*nz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(ny,2)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
    10*ny*ty*(4*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + nx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   2*(5*ty*(2*tx*tz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 2*tx*tz*pow(nz,2)*(-pow(tx,2) + pow(tz,2)) + 
         nx*nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)))),
   5*ty*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    2*ny*(nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       10*nx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    5*ty*(8*nx*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   tz*pow(nz,2)*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
    20*nx*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    10*ny*ty*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(ny,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)) + 
    2*pow(nx,2)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)),
   5*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    10*tx*pow(ny,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    2*nx*nz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
    5*tx*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    20*ny*ty*(-2*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + nx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(5*tx*ty*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + ny*(nx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         5*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      5*ty*(-(tx*tz*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + tx*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         nx*nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   10*ty*pow(ny,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    2*ny*(nz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*nx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    ty*(40*nx*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
       pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))),
   tz*pow(nz,2)*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
    10*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    20*ny*ty*(2*nx*tx*tz*(pow(ty,2) - pow(tz,2)) + nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) + 
    2*pow(ny,2)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)),
   5*tx*pow(ny,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    nz*(2*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*ny*ty*(20*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   2*(10*tx*ty*tz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + ny*(nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nz*ty*(10*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   -2*ny*nz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ty*pow(nz,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(ny,2)*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   -(tz*pow(nz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 2*ny*nz*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(ny,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)),
   pow(nx,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    2*nx*nz*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   pow(nx,2)*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*nx*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       2*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   2*nx*nz*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   pow(ny,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(nx,2)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,2)*(4*tx*(-(ry*ty*(pow(tx,2) - 3*pow(tz,2))) + rz*tz*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) - 
       rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    2*nx*nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
       rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    2*ny*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nx*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   2*(2*pow(nx,2)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
      2*pow(nz,2)*(-3*rx*ty*tz*pow(tx,2) - ry*tz*pow(tx,3) - rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + ry*tx*pow(tz,3) + rx*ty*pow(tz,3)) + 
      nx*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(nz*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   -8*nx*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
       rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + pow(ny,2)*
     (4*ty*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(nx,2)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
       ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    2*ny*(2*nx*(-2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   pow(ny,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    4*nx*nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(nx,2)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
       rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
       rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    2*ny*(4*nx*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(ny,2)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
       rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    2*nx*nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    4*ny*(2*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       nx*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   4*(-(pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
           rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
      pow(nx,2)*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
      pow(ny,2)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))) + 
      nx*nz*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
            rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
            rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   -8*nx*nz*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(ny,2)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,2)*(4*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    2*ny*(-(nx*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(ny,2)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
       rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    4*ny*(2*nx*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nz*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   pow(ny,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*ny*(-4*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*(2*pow(ny,2)*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
      nz*(-2*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      ny*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    2*ny*nz*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   2*ny*nz*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(ny,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   -2*nz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + nx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   6*ty*(-(nz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + nx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    ny*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   2*nx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   6*ny*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
    2*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
    3*nx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   6*nx*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 6*nz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    ny*(6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5)),
   2*ty*(nz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       10*nx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    3*ny*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*tz*(2*nx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
       3*ny*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    3*nz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*nz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
    20*ny*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    3*nx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(5*nz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      tz*(ny*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
         nx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   3*ny*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*ty*(nz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       3*nx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   3*nz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*tz*(3*nx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*ny*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -6*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ny*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    nx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   6*ny*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*nz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    2*nx*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)),
   -2*nz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + ny*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*ny*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + nz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   nx*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    nz*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   ny*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    nx*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   nz*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    nx*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   nx*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    ny*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))
    ,ny*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    nx*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   ny*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    10*nx*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   nz*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    ny*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    2*nx*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   nx*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    10*ny*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(ny*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      nx*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      5*nz*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   ny*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    nz*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))\
     + nx*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   nz*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nx*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    2*ny*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   nx*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    ny*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    ny*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   ny*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    nz*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   nz*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    ny*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6),
   7*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   7*tz*pow(tx,6) - 35*pow(tx,4)*pow(tz,3) + 21*pow(tx,2)*pow(tz,5) - pow(tz,7),
   7*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   14*tx*ty*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)),
   7*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   35*pow(tx,4)*(3*tz*pow(ty,2) - pow(tz,3)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,3) - pow(tz,5)) + 21*pow(ty,2)*pow(tz,5) - 3*pow(tz,7),
   7*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   28*tx*ty*tz*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),
   7*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)),
   -35*pow(ty,4)*pow(tz,3) + 42*pow(ty,2)*pow(tz,5) + 21*pow(tx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) - 3*pow(tz,7),
   7*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),14*tx*ty*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)),
   pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6),
   7*tz*pow(ty,6) - 35*pow(ty,4)*pow(tz,3) + 21*pow(ty,2)*pow(tz,5) - pow(tz,7),
   -2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
    2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
    3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   6*rx*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 6*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    ry*(6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5)),
   2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
       3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
    20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
         rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   6*ry*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    2*rx*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)),
   -2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6));


    

    degree = 8;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,8);

    tensor_project[degree].P << pow(nx,8) - 28*pow(nx,6)*pow(nz,2) + 70*pow(nx,4)*pow(nz,4) - 28*pow(nx,2)*pow(nz,6) + pow(nz,8),
   8*nx*ny*(pow(nx,6) - 21*pow(nx,4)*pow(nz,2) + 35*pow(nx,2)*pow(nz,4) - 7*pow(nz,6)),
   8*nx*nz*(pow(nx,6) - 7*pow(nx,4)*pow(nz,2) + 7*pow(nx,2)*pow(nz,4) - pow(nz,6)),
   4*(7*pow(nx,6)*(pow(ny,2) - pow(nz,2)) - 35*pow(nx,4)*(3*pow(ny,2)*pow(nz,2) - pow(nz,4)) + 21*pow(nx,2)*(5*pow(ny,2)*pow(nz,4) - pow(nz,6)) - 
      7*pow(ny,2)*pow(nz,6) + pow(nz,8)),-8*ny*nz*(-7*pow(nx,6) + 35*pow(nx,4)*pow(nz,2) - 21*pow(nx,2)*pow(nz,4) + pow(nz,6)),
   56*nx*ny*(pow(nx,4)*(pow(ny,2) - 3*pow(nz,2)) + 10*pow(nx,2)*pow(nz,2)*(-pow(ny,2) + pow(nz,2)) + 5*pow(ny,2)*pow(nz,4) - 3*pow(nz,6)),
   8*nx*nz*(7*pow(nx,4)*(3*pow(ny,2) - pow(nz,2)) + 14*pow(nx,2)*pow(nz,2)*(-5*pow(ny,2) + pow(nz,2)) + 21*pow(ny,2)*pow(nz,4) - 3*pow(nz,6)),
   70*pow(ny,4)*pow(nz,4) + 70*pow(nx,4)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 84*pow(ny,2)*pow(nz,6) - 
    84*pow(nx,2)*(5*pow(ny,4)*pow(nz,2) - 10*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 6*pow(nz,8),
   8*ny*nz*(35*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 7*pow(ny,2)*pow(nz,4) + pow(nx,2)*(-70*pow(ny,2)*pow(nz,2) + 42*pow(nz,4)) - 3*pow(nz,6)),
   56*nx*ny*(-3*pow(ny,4)*pow(nz,2) + 10*pow(ny,2)*pow(nz,4) + pow(nx,2)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) - 3*pow(nz,6)),
   8*nx*nz*(-35*pow(ny,4)*pow(nz,2) + 42*pow(ny,2)*pow(nz,4) + 7*pow(nx,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 3*pow(nz,6)),
   4*(-7*pow(ny,6)*pow(nz,2) + 35*pow(ny,4)*pow(nz,4) + 7*pow(nx,2)*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) - 
      21*pow(ny,2)*pow(nz,6) + pow(nz,8)),8*ny*nz*(-7*pow(ny,4)*pow(nz,2) + 14*pow(ny,2)*pow(nz,4) + 
      7*pow(nx,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) - 3*pow(nz,6)),
   8*nx*ny*(pow(ny,6) - 21*pow(ny,4)*pow(nz,2) + 35*pow(ny,2)*pow(nz,4) - 7*pow(nz,6)),
   -8*nx*nz*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6)),
   pow(ny,8) - 28*pow(ny,6)*pow(nz,2) + 70*pow(ny,4)*pow(nz,4) - 28*pow(ny,2)*pow(nz,6) + pow(nz,8),
   8*ny*nz*(pow(ny,6) - 7*pow(ny,4)*pow(nz,2) + 7*pow(ny,2)*pow(nz,4) - pow(nz,6)),
   -7*nz*tz*pow(nx,6) + tx*pow(nx,7) - 21*tx*pow(nx,5)*pow(nz,2) + 35*tz*pow(nx,4)*pow(nz,3) + 35*tx*pow(nx,3)*pow(nz,4) - 
    21*tz*pow(nx,2)*pow(nz,5) - 7*nx*tx*pow(nz,6) + tz*pow(nz,7),
   -21*nz*(nz*ty + 2*ny*tz)*pow(nx,5) + 7*ny*tx*pow(nx,6) + ty*pow(nx,7) - 105*ny*tx*pow(nx,4)*pow(nz,2) + 
    35*(nz*ty + 4*ny*tz)*pow(nx,3)*pow(nz,3) + 105*ny*tx*pow(nx,2)*pow(nz,4) - 7*nx*(nz*ty + 6*ny*tz)*pow(nz,5) - 7*ny*tx*pow(nz,6),
   7*nz*tx*pow(nx,6) + tz*pow(nx,7) - 21*tz*pow(nx,5)*pow(nz,2) - 35*tx*pow(nx,4)*pow(nz,3) + 35*tz*pow(nx,3)*pow(nz,4) + 
    21*tx*pow(nx,2)*pow(nz,5) - 7*nx*tz*pow(nz,6) - tx*pow(nz,7),
   7*(ny*ty - nz*tz)*pow(nx,6) + 21*tx*pow(nx,5)*(pow(ny,2) - pow(nz,2)) + 70*tx*pow(nx,3)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    35*nz*pow(nx,4)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 21*pow(nx,2)*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) - 
    21*nx*tx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4) + (-7*ny*nz*ty - 21*tz*pow(ny,2) + 4*tz*pow(nz,2))*pow(nz,5),
   42*ny*nz*tx*pow(nx,5) + 7*(nz*ty + ny*tz)*pow(nx,6) - 35*(nz*ty + 3*ny*tz)*pow(nx,4)*pow(nz,2) - 140*ny*tx*pow(nx,3)*pow(nz,3) + 
    21*(nz*ty + 5*ny*tz)*pow(nx,2)*pow(nz,4) + 42*nx*ny*tx*pow(nz,5) - (nz*ty + 7*ny*tz)*pow(nz,6),
   7*(5*ny*tx*pow(nx,4)*(pow(ny,2) - 3*pow(nz,2)) - 30*ny*tx*pow(nx,2)*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
      3*pow(nx,5)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
      nx*pow(nz,3)*(15*nz*ty*pow(ny,2) + 20*tz*pow(ny,3) - 18*ny*tz*pow(nz,2) - 3*ty*pow(nz,3)) + 
      10*nz*pow(nx,3)*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + tx*(5*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,4)),
   -35*nz*tx*pow(nx,4)*(-3*pow(ny,2) + pow(nz,2)) + 21*pow(nx,5)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 
    70*pow(nx,3)*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 42*tx*pow(nx,2)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
    21*nx*(2*ny*nz*ty + 5*tz*pow(ny,2) - tz*pow(nz,2))*pow(nz,4) - 3*tx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,5),
   35*pow(nx,4)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
    21*nx*tx*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 35*tx*pow(nx,3)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    21*nz*pow(nx,2)*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
    pow(nz,3)*(35*nz*ty*pow(ny,3) + 35*tz*pow(ny,4) - 63*tz*pow(ny,2)*pow(nz,2) - 21*ny*ty*pow(nz,3) + 6*tz*pow(nz,4)),
   140*ny*nz*tx*pow(nx,3)*(pow(ny,2) - pow(nz,2)) + 28*nx*tx*(-5*pow(ny,3) + 3*ny*pow(nz,2))*pow(nz,3) + 
    35*pow(nx,4)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
    42*pow(nx,2)*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
    (21*nz*ty*pow(ny,2) + 35*tz*pow(ny,3) - 21*ny*tz*pow(nz,2) - 3*ty*pow(nz,3))*pow(nz,4),
   7*(ny*tx*pow(nz,2)*(-3*pow(ny,4) + 10*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 3*ny*tx*pow(nx,2)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
      5*pow(nx,3)*(-4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4)) - 
      nx*nz*(15*nz*ty*pow(ny,4) + 6*tz*pow(ny,5) - 40*tz*pow(ny,3)*pow(nz,2) - 30*ty*pow(ny,2)*pow(nz,3) + 18*ny*tz*pow(nz,4) + 3*ty*pow(nz,5))),
   tx*pow(nz,3)*(-35*pow(ny,4) + 42*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 21*nz*tx*pow(nx,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    35*pow(nx,3)*(4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4)) - 
    7*nx*pow(nz,2)*(20*nz*ty*pow(ny,3) + 15*tz*pow(ny,4) - 30*tz*pow(ny,2)*pow(nz,2) - 12*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)),
   21*pow(nx,2)*(-5*nz*tz*pow(ny,4) + ty*pow(ny,5) - 10*ty*pow(ny,3)*pow(nz,2) + 10*tz*pow(ny,2)*pow(nz,3) + 5*ny*ty*pow(nz,4) - tz*pow(nz,5)) + 
    7*nx*tx*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) + 
    nz*(-21*nz*ty*pow(ny,5) - 7*tz*pow(ny,6) + 70*tz*pow(ny,4)*pow(nz,2) + 70*ty*pow(ny,3)*pow(nz,3) - 63*tz*pow(ny,2)*pow(nz,4) - 
       21*ny*ty*pow(nz,5) + 4*tz*pow(nz,6)),14*nx*ny*nz*tx*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
    21*pow(nx,2)*(5*nz*ty*pow(ny,4) + tz*pow(ny,5) - 10*tz*pow(ny,3)*pow(nz,2) - 10*ty*pow(ny,2)*pow(nz,3) + 5*ny*tz*pow(nz,4) + ty*pow(nz,5)) - 
    pow(nz,2)*(35*nz*ty*pow(ny,4) + 21*tz*pow(ny,5) - 70*tz*pow(ny,3)*pow(nz,2) - 42*ty*pow(ny,2)*pow(nz,3) + 21*ny*tz*pow(nz,4) + 3*ty*pow(nz,5)),
   -21*nz*(nz*tx + 2*nx*tz)*pow(ny,5) + 7*nx*ty*pow(ny,6) + tx*pow(ny,7) - 105*nx*ty*pow(ny,4)*pow(nz,2) + 
    35*(nz*tx + 4*nx*tz)*pow(ny,3)*pow(nz,3) + 105*nx*ty*pow(ny,2)*pow(nz,4) - 7*ny*(nz*tx + 6*nx*tz)*pow(nz,5) - 7*nx*ty*pow(nz,6),
   42*nx*nz*ty*pow(ny,5) + 7*(nz*tx + nx*tz)*pow(ny,6) - 35*(nz*tx + 3*nx*tz)*pow(ny,4)*pow(nz,2) - 140*nx*ty*pow(ny,3)*pow(nz,3) + 
    21*(nz*tx + 5*nx*tz)*pow(ny,2)*pow(nz,4) + 42*nx*ny*ty*pow(nz,5) - (nz*tx + 7*nx*tz)*pow(nz,6),
   -7*nz*tz*pow(ny,6) + ty*pow(ny,7) - 21*ty*pow(ny,5)*pow(nz,2) + 35*tz*pow(ny,4)*pow(nz,3) + 35*ty*pow(ny,3)*pow(nz,4) - 
    21*tz*pow(ny,2)*pow(nz,5) - 7*ny*ty*pow(nz,6) + tz*pow(nz,7),
   7*nz*ty*pow(ny,6) + tz*pow(ny,7) - 21*tz*pow(ny,5)*pow(nz,2) - 35*ty*pow(ny,4)*pow(nz,3) + 35*tz*pow(ny,3)*pow(nz,4) + 
    21*ty*pow(ny,2)*pow(nz,5) - 7*ny*tz*pow(nz,6) - ty*pow(nz,7),
   -7*nz*rz*pow(nx,6) + rx*pow(nx,7) - 21*rx*pow(nx,5)*pow(nz,2) + 35*rz*pow(nx,4)*pow(nz,3) + 35*rx*pow(nx,3)*pow(nz,4) - 
    21*rz*pow(nx,2)*pow(nz,5) - 7*nx*rx*pow(nz,6) + rz*pow(nz,7),
   -21*nz*(nz*ry + 2*ny*rz)*pow(nx,5) + 7*ny*rx*pow(nx,6) + ry*pow(nx,7) - 105*ny*rx*pow(nx,4)*pow(nz,2) + 
    35*(nz*ry + 4*ny*rz)*pow(nx,3)*pow(nz,3) + 105*ny*rx*pow(nx,2)*pow(nz,4) - 7*nx*(nz*ry + 6*ny*rz)*pow(nz,5) - 7*ny*rx*pow(nz,6),
   7*nz*rx*pow(nx,6) + rz*pow(nx,7) - 21*rz*pow(nx,5)*pow(nz,2) - 35*rx*pow(nx,4)*pow(nz,3) + 35*rz*pow(nx,3)*pow(nz,4) + 
    21*rx*pow(nx,2)*pow(nz,5) - 7*nx*rz*pow(nz,6) - rx*pow(nz,7),
   7*(ny*ry - nz*rz)*pow(nx,6) + 21*rx*pow(nx,5)*(pow(ny,2) - pow(nz,2)) + 70*rx*pow(nx,3)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
    35*nz*pow(nx,4)*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2)) + 21*pow(nx,2)*(5*ny*nz*ry + 10*rz*pow(ny,2) - 3*rz*pow(nz,2))*pow(nz,3) - 
    21*nx*rx*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4) + (-7*ny*nz*ry - 21*rz*pow(ny,2) + 4*rz*pow(nz,2))*pow(nz,5),
   42*ny*nz*rx*pow(nx,5) + 7*(nz*ry + ny*rz)*pow(nx,6) - 35*(nz*ry + 3*ny*rz)*pow(nx,4)*pow(nz,2) - 140*ny*rx*pow(nx,3)*pow(nz,3) + 
    21*(nz*ry + 5*ny*rz)*pow(nx,2)*pow(nz,4) + 42*nx*ny*rx*pow(nz,5) - (nz*ry + 7*ny*rz)*pow(nz,6),
   7*(5*ny*rx*pow(nx,4)*(pow(ny,2) - 3*pow(nz,2)) - 30*ny*rx*pow(nx,2)*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
      3*pow(nx,5)*(-2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2)) + 
      nx*pow(nz,3)*(15*nz*ry*pow(ny,2) + 20*rz*pow(ny,3) - 18*ny*rz*pow(nz,2) - 3*ry*pow(nz,3)) + 
      10*nz*pow(nx,3)*(-3*nz*ry*pow(ny,2) - 2*rz*pow(ny,3) + 4*ny*rz*pow(nz,2) + ry*pow(nz,3)) + rx*(5*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,4)),
   -35*nz*rx*pow(nx,4)*(-3*pow(ny,2) + pow(nz,2)) + 21*pow(nx,5)*(2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2)) + 
    70*pow(nx,3)*pow(nz,2)*(-2*ny*nz*ry - 3*rz*pow(ny,2) + rz*pow(nz,2)) + 42*rx*pow(nx,2)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
    21*nx*(2*ny*nz*ry + 5*rz*pow(ny,2) - rz*pow(nz,2))*pow(nz,4) - 3*rx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,5),
   35*pow(nx,4)*(-3*nz*rz*pow(ny,2) + ry*pow(ny,3) - 3*ny*ry*pow(nz,2) + rz*pow(nz,3)) - 
    21*nx*rx*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 35*rx*pow(nx,3)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    21*nz*pow(nx,2)*(10*nz*ry*pow(ny,3) + 5*rz*pow(ny,4) - 20*rz*pow(ny,2)*pow(nz,2) - 10*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)) + 
    pow(nz,3)*(35*nz*ry*pow(ny,3) + 35*rz*pow(ny,4) - 63*rz*pow(ny,2)*pow(nz,2) - 21*ny*ry*pow(nz,3) + 6*rz*pow(nz,4)),
   140*ny*nz*rx*pow(nx,3)*(pow(ny,2) - pow(nz,2)) + 28*nx*rx*(-5*pow(ny,3) + 3*ny*pow(nz,2))*pow(nz,3) + 
    35*pow(nx,4)*(3*nz*ry*pow(ny,2) + rz*pow(ny,3) - 3*ny*rz*pow(nz,2) - ry*pow(nz,3)) + 
    42*pow(nx,2)*pow(nz,2)*(-5*nz*ry*pow(ny,2) - 5*rz*pow(ny,3) + 5*ny*rz*pow(nz,2) + ry*pow(nz,3)) + 
    (21*nz*ry*pow(ny,2) + 35*rz*pow(ny,3) - 21*ny*rz*pow(nz,2) - 3*ry*pow(nz,3))*pow(nz,4),
   7*(ny*rx*pow(nz,2)*(-3*pow(ny,4) + 10*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 3*ny*rx*pow(nx,2)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
      5*pow(nx,3)*(-4*nz*rz*pow(ny,3) + ry*pow(ny,4) - 6*ry*pow(ny,2)*pow(nz,2) + 4*ny*rz*pow(nz,3) + ry*pow(nz,4)) - 
      nx*nz*(15*nz*ry*pow(ny,4) + 6*rz*pow(ny,5) - 40*rz*pow(ny,3)*pow(nz,2) - 30*ry*pow(ny,2)*pow(nz,3) + 18*ny*rz*pow(nz,4) + 3*ry*pow(nz,5))),
   rx*pow(nz,3)*(-35*pow(ny,4) + 42*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 21*nz*rx*pow(nx,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
    35*pow(nx,3)*(4*nz*ry*pow(ny,3) + rz*pow(ny,4) - 6*rz*pow(ny,2)*pow(nz,2) - 4*ny*ry*pow(nz,3) + rz*pow(nz,4)) - 
    7*nx*pow(nz,2)*(20*nz*ry*pow(ny,3) + 15*rz*pow(ny,4) - 30*rz*pow(ny,2)*pow(nz,2) - 12*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)),
   21*pow(nx,2)*(-5*nz*rz*pow(ny,4) + ry*pow(ny,5) - 10*ry*pow(ny,3)*pow(nz,2) + 10*rz*pow(ny,2)*pow(nz,3) + 5*ny*ry*pow(nz,4) - rz*pow(nz,5)) + 
    7*nx*rx*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) + 
    nz*(-21*nz*ry*pow(ny,5) - 7*rz*pow(ny,6) + 70*rz*pow(ny,4)*pow(nz,2) + 70*ry*pow(ny,3)*pow(nz,3) - 63*rz*pow(ny,2)*pow(nz,4) - 
       21*ny*ry*pow(nz,5) + 4*rz*pow(nz,6)),14*nx*ny*nz*rx*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
    21*pow(nx,2)*(5*nz*ry*pow(ny,4) + rz*pow(ny,5) - 10*rz*pow(ny,3)*pow(nz,2) - 10*ry*pow(ny,2)*pow(nz,3) + 5*ny*rz*pow(nz,4) + ry*pow(nz,5)) - 
    pow(nz,2)*(35*nz*ry*pow(ny,4) + 21*rz*pow(ny,5) - 70*rz*pow(ny,3)*pow(nz,2) - 42*ry*pow(ny,2)*pow(nz,3) + 21*ny*rz*pow(nz,4) + 3*ry*pow(nz,5)),
   -21*nz*(nz*rx + 2*nx*rz)*pow(ny,5) + 7*nx*ry*pow(ny,6) + rx*pow(ny,7) - 105*nx*ry*pow(ny,4)*pow(nz,2) + 
    35*(nz*rx + 4*nx*rz)*pow(ny,3)*pow(nz,3) + 105*nx*ry*pow(ny,2)*pow(nz,4) - 7*ny*(nz*rx + 6*nx*rz)*pow(nz,5) - 7*nx*ry*pow(nz,6),
   42*nx*nz*ry*pow(ny,5) + 7*(nz*rx + nx*rz)*pow(ny,6) - 35*(nz*rx + 3*nx*rz)*pow(ny,4)*pow(nz,2) - 140*nx*ry*pow(ny,3)*pow(nz,3) + 
    21*(nz*rx + 5*nx*rz)*pow(ny,2)*pow(nz,4) + 42*nx*ny*ry*pow(nz,5) - (nz*rx + 7*nx*rz)*pow(nz,6),
   -7*nz*rz*pow(ny,6) + ry*pow(ny,7) - 21*ry*pow(ny,5)*pow(nz,2) + 35*rz*pow(ny,4)*pow(nz,3) + 35*ry*pow(ny,3)*pow(nz,4) - 
    21*rz*pow(ny,2)*pow(nz,5) - 7*ny*ry*pow(nz,6) + rz*pow(nz,7),
   7*nz*ry*pow(ny,6) + rz*pow(ny,7) - 21*rz*pow(ny,5)*pow(nz,2) - 35*ry*pow(ny,4)*pow(nz,3) + 35*rz*pow(ny,3)*pow(nz,4) + 
    21*ry*pow(ny,2)*pow(nz,5) - 7*ny*rz*pow(nz,6) - ry*pow(nz,7),
   -12*nz*tx*tz*pow(nx,5) + 40*tx*tz*pow(nx,3)*pow(nz,3) - 12*nx*tx*tz*pow(nz,5) + pow(nx,6)*(pow(tx,2) - pow(tz,2)) - 
    15*pow(nx,4)*pow(nz,2)*(pow(tx,2) - pow(tz,2)) + 15*pow(nx,2)*pow(nz,4)*(pow(tx,2) - pow(tz,2)) + pow(nz,6)*(-pow(tx,2) + pow(tz,2)),
   2*(-15*nz*tx*(nz*ty + 2*ny*tz)*pow(nx,4) + tx*ty*pow(nx,6) + 15*tx*(nz*ty + 4*ny*tz)*pow(nx,2)*pow(nz,3) - tx*(nz*ty + 6*ny*tz)*pow(nz,5) + 
      10*pow(nx,3)*pow(nz,2)*(2*nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 3*nx*pow(nz,4)*(2*nz*ty*tz - 5*ny*pow(tx,2) + 5*ny*pow(tz,2)) - 
      3*pow(nx,5)*(2*nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   2*tx*tz*pow(nx,6) - 30*tx*tz*pow(nx,4)*pow(nz,2) + 30*tx*tz*pow(nx,2)*pow(nz,4) - 2*tx*tz*pow(nz,6) + 6*nz*pow(nx,5)*(pow(tx,2) - pow(tz,2)) - 
    20*pow(nx,3)*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + 6*nx*pow(nz,5)*(pow(tx,2) - pow(tz,2)),
   12*tx*(ny*ty - nz*tz)*pow(nx,5) + 40*nz*tx*pow(nx,3)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 
    12*nx*tx*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) - 
    pow(nz,4)*(12*ny*nz*ty*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 15*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,2)*pow(nz,2)*(8*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,4)*(-4*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,6)*(pow(ty,2) - pow(tz,2)),2*(6*tx*(nz*ty + ny*tz)*pow(nx,5) + ty*tz*pow(nx,6) - 20*tx*(nz*ty + 3*ny*tz)*pow(nx,3)*pow(nz,2) + 
      6*nx*tx*(nz*ty + 5*ny*tz)*pow(nz,4) + 15*pow(nx,2)*pow(nz,3)*(nz*ty*tz - 2*ny*pow(tx,2) + 2*ny*pow(tz,2)) - 
      pow(nz,5)*(nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 15*nz*pow(nx,4)*(nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   2*(15*tx*pow(nx,4)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
      tx*pow(nz,3)*(15*nz*ty*pow(ny,2) + 20*tz*pow(ny,3) - 18*ny*tz*pow(nz,2) - 3*ty*pow(nz,3)) + 
      30*nz*tx*pow(nx,2)*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      3*nx*pow(nz,2)*(20*nz*ty*tz*pow(ny,2) - 6*ty*tz*pow(nz,3) + 5*ny*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 
         10*pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 10*pow(nx,3)*
       (-6*nz*ty*tz*pow(ny,2) + 4*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2))) - 
      3*pow(nx,5)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2)))),
   2*(15*tx*pow(nx,4)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 30*tx*pow(nx,2)*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 
      3*tx*(2*ny*nz*ty + 5*tz*pow(ny,2) - tz*pow(nz,2))*pow(nz,4) + 
      3*nx*pow(nz,3)*(10*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 10*pow(ny,2)*(pow(tx,2) - pow(tz,2))) - 
      10*nz*pow(nx,3)*(6*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
      3*pow(nx,5)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2)))),
   40*tx*pow(nx,3)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
    12*nx*nz*tx*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
    pow(nz,2)*(40*nz*ty*tz*pow(ny,3) - 36*ny*ty*tz*pow(nz,3) + 15*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 
       3*pow(nz,4)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 15*pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,2)*(-8*nz*ty*tz*pow(ny,3) + 16*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       6*pow(ny,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,4)*(-4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2))),
   2*(20*tx*pow(nx,3)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      12*nx*tx*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      pow(nz,3)*(15*nz*ty*tz*pow(ny,2) - 3*ty*tz*pow(nz,3) + 3*ny*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 
         10*pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 30*nz*pow(nx,2)*
       (-3*nz*ty*tz*pow(ny,2) + ty*tz*pow(nz,3) - ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 
      15*pow(nx,4)*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2)))),
   2*(15*tx*pow(nx,2)*(-4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4)) - 
      nz*tx*(15*nz*ty*pow(ny,4) + 6*tz*pow(ny,5) - 40*tz*pow(ny,3)*pow(nz,2) - 30*ty*pow(ny,2)*pow(nz,3) + 18*ny*tz*pow(nz,4) + 3*ty*pow(nz,5)) + 
      3*nx*(-10*nz*ty*tz*pow(ny,4) + 40*ty*tz*pow(ny,2)*pow(nz,3) - 6*ty*tz*pow(nz,5) + 5*ny*pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
         10*pow(ny,3)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,5)*(pow(tx,2) - pow(tz,2))) + 
      10*pow(nx,3)*(-6*nz*ty*tz*pow(ny,2) + 2*ty*tz*pow(nz,3) + pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 3*ny*pow(nz,2)*(-pow(ty,2) + pow(tz,2)))),
   2*(tx*pow(nz,2)*(-20*nz*ty*pow(ny,3) - 15*tz*pow(ny,4) + 30*tz*pow(ny,2)*pow(nz,2) + 12*ny*ty*pow(nz,3) - 3*tz*pow(nz,4)) + 
      15*tx*pow(nx,2)*(4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4)) + 
      3*nx*nz*(-20*nz*ty*tz*pow(ny,3) + 20*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
         10*pow(ny,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 5*pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
      10*pow(nx,3)*(2*ty*tz*pow(ny,3) - 6*ny*ty*tz*pow(nz,2) + 3*nz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,3)*(-pow(ty,2) + pow(tz,2)))),
   12*ty*(nx*tx - nz*tz)*pow(ny,5) + 40*nz*ty*pow(ny,3)*(-3*nx*nz*tx - 3*tz*pow(nx,2) + 2*tz*pow(nz,2)) + 
    12*ny*ty*(5*nx*nz*tx + 10*tz*pow(nx,2) - 3*tz*pow(nz,2))*pow(nz,3) + pow(ny,6)*(pow(tx,2) - pow(tz,2)) - 
    15*pow(ny,4)*(4*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    15*pow(ny,2)*pow(nz,2)*(8*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 6*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) - 
    pow(nz,4)*(12*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) + 15*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   2*(15*ty*pow(ny,4)*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + 30*ty*pow(ny,2)*pow(nz,2)*(-2*nx*nz*tx - 3*tz*pow(nx,2) + tz*pow(nz,2)) + 
      3*ty*(2*nx*nz*tx + 5*tz*pow(nx,2) - tz*pow(nz,2))*pow(nz,4) + 3*pow(ny,5)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
      10*nz*pow(ny,3)*(6*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 3*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
      3*ny*pow(nz,3)*(10*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 10*pow(nx,2)*(-pow(ty,2) + pow(tz,2)))),
   2*(-15*nz*ty*(nz*tx + 2*nx*tz)*pow(ny,4) + tx*ty*pow(ny,6) + 15*ty*(nz*tx + 4*nx*tz)*pow(ny,2)*pow(nz,3) - ty*(nz*tx + 6*nx*tz)*pow(nz,5) + 
      10*pow(ny,3)*pow(nz,2)*(2*nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 3*ny*pow(nz,4)*(2*nz*tx*tz - 5*nx*pow(ty,2) + 5*nx*pow(tz,2)) - 
      3*pow(ny,5)*(2*nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   2*(6*ty*(nz*tx + nx*tz)*pow(ny,5) + tx*tz*pow(ny,6) - 20*ty*(nz*tx + 3*nx*tz)*pow(ny,3)*pow(nz,2) + 6*ny*ty*(nz*tx + 5*nx*tz)*pow(nz,4) + 
      15*pow(ny,2)*pow(nz,3)*(nz*tx*tz - 2*nx*pow(ty,2) + 2*nx*pow(tz,2)) - pow(nz,5)*(nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 
      15*nz*pow(ny,4)*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   -12*nz*ty*tz*pow(ny,5) + 40*ty*tz*pow(ny,3)*pow(nz,3) - 12*ny*ty*tz*pow(nz,5) + pow(ny,6)*(pow(ty,2) - pow(tz,2)) - 
    15*pow(ny,4)*pow(nz,2)*(pow(ty,2) - pow(tz,2)) + 15*pow(ny,2)*pow(nz,4)*(pow(ty,2) - pow(tz,2)) + pow(nz,6)*(-pow(ty,2) + pow(tz,2)),
   2*ty*tz*pow(ny,6) - 30*ty*tz*pow(ny,4)*pow(nz,2) + 30*ty*tz*pow(ny,2)*pow(nz,4) - 2*ty*tz*pow(nz,6) + 6*nz*pow(ny,5)*(pow(ty,2) - pow(tz,2)) - 
    20*pow(ny,3)*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + 6*ny*pow(nz,5)*(pow(ty,2) - pow(tz,2)),
   -6*nz*(rz*tx + rx*tz)*pow(nx,5) + (rx*tx - rz*tz)*pow(nx,6) - 15*(rx*tx - rz*tz)*pow(nx,4)*pow(nz,2) + 20*(rz*tx + rx*tz)*pow(nx,3)*pow(nz,3) + 
    15*(rx*tx - rz*tz)*pow(nx,2)*pow(nz,4) - 6*nx*(rz*tx + rx*tz)*pow(nz,5) + (-(rx*tx) + rz*tz)*pow(nz,6),
   -15*nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz)*pow(nx,4) - 6*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,5) + 
    (ry*tx + rx*ty)*pow(nx,6) + 20*(-3*ny*rx*tx + nz*rz*ty + nz*ry*tz + 3*ny*rz*tz)*pow(nx,3)*pow(nz,2) + 
    15*(nz*ry*tx + 4*ny*rz*tx + nz*rx*ty + 4*ny*rx*tz)*pow(nx,2)*pow(nz,3) - 6*nx*(-5*ny*rx*tx + nz*rz*ty + nz*ry*tz + 5*ny*rz*tz)*pow(nz,4) - 
    (nz*ry*tx + 6*ny*rz*tx + nz*rx*ty + 6*ny*rx*tz)*pow(nz,5),
   6*nz*(rx*tx - rz*tz)*pow(nx,5) + (rz*tx + rx*tz)*pow(nx,6) - 15*(rz*tx + rx*tz)*pow(nx,4)*pow(nz,2) - 20*(rx*tx - rz*tz)*pow(nx,3)*pow(nz,3) + 
    15*(rz*tx + rx*tz)*pow(nx,2)*pow(nz,4) + 6*nx*(rx*tx - rz*tz)*pow(nz,5) - (rz*tx + rx*tz)*pow(nz,6),
   6*(ny*(ry*tx + rx*ty) - nz*(rz*tx + rx*tz))*pow(nx,5) + (ry*ty - rz*tz)*pow(nx,6) - 
    20*nz*pow(nx,3)*(3*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(nz,2)) + 
    15*pow(nx,2)*pow(nz,2)*(4*ny*nz*(rz*ty + ry*tz) + (-6*rx*tx + 6*rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2)) + 
    15*pow(nx,4)*(-2*ny*nz*(rz*ty + ry*tz) + (rx*tx - rz*tz)*pow(ny,2) - (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    6*nx*(5*ny*nz*(ry*tx + rx*ty) + 10*(rz*tx + rx*tz)*pow(ny,2) - 3*(rz*tx + rx*tz)*pow(nz,2))*pow(nz,3) - 
    (6*ny*nz*(rz*ty + ry*tz) - 15*(rx*tx - rz*tz)*pow(ny,2) + (3*rx*tx + ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,4),
   -15*nz*(-2*ny*rx*tx + nz*rz*ty + nz*ry*tz + 2*ny*rz*tz)*pow(nx,4) + 6*(nz*ry*tx + ny*rz*tx + nz*rx*ty + ny*rx*tz)*pow(nx,5) + 
    (rz*ty + ry*tz)*pow(nx,6) - 20*(nz*ry*tx + 3*ny*rz*tx + nz*rx*ty + 3*ny*rx*tz)*pow(nx,3)*pow(nz,2) + 
    15*(-4*ny*rx*tx + nz*rz*ty + nz*ry*tz + 4*ny*rz*tz)*pow(nx,2)*pow(nz,3) + 6*nx*(nz*ry*tx + 5*ny*rz*tx + nz*rx*ty + 5*ny*rx*tz)*pow(nz,4) - 
    (-6*ny*rx*tx + nz*rz*ty + nz*ry*tz + 6*ny*rz*tz)*pow(nz,5),
   -6*(-(ny*ry*ty) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,5) + 
    15*pow(nx,4)*(-2*ny*nz*(rz*tx + rx*tz) + (ry*tx + rx*ty)*pow(ny,2) - (ry*tx + rx*ty)*pow(nz,2)) + 
    pow(nz,3)*(15*nz*(ry*tx + rx*ty)*pow(ny,2) + 20*(rz*tx + rx*tz)*pow(ny,3) - 18*ny*(rz*tx + rx*tz)*pow(nz,2) - 3*(ry*tx + rx*ty)*pow(nz,3)) + 
    30*nz*pow(nx,2)*(-3*nz*(ry*tx + rx*ty)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(ny,3) + 4*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
    20*pow(nx,3)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (rx*tx - rz*tz)*pow(ny,3) - 3*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + 
       2*(rz*ty + ry*tz)*pow(nz,3)) - 6*nx*pow(nz,2)*(-10*nz*(rz*ty + ry*tz)*pow(ny,2) + 10*(rx*tx - rz*tz)*pow(ny,3) - 
       5*ny*(2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2) + 3*(rz*ty + ry*tz)*pow(nz,3)),
   6*(nz*ry*ty + ny*rz*ty + ny*ry*tz - nz*rz*tz)*pow(nx,5) + 
    15*pow(nx,4)*(2*ny*nz*(ry*tx + rx*ty) + (rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    30*pow(nx,2)*pow(nz,2)*(2*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    20*nz*pow(nx,3)*(3*ny*nz*(rz*ty + ry*tz) + (-3*rx*tx + 3*rz*tz)*pow(ny,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    6*nx*(5*ny*nz*(rz*ty + ry*tz) - 10*(rx*tx - rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3) + 
    3*(2*ny*nz*(ry*tx + rx*ty) + 5*(rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2))*pow(nz,4),
   15*pow(nx,4)*(-2*ny*nz*(rz*ty + ry*tz) + (ry*ty - rz*tz)*pow(ny,2) + (-(ry*ty) + rz*tz)*pow(nz,2)) + 
    20*pow(nx,3)*(-3*nz*(rz*tx + rx*tz)*pow(ny,2) + (ry*tx + rx*ty)*pow(ny,3) - 3*ny*(ry*tx + rx*ty)*pow(nz,2) + (rz*tx + rx*tz)*pow(nz,3)) - 
    6*nx*nz*(10*nz*(ry*tx + rx*ty)*pow(ny,3) + 5*(rz*tx + rx*tz)*pow(ny,4) - 20*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       10*ny*(ry*tx + rx*ty)*pow(nz,3) + 3*(rz*tx + rx*tz)*pow(nz,4)) + 
    15*pow(nx,2)*(-4*nz*(rz*ty + ry*tz)*pow(ny,3) + (rx*tx - rz*tz)*pow(ny,4) - 6*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,2)*pow(nz,2) + 
       8*ny*(rz*ty + ry*tz)*pow(nz,3) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4)) + 
    pow(nz,2)*(20*nz*(rz*ty + ry*tz)*pow(ny,3) - 15*(rx*tx - rz*tz)*pow(ny,4) + 15*(2*rx*tx + ry*ty - 3*rz*tz)*pow(ny,2)*pow(nz,2) - 
       18*ny*(rz*ty + ry*tz)*pow(nz,3) - 3*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,4)),
   15*pow(nx,4)*(2*ny*nz*(ry*ty - rz*tz) + (rz*ty + ry*tz)*pow(ny,2) - (rz*ty + ry*tz)*pow(nz,2)) + 
    20*pow(nx,3)*(3*nz*(ry*tx + rx*ty)*pow(ny,2) + (rz*tx + rx*tz)*pow(ny,3) - 3*ny*(rz*tx + rx*tz)*pow(nz,2) - (ry*tx + rx*ty)*pow(nz,3)) + 
    12*nx*pow(nz,2)*(-5*nz*(ry*tx + rx*ty)*pow(ny,2) - 5*(rz*tx + rx*tz)*pow(ny,3) + 5*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
    pow(nz,3)*(15*nz*(rz*ty + ry*tz)*pow(ny,2) - 20*(rx*tx - rz*tz)*pow(ny,3) + 6*ny*(2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2) - 
       3*(rz*ty + ry*tz)*pow(nz,3)) + 30*nz*pow(nx,2)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + 2*(rx*tx - rz*tz)*pow(ny,3) - 
       2*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3)),
   20*pow(nx,3)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (ry*ty - rz*tz)*pow(ny,3) + 3*ny*(-(ry*ty) + rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3)) + 
    15*pow(nx,2)*(-4*nz*(rz*tx + rx*tz)*pow(ny,3) + (ry*tx + rx*ty)*pow(ny,4) - 6*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,2) + 
       4*ny*(rz*tx + rx*tz)*pow(nz,3) + (ry*tx + rx*ty)*pow(nz,4)) - 
    nz*(15*nz*(ry*tx + rx*ty)*pow(ny,4) + 6*(rz*tx + rx*tz)*pow(ny,5) - 40*(rz*tx + rx*tz)*pow(ny,3)*pow(nz,2) - 
       30*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,3) + 18*ny*(rz*tx + rx*tz)*pow(nz,4) + 3*(ry*tx + rx*ty)*pow(nz,5)) + 
    6*nx*(-5*nz*(rz*ty + ry*tz)*pow(ny,4) + (rx*tx - rz*tz)*pow(ny,5) - 10*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,3)*pow(nz,2) + 
       20*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,3) + 5*ny*(rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4) - 3*(rz*ty + ry*tz)*pow(nz,5)),
   20*pow(nx,3)*(3*nz*(ry*ty - rz*tz)*pow(ny,2) + (rz*ty + ry*tz)*pow(ny,3) - 3*ny*(rz*ty + ry*tz)*pow(nz,2) + (-(ry*ty) + rz*tz)*pow(nz,3)) + 
    15*pow(nx,2)*(4*nz*(ry*tx + rx*ty)*pow(ny,3) + (rz*tx + rx*tz)*pow(ny,4) - 6*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       4*ny*(ry*tx + rx*ty)*pow(nz,3) + (rz*tx + rx*tz)*pow(nz,4)) - 
    pow(nz,2)*(20*nz*(ry*tx + rx*ty)*pow(ny,3) + 15*(rz*tx + rx*tz)*pow(ny,4) - 30*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       12*ny*(ry*tx + rx*ty)*pow(nz,3) + 3*(rz*tx + rx*tz)*pow(nz,4)) + 
    6*nx*nz*(-10*nz*(rz*ty + ry*tz)*pow(ny,3) + 5*(rx*tx - rz*tz)*pow(ny,4) - 10*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,2)*pow(nz,2) + 
       10*ny*(rz*ty + ry*tz)*pow(nz,3) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4)),
   6*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz))*pow(ny,5) + (rx*tx - rz*tz)*pow(ny,6) - 
    20*nz*pow(ny,3)*(3*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - 2*(rz*ty + ry*tz)*pow(nz,2)) + 
    15*pow(ny,2)*pow(nz,2)*(4*nx*nz*(rz*tx + rx*tz) + 6*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2)) - 
    15*pow(ny,4)*(2*nx*nz*(rz*tx + rx*tz) + (-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    6*ny*(5*nx*nz*(ry*tx + rx*ty) + 10*(rz*ty + ry*tz)*pow(nx,2) - 3*(rz*ty + ry*tz)*pow(nz,2))*pow(nz,3) - 
    (6*nx*nz*(rz*tx + rx*tz) + 15*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 3*ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,4),
   6*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz)*pow(ny,5) + 
    15*pow(ny,4)*(2*nx*nz*(ry*tx + rx*ty) + (rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    30*pow(ny,2)*pow(nz,2)*(2*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    20*nz*pow(ny,3)*(3*nx*nz*(rz*tx + rx*tz) + 3*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    6*ny*(5*nx*nz*(rz*tx + rx*tz) + 10*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3) + 
    3*(2*nx*nz*(ry*tx + rx*ty) + 5*(rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2))*pow(nz,4),
   -15*nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz)*pow(ny,4) - 6*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz)*pow(ny,5) + 
    (ry*tx + rx*ty)*pow(ny,6) + 20*(nz*rz*tx - 3*nx*ry*ty + nz*rx*tz + 3*nx*rz*tz)*pow(ny,3)*pow(nz,2) + 
    15*(nz*ry*tx + nz*rx*ty + 4*nx*rz*ty + 4*nx*ry*tz)*pow(ny,2)*pow(nz,3) - 6*ny*(nz*rz*tx - 5*nx*ry*ty + nz*rx*tz + 5*nx*rz*tz)*pow(nz,4) - 
    (nz*ry*tx + nz*rx*ty + 6*nx*rz*ty + 6*nx*ry*tz)*pow(nz,5),
   -15*nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz)*pow(ny,4) + 6*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz)*pow(ny,5) + 
    (rz*tx + rx*tz)*pow(ny,6) - 20*(nz*ry*tx + nz*rx*ty + 3*nx*rz*ty + 3*nx*ry*tz)*pow(ny,3)*pow(nz,2) + 
    15*(nz*rz*tx - 4*nx*ry*ty + nz*rx*tz + 4*nx*rz*tz)*pow(ny,2)*pow(nz,3) + 6*ny*(nz*ry*tx + nz*rx*ty + 5*nx*rz*ty + 5*nx*ry*tz)*pow(nz,4) - 
    (nz*rz*tx - 6*nx*ry*ty + nz*rx*tz + 6*nx*rz*tz)*pow(nz,5),
   -6*nz*(rz*ty + ry*tz)*pow(ny,5) + (ry*ty - rz*tz)*pow(ny,6) - 15*(ry*ty - rz*tz)*pow(ny,4)*pow(nz,2) + 20*(rz*ty + ry*tz)*pow(ny,3)*pow(nz,3) + 
    15*(ry*ty - rz*tz)*pow(ny,2)*pow(nz,4) - 6*ny*(rz*ty + ry*tz)*pow(nz,5) + (-(ry*ty) + rz*tz)*pow(nz,6),
   6*nz*(ry*ty - rz*tz)*pow(ny,5) + (rz*ty + ry*tz)*pow(ny,6) - 15*(rz*ty + ry*tz)*pow(ny,4)*pow(nz,2) - 20*(ry*ty - rz*tz)*pow(ny,3)*pow(nz,3) + 
    15*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,4) + 6*ny*(ry*ty - rz*tz)*pow(nz,5) - (rz*ty + ry*tz)*pow(nz,6),
   -10*tx*pow(nx,3)*pow(nz,2)*(pow(tx,2) - 3*pow(tz,2)) + 5*nx*tx*pow(nz,4)*(pow(tx,2) - 3*pow(tz,2)) + 
    5*nz*tz*pow(nx,4)*(-3*pow(tx,2) + pow(tz,2)) - 10*tz*pow(nx,2)*pow(nz,3)*(-3*pow(tx,2) + pow(tz,2)) + tz*pow(nz,5)*(-3*pow(tx,2) + pow(tz,2)) + 
    pow(nx,5)*(pow(tx,3) - 3*tx*pow(tz,2)),5*tx*pow(nx,4)*(-6*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 
    tx*pow(nz,4)*(-6*nz*ty*tz + 5*ny*(pow(tx,2) - 3*pow(tz,2))) + 3*ty*pow(nx,5)*(pow(tx,2) - pow(tz,2)) + 
    30*tx*pow(nx,2)*pow(nz,2)*(2*nz*ty*tz - ny*pow(tx,2) + 3*ny*pow(tz,2)) + 
    5*nx*pow(nz,3)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - 4*ny*tz*(-3*pow(tx,2) + pow(tz,2))) - 
    10*nz*pow(nx,3)*(6*ny*tz*pow(tx,2) + 3*nz*ty*(pow(tx,2) - pow(tz,2)) - 2*ny*pow(tz,3)),
   5*nz*tx*pow(nx,4)*(pow(tx,2) - 3*pow(tz,2)) - 10*tx*pow(nx,2)*pow(nz,3)*(pow(tx,2) - 3*pow(tz,2)) + tx*pow(nz,5)*(pow(tx,2) - 3*pow(tz,2)) + 
    10*tz*pow(nx,3)*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) - 5*nx*tz*pow(nz,4)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,5)*(3*tz*pow(tx,2) - pow(tz,3)),
   5*nx*tx*pow(nz,2)*(24*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    10*tx*pow(nx,3)*(-12*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    3*tx*pow(nx,5)*(pow(ty,2) - pow(tz,2)) + 5*pow(nx,4)*(3*ny*ty*(pow(tx,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    pow(nz,3)*(15*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2)) + 
       10*pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 30*nz*pow(nx,2)*
     (tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 3*ny*nz*ty*(-pow(tx,2) + pow(tz,2)) + pow(ny,2)*(-3*tz*pow(tx,2) + pow(tz,3))),
   6*tx*ty*tz*pow(nx,5) + 20*nz*tx*pow(nx,3)*(-3*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 
    10*nx*tx*pow(nz,3)*(3*nz*ty*tz - 2*ny*pow(tx,2) + 6*ny*pow(tz,2)) + 
    pow(nz,4)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - 5*ny*tz*(-3*pow(tx,2) + pow(tz,2))) - 
    30*pow(nx,2)*pow(nz,2)*(nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2))) + 
    5*pow(nx,4)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2))),
   tx*pow(nz,2)*(60*nz*ty*tz*pow(ny,2) - 18*ty*tz*pow(nz,3) + 5*ny*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 
       10*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + 10*tx*pow(nx,2)*
     (-18*nz*ty*tz*pow(ny,2) + 12*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + 
    pow(nx,5)*(pow(ty,3) - 3*ty*pow(tz,2)) + 5*nx*nz*(ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 
       18*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 12*ny*tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
       4*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) - 15*tx*pow(nx,4)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2))) + 
    10*pow(nx,3)*(-(ty*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       2*ny*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))),
   tx*pow(nz,3)*(30*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 10*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) - 
    10*nz*tx*pow(nx,2)*(18*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    15*tx*pow(nx,4)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2))) - 
    15*nx*pow(nz,2)*(4*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)) + 
       pow(ny,2)*(6*tz*pow(tx,2) - 2*pow(tz,3))) + 10*pow(nx,3)*
     (6*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
    pow(nx,5)*(3*tz*pow(ty,2) - pow(tz,3)),5*nx*tx*(-24*nz*ty*tz*pow(ny,3) + 48*ny*ty*tz*pow(nz,3) + 
       pow(nz,4)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 6*pow(ny,2)*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
       pow(ny,4)*(pow(tx,2) - 3*pow(tz,2))) + 30*tx*pow(nx,3)*
     (-4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2))) + 
    30*pow(nx,2)*(-(ny*ty*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 
       tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + nz*tz*pow(ny,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    nz*(5*ny*ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 30*nz*ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 
       30*tz*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 5*tz*pow(ny,4)*(-3*pow(tx,2) + pow(tz,2)) + 
       3*tz*pow(nz,4)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 5*pow(nx,4)*(nz*tz*(-3*pow(ty,2) + pow(tz,2)) + ny*(pow(ty,3) - 3*ty*pow(tz,2))),
   20*nx*nz*tx*(-9*nz*ty*tz*pow(ny,2) + 3*ty*tz*pow(nz,3) - ny*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
       pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + 60*tx*pow(nx,3)*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2))) + 
    pow(nz,2)*(ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 30*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       15*ny*tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 10*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) + 
    10*pow(nx,2)*(-(ty*pow(nz,3)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 9*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       3*ny*tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,3)*(3*tz*pow(tx,2) - pow(tz,3))) + 
    5*pow(nx,4)*(3*ny*tz*pow(ty,2) + nz*pow(ty,3) - 3*nz*ty*pow(tz,2) - ny*pow(tz,3)),
   ty*pow(nz,2)*(60*nz*tx*tz*pow(nx,2) - 18*tx*tz*pow(nz,3) + 5*nx*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 
       10*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 10*ty*pow(ny,2)*
     (-18*nz*tx*tz*pow(nx,2) + 12*tx*tz*pow(nz,3) - 3*nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(ny,5)*(pow(tx,3) - 3*tx*pow(tz,2)) - 15*ty*pow(ny,4)*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))) + 
    5*ny*nz*(tx*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 18*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       12*nx*tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 4*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
    10*pow(ny,3)*(tx*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 2*nx*nz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       3*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),20*ny*nz*ty*
     (-9*nz*tx*tz*pow(nx,2) + 3*tx*tz*pow(nz,3) - nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 
    60*ty*pow(ny,3)*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))) + 
    pow(nz,2)*(tx*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 30*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       15*nx*tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 10*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
    10*pow(ny,2)*(tx*pow(nz,3)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 9*nz*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    5*pow(ny,4)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)),
   5*ny*ty*pow(nz,2)*(24*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
    10*ty*pow(ny,3)*(12*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
    3*ty*pow(ny,5)*(pow(tx,2) - pow(tz,2)) + 5*pow(ny,4)*(3*nx*tx*(pow(ty,2) - pow(tz,2)) + nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    pow(nz,3)*(15*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2)) + 
       10*pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + 30*nz*pow(ny,2)*
     (tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 3*nx*nz*tx*(-pow(ty,2) + pow(tz,2)) + pow(nx,2)*(-3*tz*pow(ty,2) + pow(tz,3))),
   ty*pow(nz,3)*(30*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 10*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
    10*nz*ty*pow(ny,2)*(18*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 3*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
    15*ty*pow(ny,4)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
    15*ny*pow(nz,2)*(4*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)) + 
       pow(nx,2)*(6*tz*pow(ty,2) - 2*pow(tz,3))) + 10*pow(ny,3)*
     (6*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + 
    pow(ny,5)*(3*tz*pow(tx,2) - pow(tz,3)),5*ty*pow(ny,4)*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
    ty*pow(nz,4)*(-6*nz*tx*tz + 5*nx*(pow(ty,2) - 3*pow(tz,2))) + 3*tx*pow(ny,5)*(pow(ty,2) - pow(tz,2)) + 
    30*ty*pow(ny,2)*pow(nz,2)*(2*nz*tx*tz - nx*pow(ty,2) + 3*nx*pow(tz,2)) + 
    5*ny*pow(nz,3)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - 4*nx*tz*(-3*pow(ty,2) + pow(tz,2))) - 
    10*nz*pow(ny,3)*(6*nx*tz*pow(ty,2) + 3*nz*tx*(pow(ty,2) - pow(tz,2)) - 2*nx*pow(tz,3)),
   6*tx*ty*tz*pow(ny,5) + 20*nz*ty*pow(ny,3)*(-3*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
    10*ny*ty*pow(nz,3)*(3*nz*tx*tz - 2*nx*pow(ty,2) + 6*nx*pow(tz,2)) + 
    pow(nz,4)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - 5*nx*tz*(-3*pow(ty,2) + pow(tz,2))) - 
    30*pow(ny,2)*pow(nz,2)*(nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2))) + 
    5*pow(ny,4)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2))),
   -10*ty*pow(ny,3)*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2)) + 5*ny*ty*pow(nz,4)*(pow(ty,2) - 3*pow(tz,2)) + 
    5*nz*tz*pow(ny,4)*(-3*pow(ty,2) + pow(tz,2)) - 10*tz*pow(ny,2)*pow(nz,3)*(-3*pow(ty,2) + pow(tz,2)) + tz*pow(nz,5)*(-3*pow(ty,2) + pow(tz,2)) + 
    pow(ny,5)*(pow(ty,3) - 3*ty*pow(tz,2)),5*nz*ty*pow(ny,4)*(pow(ty,2) - 3*pow(tz,2)) - 10*ty*pow(ny,2)*pow(nz,3)*(pow(ty,2) - 3*pow(tz,2)) + 
    ty*pow(nz,5)*(pow(ty,2) - 3*pow(tz,2)) + 10*tz*pow(ny,3)*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) - 5*ny*tz*pow(nz,4)*(-3*pow(ty,2) + pow(tz,2)) + 
    pow(ny,5)*(3*tz*pow(ty,2) - pow(tz,3)),pow(nx,5)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - 
    5*nz*pow(nx,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 10*pow(nx,2)*pow(nz,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    10*pow(nx,3)*pow(nz,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) - 5*nx*pow(nz,4)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    pow(nz,5)*(-2*rx*tx*tz + rz*(-pow(tx,2) + pow(tz,2))),pow(nx,5)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    5*pow(nx,4)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    10*pow(nx,2)*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) - 
    pow(nz,4)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(10*rz*tx*tz - 5*rx*pow(tx,2) + 5*rx*pow(tz,2))) - 
    10*nz*pow(nx,3)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 2*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
    5*nx*pow(nz,3)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 4*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(nz,5)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + pow(nx,5)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
    10*pow(nx,3)*pow(nz,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 5*nx*pow(nz,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
    5*nz*pow(nx,4)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 10*pow(nx,2)*pow(nz,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),
   10*pow(nx,3)*(-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) - pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + pow(nx,5)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    5*nx*pow(nz,2)*(8*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*ry*tx*ty - 6*rz*tx*tz + 2*rx*pow(tx,2) + rx*pow(ty,2) - 3*rx*pow(tz,2)) + 
       pow(ny,2)*(12*rz*tx*tz - 6*rx*pow(tx,2) + 6*rx*pow(tz,2))) + 
    5*pow(nx,4)*(-(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       ny*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    pow(nz,3)*(-5*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
       pow(nz,2)*(6*rx*tx*tz + 2*ry*ty*tz + 3*rz*pow(tx,2) + rz*pow(ty,2) - 4*rz*pow(tz,2)) - 10*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))
       ) + 10*nz*pow(nx,2)*(3*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
    ,2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,5) - 20*nz*pow(nx,3)*
     (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*rz*tx*tz - rx*pow(tx,2) + rx*pow(tz,2))) + 
    10*nx*pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(4*rz*tx*tz - 2*rx*pow(tx,2) + 2*rx*pow(tz,2))) + 
    5*pow(nx,4)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
    10*pow(nx,2)*pow(nz,2)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 3*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
    pow(nz,4)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 5*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   10*pow(nx,2)*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*ny*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    pow(nx,5)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + 5*pow(nx,4)*
     (-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    10*pow(nx,3)*(-(pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
       2*ny*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)))
      + 5*nx*nz*(pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       4*ny*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,2)*(20*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) - 6*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       5*ny*pow(nz,2)*(2*tx*(ry*ty - 3*rz*tz) + rx*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 10*pow(ny,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))
       ),pow(nx,5)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    5*pow(nx,4)*(2*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    10*nz*pow(nx,2)*(6*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) + 
    10*pow(nx,3)*(-(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       pow(ny,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 2*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    5*nx*pow(nz,2)*(4*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
      + pow(nz,3)*(10*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*tx*(ry*ty - 3*rz*tz) + rx*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       10*pow(ny,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))),
   5*nx*(-8*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) + 16*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       pow(nz,4)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       6*pow(ny,2)*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,4)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    10*pow(nx,3)*(-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(ny,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
       pow(nz,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    10*pow(nx,2)*(pow(nz,3)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       3*ny*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       3*nz*pow(ny,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    nz*(5*ny*pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       10*pow(ny,2)*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       3*pow(nz,4)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 5*pow(ny,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
       10*nz*pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    5*pow(nx,4)*(nz*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + ny*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   20*nx*nz*(-3*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       ny*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    20*pow(nx,3)*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
       ny*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    pow(nz,2)*(pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       5*ny*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       10*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 10*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    10*pow(nx,2)*(-(pow(nz,3)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
       3*ny*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
       3*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    5*pow(nx,4)*(nz*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2)) + ny*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   10*pow(ny,2)*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
    pow(ny,5)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) + 5*ny*nz*
     (pow(nz,3)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
       4*nx*pow(nz,2)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 6*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    10*pow(ny,3)*(pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))
       ) + 5*pow(ny,4)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,2)*(20*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - 6*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       5*nx*pow(nz,2)*(2*ty*(rx*tx - 3*rz*tz) + ry*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 10*pow(nx,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))
       ),20*ny*nz*(-3*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
    pow(nz,2)*(pow(nz,3)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
       5*nx*pow(nz,2)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       10*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 10*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    20*pow(ny,3)*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
       nx*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    5*pow(ny,4)*(nz*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + nx*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
    10*pow(ny,2)*(pow(nz,3)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*nx*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*nz*pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,3)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,5)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) - 
    5*pow(ny,4)*(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       nx*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    5*ny*pow(nz,2)*(8*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 6*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2)))
      - pow(nz,3)*(-5*nx*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 6*ry*ty*tz + rz*pow(tx,2) + 3*rz*pow(ty,2) - 4*rz*pow(tz,2)) - 10*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))
       ) + 10*nz*pow(ny,2)*(3*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 10*pow(ny,3)*(4*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,5)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 5*pow(ny,4)*
     (2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,3)*(10*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 
       10*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    5*ny*pow(nz,2)*(4*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 10*nz*pow(ny,2)*(6*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))) - 
    10*pow(ny,3)*(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,2)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,5)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) - 
    5*pow(ny,4)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    10*pow(ny,2)*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 3*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
    pow(nz,4)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 5*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
    10*nz*pow(ny,3)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 2*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
    5*ny*pow(nz,3)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 4*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,5) - 20*nz*pow(ny,3)*
     (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    10*ny*pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 2*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    5*pow(ny,4)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) - 
    10*pow(ny,2)*pow(nz,2)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 3*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
    pow(nz,4)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 5*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   pow(ny,5)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 5*nz*pow(ny,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    10*pow(ny,2)*pow(nz,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 10*pow(ny,3)*pow(nz,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) - 
    5*ny*pow(nz,4)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + pow(nz,5)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))),
   pow(nz,5)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) + pow(ny,5)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    10*pow(ny,3)*pow(nz,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 5*ny*pow(nz,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    5*nz*pow(ny,4)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + 10*pow(ny,2)*pow(nz,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),
   -16*nz*tx*tz*pow(nx,3)*(pow(tx,2) - pow(tz,2)) + 16*nx*tx*tz*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + 
    pow(nx,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 6*pow(nx,2)*pow(nz,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    pow(nz,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   4*(-6*nz*tx*pow(nx,2)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*ny*tz*(pow(tx,2) - pow(tz,2))) + 
      tx*pow(nz,3)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 4*ny*tz*(pow(tx,2) - pow(tz,2))) + tx*ty*pow(nx,4)*(pow(tx,2) - 3*pow(tz,2)) + 
      nx*pow(nz,2)*(-4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) - 3*ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nx,3)*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*(tx*tz*pow(nx,4)*(pow(tx,2) - pow(tz,2)) + tx*tz*pow(nz,4)*(pow(tx,2) - pow(tz,2)) + 6*tx*tz*pow(nx,2)*pow(nz,2)*(-pow(tx,2) + pow(tz,2)) + 
      nz*pow(nx,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - nx*pow(nz,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   2*(8*tx*pow(nx,3)*(ny*ty*(pow(tx,2) - 3*pow(tz,2)) + nz*tz*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
      8*nx*nz*tx*(3*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
      pow(nx,4)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      pow(nz,2)*(-8*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) - 3*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
      3*pow(nx,2)*(8*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(4*tx*pow(nx,3)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + ny*tz*(pow(tx,2) - pow(tz,2))) - 
      4*nx*tx*pow(nz,2)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*ny*tz*(pow(tx,2) - pow(tz,2))) - ty*tz*pow(nx,4)*(-3*pow(tx,2) + pow(tz,2)) + 
      pow(nz,3)*(-(nz*ty*tz*(-3*pow(tx,2) + pow(tz,2))) - ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      3*nz*pow(nx,2)*(2*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*(6*tx*pow(nx,2)*(-(ty*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) - 
         2*ny*nz*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + tx*ty*pow(nx,4)*(pow(ty,2) - 3*pow(tz,2)) + 
      nz*tx*(ty*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 6*nz*ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) + 
         4*ny*tz*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + pow(ny,3)*(-4*tz*pow(tx,2) + 4*pow(tz,3))) + 
      2*pow(nx,3)*(-2*nz*ty*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         ny*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nx*(4*ty*tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 12*nz*ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2)) + 
         pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         3*ny*pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(-(tx*tz*pow(nx,4)*(-3*pow(ty,2) + pow(tz,2))) + 6*tx*pow(nx,2)*
       (2*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      tx*pow(nz,2)*(-4*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + 
         pow(ny,2)*(-6*tz*pow(tx,2) + 6*pow(tz,3))) + 2*pow(nx,3)*
       (-2*ny*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nx*nz*(12*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   16*ny*ty*(-3*nx*tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
       tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 3*nz*tz*pow(nx,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
    16*nx*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 16*nz*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 
    16*ty*pow(ny,3)*(nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(ny,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + pow(nx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    6*pow(nx,2)*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(nz,4)*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) - 
    6*pow(ny,2)*(8*nx*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   4*(ty*(-4*nx*tx*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 4*nz*tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
         tz*pow(nz,4)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*tz*pow(nx,2)*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         tz*pow(nx,4)*(pow(ty,2) - pow(tz,2))) + 6*ty*pow(ny,2)*
       (2*nx*nz*tx*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
      pow(ny,3)*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
      ny*(12*nx*tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 4*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) - 
         6*nz*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,3)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(-6*ty*pow(ny,2)*(tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         2*nx*nz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + tx*ty*pow(ny,4)*(pow(tx,2) - 3*pow(tz,2)) + 
      nz*ty*(tx*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         4*nx*tz*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 4*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 
      2*pow(ny,3)*(-2*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         nx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(4*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 12*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         3*nx*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(tx*tz*pow(ny,4)*(pow(tx,2) - pow(tz,2)) - 4*ny*ty*(tx*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 
         3*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 4*ty*pow(ny,3)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)) + 
      6*pow(ny,2)*(-(tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) - tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         nx*nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nz*(tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 6*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         nx*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(8*ty*pow(ny,3)*(nx*tx*(pow(ty,2) - 3*pow(tz,2)) + nz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) - 
      8*ny*nz*ty*(3*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
      pow(ny,4)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      pow(nz,2)*(-8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) - 
      3*pow(ny,2)*(-8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(-(ty*tz*pow(ny,4)*(-3*pow(tx,2) + pow(tz,2))) + 6*ty*pow(ny,2)*
       (2*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
      ty*pow(nz,2)*(-4*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
         pow(nx,2)*(-6*tz*pow(ty,2) + 6*pow(tz,3))) + 2*pow(ny,3)*
       (-2*nx*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ny*nz*(12*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(-6*nz*ty*pow(ny,2)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 2*nx*tz*(pow(ty,2) - pow(tz,2))) + 
      ty*pow(nz,3)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 4*nx*tz*(pow(ty,2) - pow(tz,2))) + tx*ty*pow(ny,4)*(pow(ty,2) - 3*pow(tz,2)) + 
      ny*pow(nz,2)*(-4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 3*nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(ny,3)*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   4*(4*ty*pow(ny,3)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + nx*tz*(pow(ty,2) - pow(tz,2))) - 
      4*ny*ty*pow(nz,2)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*nx*tz*(pow(ty,2) - pow(tz,2))) - tx*tz*pow(ny,4)*(-3*pow(ty,2) + pow(tz,2)) + 
      pow(nz,3)*(-(nz*tx*tz*(-3*pow(ty,2) + pow(tz,2))) - nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      3*nz*pow(ny,2)*(2*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   -16*nz*ty*tz*pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 16*ny*ty*tz*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + 
    pow(ny,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 6*pow(ny,2)*pow(nz,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    pow(nz,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*(ty*tz*pow(ny,4)*(pow(ty,2) - pow(tz,2)) + ty*tz*pow(nz,4)*(pow(ty,2) - pow(tz,2)) + 6*ty*tz*pow(ny,2)*pow(nz,2)*(-pow(ty,2) + pow(tz,2)) + 
      nz*pow(ny,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ny*pow(nz,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   pow(nx,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    6*pow(nx,2)*pow(nz,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nz,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    4*nz*pow(nx,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    4*nx*pow(nz,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)),
   pow(nx,4)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    6*nz*pow(nx,2)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    pow(nz,3)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       4*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    4*pow(nx,3)*(ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) - 
    4*nx*pow(nz,2)*(3*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   4*nz*pow(nx,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    4*nx*pow(nz,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nx,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
    6*pow(nx,2)*pow(nz,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    pow(nz,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)),
   pow(nx,4)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    pow(nz,2)*(pow(nz,2)*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) + 3*ry*ty*(pow(tx,2) - pow(tz,2)) + 
          rz*tz*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2))) + 
       4*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) - 
       6*pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    6*pow(nx,2)*(pow(nz,2)*(-(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
          3*ry*ty*(-pow(tx,2) + pow(tz,2))) - 2*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    4*pow(nx,3)*(-(nz*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2)))) + 
       ny*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    4*nx*nz*(pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       3*ny*nz*(6*rz*tx*ty*tz - 3*rx*ty*pow(tx,2) - ry*pow(tx,3) + 3*ry*tx*pow(tz,2) + 3*rx*ty*pow(tz,2)) + 
       pow(ny,2)*(-9*rx*tz*pow(tx,2) - 3*rz*pow(tx,3) + 9*rz*tx*pow(tz,2) + 3*rx*pow(tz,3))),
   pow(nx,4)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    4*pow(nx,3)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) - 
    4*nx*pow(nz,2)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       3*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    6*nz*pow(nx,2)*(2*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) + 
    pow(nz,3)*(nz*(6*rx*tx*ty*tz + 3*rz*ty*pow(tx,2) + 3*ry*tz*pow(tx,2) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) - 
       4*ny*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3))),
   pow(nx,4)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    4*pow(nx,3)*(-(nz*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
       ny*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)))) + 
    6*pow(nx,2)*(-(pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)))) - 
       2*ny*nz*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2))) + 
       pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    4*nx*(pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       3*nz*pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       3*ny*pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
          3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2))))) + 
    nz*(pow(nz,3)*(ty*(-18*rz*tx*tz + rx*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       4*ny*pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       6*nz*pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       4*pow(ny,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))),
   pow(nx,4)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    4*pow(nx,3)*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       ny*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2)))) - 
    4*nx*nz*(3*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) - 
       3*pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + pow(nz,2)*(pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       4*ny*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       6*pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    6*pow(nx,2)*(-(pow(nz,2)*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2)))) + 
       2*ny*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))),
   4*nx*pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) + 
    pow(nz,4)*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + 
       3*rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    4*nz*pow(nx,3)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    pow(ny,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nx,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    6*pow(nx,2)*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
       3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))) + 
    4*ny*(pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       3*nx*pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,3)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
       3*nz*pow(nx,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) - 
    6*pow(ny,2)*(2*nx*nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + 4*pow(ny,3)*(nx*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(nz,4)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
    6*pow(nx,2)*pow(nz,2)*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) - 
    4*nx*pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
    4*nz*pow(nx,3)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) - 
    4*ny*(3*nx*pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       3*nz*pow(nx,2)*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nx,3)*(-6*ry*tx*ty*tz + rx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*rz*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,3)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + pow(nx,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    4*pow(ny,3)*(nz*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nx*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    6*pow(ny,2)*(pow(nx,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       2*nx*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(ny,4)*(ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))) - 
    4*pow(ny,3)*(nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       nx*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2)))) - 
    6*pow(ny,2)*(pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       2*nx*nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) + 
    4*ny*(pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
       3*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       3*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) + 
    nz*(4*nx*pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,3)*(ry*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + ty*(-18*rz*tx*tz + 3*rx*pow(tx,2) + 2*rx*pow(ty,2) - 9*rx*pow(tz,2))) - 
       6*nz*pow(nx,2)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       4*pow(nx,3)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   -6*pow(ny,2)*(pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       2*nx*nz*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nx,2)*(-6*ry*tx*ty*tz + rx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)))) + 
    4*pow(ny,3)*(nx*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    nz*(pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
       6*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       4*pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       4*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) + 
    pow(ny,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
    4*ny*(pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       3*nz*pow(nx,2)*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*nx*pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,3)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,4)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    4*pow(ny,3)*(nx*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) + 
    pow(nz,2)*(pow(nz,2)*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
          rz*tz*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2))) + 
       4*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       6*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2)))) - 
    6*pow(ny,2)*(2*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     + 4*ny*nz*(pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       3*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,4)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    4*pow(ny,3)*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       nx*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)))) - 
    4*ny*nz*(3*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       3*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     - 6*pow(ny,2)*(2*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))) + 
    pow(nz,2)*(pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 
       4*nx*nz*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       6*pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,4)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
    6*nz*pow(ny,2)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    pow(nz,3)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       4*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    4*ny*pow(nz,2)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       3*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))) - 
    4*pow(ny,3)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,4)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    4*pow(ny,3)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) - 
    4*ny*pow(nz,2)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       3*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    pow(nz,3)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       4*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))) - 
    6*nz*pow(ny,2)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       2*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    6*pow(ny,2)*pow(nz,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(nz,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    4*nz*pow(ny,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    4*ny*pow(nz,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   4*nz*pow(ny,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    4*ny*pow(nz,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(ny,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) - 
    6*pow(ny,2)*pow(nz,2)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    pow(nz,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   -3*nz*tz*pow(nx,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + tz*pow(nz,3)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    3*nx*tx*pow(nz,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + pow(nx,3)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   5*ty*pow(nx,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    3*nx*nz*(2*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    tx*pow(nz,2)*(20*nz*ty*tz*(pow(tx,2) - pow(tz,2)) - 3*ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
    3*tx*pow(nx,2)*(20*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))),
   -3*nx*tz*pow(nz,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 3*nz*tx*pow(nx,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
    tx*pow(nz,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + pow(nx,3)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)),
   10*tx*pow(nx,3)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    3*pow(nx,2)*(nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ny*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    3*nx*tx*(40*ny*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) - pow(nz,2)*
        (pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       pow(ny,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
    nz*(-15*ny*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       2*tz*pow(nz,2)*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
       3*pow(ny,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5))),
   20*tx*ty*tz*pow(nx,3)*(pow(tx,2) - pow(tz,2)) + 3*pow(nx,2)*
     (ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(3*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    6*nx*nz*tx*(10*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))),
   pow(ny,3)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
    15*ty*pow(ny,2)*(4*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + nx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*(tx*pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) - 
       10*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*nx*nz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    5*ty*(4*tx*tz*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 12*nz*tx*tz*pow(nx,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       2*pow(nx,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       3*nx*pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   -(tx*pow(nz,3)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)))) + 
    2*tz*pow(nx,3)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    30*nz*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*nx*tz*pow(nz,2)*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    30*ny*ty*(2*tx*tz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 2*tx*tz*pow(nz,2)*(-pow(tx,2) + pow(tz,2)) + 
       nx*nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    3*pow(ny,2)*(nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   5*ty*pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*tx*pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*nz*tz*pow(nx,2)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
    15*nx*tx*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    tz*pow(nz,3)*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
    3*pow(ny,2)*(nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       10*nx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    15*ny*ty*(8*nx*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))),
   15*ty*pow(ny,2)*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    5*ty*(12*nx*tx*tz*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 4*tx*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2)) - 
       6*nz*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,3)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(ny,3)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)) + 
    3*ny*(tz*pow(nz,2)*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       20*nx*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*pow(nx,2)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5))),
   10*tx*pow(ny,3)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    30*ty*pow(ny,2)*(-2*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       nx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*(-5*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*nx*nz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       5*tx*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    ty*(20*tx*tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 60*nz*tx*tz*pow(nx,2)*(-pow(ty,2) + pow(tz,2)) + 
       pow(nx,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
       3*nx*pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))),
   20*tx*ty*tz*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + tz*pow(nx,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    15*nz*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*nx*tz*pow(nz,2)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
    5*tx*pow(nz,3)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    6*pow(ny,2)*(nx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       5*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    60*ny*ty*(-(tx*tz*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + tx*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       nx*nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))),
   10*ty*pow(ny,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    3*pow(ny,2)*(nz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*nx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*ty*(40*nx*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
    nz*(-15*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*tz*pow(nz,2)*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
       3*pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5))),
   30*ty*pow(ny,2)*(2*nx*tx*tz*(pow(ty,2) - pow(tz,2)) + nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nz*ty*(60*nx*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + 3*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
       pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
    2*pow(ny,3)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)) + 
    3*ny*(tz*pow(nz,2)*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       10*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5))),
   5*tx*pow(ny,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*ny*nz*(2*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    ty*pow(nz,2)*(20*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - 3*nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
    3*ty*pow(ny,2)*(20*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   20*tx*ty*tz*pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 3*pow(ny,2)*
     (nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(3*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    6*ny*nz*ty*(10*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   -3*nz*tz*pow(ny,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + tz*pow(nz,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    3*ny*ty*pow(nz,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + pow(ny,3)*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   -3*ny*tz*pow(nz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 3*nz*ty*pow(ny,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
    ty*pow(nz,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + pow(ny,3)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)),
   pow(nx,3)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    3*nz*pow(nx,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   pow(nx,3)*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,2)*(4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       3*ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(nx,2)*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*nz*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       2*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   3*nz*pow(nx,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,3)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   2*pow(nx,3)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nz*(-3*ny*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
       3*pow(ny,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       2*pow(nz,2)*(2*tz*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
    3*nx*(-8*ny*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       pow(ny,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tx*(-(ry*ty*(pow(tx,2) - 3*pow(tz,2))) + rz*tz*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) - 
          rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    3*pow(nx,2)*(ny*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
       nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*pow(nx,3)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
    6*nx*nz*(2*nz*(-3*rx*ty*tz*pow(tx,2) - ry*tz*pow(tx,3) - rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + ry*tx*pow(tz,3) + rx*ty*pow(tz,3)) + 
       ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(nx,2)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       3*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*pow(nz,3)*(tz*(rx*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) + 
       rz*tx*ty*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 12*nz*pow(nx,2)*
     (tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
       rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + pow(ny,3)*
     (4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(nx,3)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
       ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    3*pow(ny,2)*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nx*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    3*ny*(2*pow(nx,2)*(-2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tx*(ry*ty*(pow(tx,2) - 3*pow(tz,2)) + rz*tz*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
          rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       2*nx*nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   6*nz*pow(nx,2)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(nx,3)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
       rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(4*tx*(-(ry*ty*(pow(tx,2) - 3*pow(tz,2))) + rz*tz*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) - 
       rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
       rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    6*ny*(2*pow(nx,2)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       2*pow(nz,2)*(-3*rx*ty*tz*pow(tx,2) - ry*tz*pow(tx,3) - rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + ry*tx*pow(tz,3) + rx*ty*pow(tz,3)) + 
       nx*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(ny,2)*(nz*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       nx*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   pow(ny,3)*(4*ty*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,3)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
       rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    3*nz*pow(nx,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(nz,3)*(4*tz*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
       rz*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) - 
    3*ny*(8*nx*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       2*pow(nx,2)*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
          ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    3*pow(ny,2)*(2*nx*(-2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   -12*nx*pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
       rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 4*pow(nx,3)*
     (tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(ny,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    6*nz*pow(nx,2)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,3)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
       ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    3*pow(ny,2)*(4*nx*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    3*ny*(-4*nx*nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       2*pow(nx,2)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
          rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*pow(nz,3)*(tz*(rx*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + ry*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
       rz*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 12*nz*pow(nx,2)*
     (tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(nx,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    2*pow(ny,3)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
       ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    6*pow(ny,2)*(2*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       nx*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*ny*(-(pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       2*nx*nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*pow(ny,3)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))) + 
    3*nz*pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
       rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,3)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
       rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    3*nx*pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    12*ny*(pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(-(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + ry*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*ty*tz*(-pow(ty,2) + pow(tz,2))) + 
       nx*nz*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    6*pow(ny,2)*(nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       nx*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
          rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   2*pow(ny,3)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    nz*(-3*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       3*pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*pow(nz,2)*(2*tz*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
    3*ny*(8*nx*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
          ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    3*pow(ny,2)*(-(nx*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*pow(ny,3)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
       rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    6*pow(ny,2)*(2*nx*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nz*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    nz*(-12*nx*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       3*pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
          ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    3*ny*(-2*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   pow(ny,3)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,2)*(4*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       3*nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*pow(ny,2)*(4*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*ny*nz*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   4*pow(ny,3)*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) - 
    6*ny*nz*(2*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(ny,2)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       3*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   pow(ny,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*pow(nz,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nz*pow(ny,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   3*nz*pow(ny,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(ny,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*pow(nz,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   -4*nx*nz*tx*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(nx,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    pow(nz,2)*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6)),
   2*(3*tx*ty*pow(nx,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
      nz*tx*(2*ny*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      nx*(-6*nz*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         ny*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))),
   -2*tx*tz*pow(nz,2)*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + pow(nx,2)*(6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5)) + 
    2*nx*nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   -4*nx*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
    12*ny*ty*(-(nz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + nx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    pow(nz,2)*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
       4*pow(tz,6)) + pow(ny,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    3*pow(nx,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*(3*ty*(-(tz*pow(nz,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 2*nx*nz*tx*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
         pow(nx,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5))) + 
      ny*(2*nx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
         nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))),
   2*(3*tx*ty*pow(ny,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
      ty*(-10*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         2*nx*nz*tz*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         tx*pow(nz,2)*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4))) + 
      ny*(-2*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         3*nx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   2*(2*tx*tz*pow(nx,2)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
      tx*tz*pow(nz,2)*(-3*pow(tx,4) + 30*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 9*pow(tz,4)) + 
      6*ny*ty*(nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
      pow(ny,2)*(3*tz*pow(tx,5) - 10*pow(tx,3)*pow(tz,3) + 3*tx*pow(tz,5)) + 
      3*nx*nz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   4*nx*nz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
    3*pow(nx,2)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    4*ny*ty*(nz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       10*nx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    3*pow(nz,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
       5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6)) + 
    3*pow(ny,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*(3*ty*tz*pow(ny,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
      ty*(tz*pow(nz,2)*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         20*nx*nz*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         2*pow(nx,2)*(15*tz*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,3) + 3*pow(tz,5))) + 
      ny*(4*nx*tx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
         3*nz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   2*(10*tx*ty*pow(ny,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      ny*(2*nz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
         3*nx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      ty*(3*tx*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
         2*nx*nz*tz*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         tx*pow(nz,2)*(-10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))))),
   2*(2*tx*tz*pow(ny,2)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
      tx*tz*pow(nz,2)*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
      3*tx*tz*pow(nx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      3*nx*nz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      4*ny*ty*(nx*tz*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         5*nz*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   -12*nx*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    3*pow(ny,2)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    4*ny*ty*(nz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       3*nx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) - 
    pow(nz,2)*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       4*pow(tz,6)) + pow(nx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*(2*ty*tz*pow(ny,2)*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      3*ny*(2*nx*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      ty*(tz*pow(nz,2)*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         6*nx*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + pow(nx,2)*(3*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + 3*pow(tz,5)))),
   2*(3*tx*ty*pow(ny,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
      nz*ty*(2*nx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      ny*(-6*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   2*(3*tx*tz*pow(ny,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      ny*(6*nz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 2*nx*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4))) + 
      nz*(-3*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   -4*ny*nz*ty*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(ny,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    pow(nz,2)*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)),
   -2*ty*tz*pow(nz,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + pow(ny,2)*(6*tz*pow(ty,5) - 20*pow(ty,3)*pow(tz,3) + 6*ty*pow(tz,5)) + 
    2*ny*nz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   pow(nz,2)*(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    pow(nx,2)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    2*nx*nz*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   pow(nx,2)*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    2*nx*(ny*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
       nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) - 
    nz*(2*ny*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   2*nx*nz*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    pow(nx,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    pow(nz,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   pow(nx,2)*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(rx*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
       2*rz*tz*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
    pow(ny,2)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    2*nx*nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))
      + 2*ny*(-(nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
            tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) + 
       nx*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   pow(nx,2)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    2*nx*nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    2*ny*(nz*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nx*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)))),
   pow(ny,2)*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    10*pow(nx,2)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    2*nx*nz*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    pow(nz,2)*(ry*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       5*ty*(-4*rz*tx*tz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    2*ny*(nx*(-5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 
                3*pow(tz,4))))),2*nx*nz*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(ny,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    2*pow(nx,2)*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))
      + 2*ny*(nx*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,2)*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,2)*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,2)*(-5*rx*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
       5*ry*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rz*tz*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) - 
    2*nx*nz*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    2*ny*(-10*nx*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   pow(ny,2)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    2*pow(nx,2)*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    20*nx*nz*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    2*ny*(nz*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
          5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))
         + 2*nx*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   10*pow(ny,2)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    2*nx*nz*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))\
     + pow(nx,2)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    pow(nz,2)*(5*ry*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ty*(-20*rz*tx*tz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) - 
    2*ny*(nx*(-5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   2*nx*nz*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,2)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    2*pow(ny,2)*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    4*ny*(nx*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       5*nz*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,2)*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,2)*(5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       2*rz*tz*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
       ry*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
    pow(nx,2)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    2*nx*nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    2*ny*(nz*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))
          ) - nx*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(nx,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    2*pow(ny,2)*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    pow(nz,2)*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))\
     + 2*nx*nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) + 
    2*ny*(nz*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
          5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))
         + nx*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,2)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    2*ny*(nx*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    nz*(2*nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(ny,2)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(2*nx*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) + 
    2*ny*(nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(nz,2)*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    pow(ny,2)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    2*ny*nz*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   2*ny*nz*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    pow(ny,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    pow(nz,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   nz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    nx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   7*ty*(-2*nz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       nx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   -(nx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))) + 
    nz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   7*nx*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    7*ny*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    nz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 21*pow(ty,2)*pow(tz,4) + 
       4*pow(tz,6)),tz*(14*nx*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ny*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    7*nz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   7*(ty*(-2*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         nx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
      ny*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))),
   tz*(14*ny*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nx*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    7*nz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   7*ny*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    nz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
       21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
    7*nx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   tz*(28*nx*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ny*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    7*nz*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   7*(ty*(2*nz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         nx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
      ny*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   tz*(28*ny*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
    7*nz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   7*ny*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
    7*nx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    nz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       4*pow(tz,6)),tz*(14*nx*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ny*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
    7*nz*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)),
   ty*(-14*nz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
    7*ny*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   tz*(14*ny*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    7*nz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   nz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    ny*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6)),
   -(ny*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
    nz*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6)),
   nx*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(-2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6))),
   -2*nz*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
       tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
    ny*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))),
   nz*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))),
   -(nz*(2*tz*(3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
            rx*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
         rz*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
            4*pow(tz,6)))) + ny*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
          rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))\
     + nx*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*nx*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
       tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
    nz*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))),
   -2*nz*(rz*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) + 
       tz*(ry*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)))) + 
    ny*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   ny*(6*rx*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 6*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
       ry*(6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5))) + 
    nz*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   nx*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    nz*(2*tz*(rx*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          ry*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
          5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))) + 
    ny*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   4*nx*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    nz*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   ny*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    nx*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    2*nz*(tz*(ry*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*tx*ty*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))),
   nz*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    4*ny*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    nx*(3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   ny*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    nz*(2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ry*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))
      + nx*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*nx*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    ny*(3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    nz*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   -2*nz*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    ny*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   ny*(6*ry*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       2*rx*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4))) + 
    nz*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   ny*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(-2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   nz*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8),
   8*tx*ty*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)),
   8*tx*tz*(pow(tx,6) - 7*pow(tx,4)*pow(tz,2) + 7*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   4*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
      7*pow(ty,2)*pow(tz,6) + pow(tz,8)),-8*ty*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)),
   56*tx*ty*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   8*tx*tz*(7*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 14*pow(tx,2)*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   70*pow(ty,4)*pow(tz,4) + 70*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 84*pow(ty,2)*pow(tz,6) - 
    84*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 6*pow(tz,8),
   8*ty*tz*(35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 3*pow(tz,6)),
   56*tx*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)),
   8*tx*tz*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)),
   4*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
      21*pow(ty,2)*pow(tz,6) + pow(tz,8)),8*ty*tz*(-7*pow(ty,4)*pow(tz,2) + 14*pow(ty,2)*pow(tz,4) + 
      7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 3*pow(tz,6)),
   8*tx*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)),
   -8*tx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)),
   pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8),
   8*ty*tz*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   rz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    rx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   7*ty*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    ry*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   -(rx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))) + 
    rz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   7*rx*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    7*ry*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    rz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 21*pow(ty,2)*pow(tz,4) + 
       4*pow(tz,6)),tz*(14*rx*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ry*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    7*rz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   7*(ty*(-2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         rx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
      ry*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))),
   tz*(14*ry*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rx*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    7*rz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   7*ry*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    rz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
       21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
    7*rx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   tz*(28*rx*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ry*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    7*rz*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)),
   7*(ty*(2*rz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         rx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
      ry*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   tz*(28*ry*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
    7*rz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   7*ry*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
    7*rx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    rz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       4*pow(tz,6)),tz*(14*rx*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ry*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
    7*rz*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)),
   ty*(-14*rz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
    7*ry*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   tz*(14*ry*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    7*rz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   rz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    ry*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6)),
   -(ry*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
    rz*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6));

    degree = 9;
    AssertDimension(tensor_project[degree].P.rows(),components(degree));
    AssertDimension(degree,9);

    tensor_project[degree].P << pow(nx,9) - 36*pow(nx,7)*pow(nz,2) + 126*pow(nx,5)*pow(nz,4) - 84*pow(nx,3)*pow(nz,6) + 9*nx*pow(nz,8),
   9*ny*(pow(nx,8) - 28*pow(nx,6)*pow(nz,2) + 70*pow(nx,4)*pow(nz,4) - 28*pow(nx,2)*pow(nz,6) + pow(nz,8)),
   9*nz*pow(nx,8) - 84*pow(nx,6)*pow(nz,3) + 126*pow(nx,4)*pow(nz,5) - 36*pow(nx,2)*pow(nz,7) + pow(nz,9),
   36*nx*(pow(nx,6)*(pow(ny,2) - pow(nz,2)) + 7*pow(nx,4)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 7*pow(nx,2)*(5*pow(ny,2)*pow(nz,4) - pow(nz,6)) - 
      7*pow(ny,2)*pow(nz,6) + pow(nz,8)),72*nx*ny*nz*(pow(nx,6) - 7*pow(nx,4)*pow(nz,2) + 7*pow(nx,2)*pow(nz,4) - pow(nz,6)),
   12*ny*(7*pow(nx,6)*(pow(ny,2) - 3*pow(nz,2)) + 105*pow(nx,4)*pow(nz,2)*(-pow(ny,2) + pow(nz,2)) + 
      21*pow(nx,2)*(5*pow(ny,2)*pow(nz,4) - 3*pow(nz,6)) - 7*pow(ny,2)*pow(nz,6) + 3*pow(nz,8)),
   4*(21*pow(nx,6)*(3*nz*pow(ny,2) - pow(nz,3)) - 63*pow(nx,4)*(5*pow(ny,2)*pow(nz,3) - pow(nz,5)) + 
      27*pow(nx,2)*(7*pow(ny,2)*pow(nz,5) - pow(nz,7)) - 9*pow(ny,2)*pow(nz,7) + pow(nz,9)),
   18*nx*(35*pow(ny,4)*pow(nz,4) + 7*pow(nx,4)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 42*pow(ny,2)*pow(nz,6) - 
      14*pow(nx,2)*(5*pow(ny,4)*pow(nz,2) - 10*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 3*pow(nz,8)),
   24*nx*ny*nz*(21*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 21*pow(ny,2)*pow(nz,4) + pow(nx,2)*(-70*pow(ny,2)*pow(nz,2) + 42*pow(nz,4)) - 9*pow(nz,6)),
   18*ny*(7*pow(ny,4)*pow(nz,4) + 7*pow(nx,4)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) - 14*pow(ny,2)*pow(nz,6) - 
      14*pow(nx,2)*(3*pow(ny,4)*pow(nz,2) - 10*pow(ny,2)*pow(nz,4) + 3*pow(nz,6)) + 3*pow(nz,8)),
   6*(21*pow(ny,4)*pow(nz,5) + 21*pow(nx,4)*(5*nz*pow(ny,4) - 10*pow(ny,2)*pow(nz,3) + pow(nz,5)) - 18*pow(ny,2)*pow(nz,7) - 
      6*pow(nx,2)*(35*pow(ny,4)*pow(nz,3) - 42*pow(ny,2)*pow(nz,5) + 3*pow(nz,7)) + pow(nz,9)),
   12*nx*(7*pow(nx,2)*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) + 
      3*pow(nz,2)*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6))),
   24*nx*ny*nz*(-21*pow(ny,4)*pow(nz,2) + 42*pow(ny,2)*pow(nz,4) + 7*pow(nx,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) - 9*pow(nz,6)),
   36*ny*(-(pow(ny,6)*pow(nz,2)) + 7*pow(ny,4)*pow(nz,4) + pow(nx,2)*(pow(ny,6) - 21*pow(ny,4)*pow(nz,2) + 35*pow(ny,2)*pow(nz,4) - 7*pow(nz,6)) - 
      7*pow(ny,2)*pow(nz,6) + pow(nz,8)),4*(-21*pow(ny,6)*pow(nz,3) + 63*pow(ny,4)*pow(nz,5) + 
      9*pow(nx,2)*(7*nz*pow(ny,6) - 35*pow(ny,4)*pow(nz,3) + 21*pow(ny,2)*pow(nz,5) - pow(nz,7)) - 27*pow(ny,2)*pow(nz,7) + pow(nz,9)),
   9*nx*(pow(ny,8) - 28*pow(ny,6)*pow(nz,2) + 70*pow(ny,4)*pow(nz,4) - 28*pow(ny,2)*pow(nz,6) + pow(nz,8)),
   72*nx*ny*nz*(pow(ny,6) - 7*pow(ny,4)*pow(nz,2) + 7*pow(ny,2)*pow(nz,4) - pow(nz,6)),
   pow(ny,9) - 36*pow(ny,7)*pow(nz,2) + 126*pow(ny,5)*pow(nz,4) - 84*pow(ny,3)*pow(nz,6) + 9*ny*pow(nz,8),
   9*nz*pow(ny,8) - 84*pow(ny,6)*pow(nz,3) + 126*pow(ny,4)*pow(nz,5) - 36*pow(ny,2)*pow(nz,7) + pow(nz,9),
   -8*nz*tz*pow(nx,7) + tx*pow(nx,8) - 28*tx*pow(nx,6)*pow(nz,2) + 56*tz*pow(nx,5)*pow(nz,3) + 70*tx*pow(nx,4)*pow(nz,4) - 
    56*tz*pow(nx,3)*pow(nz,5) - 28*tx*pow(nx,2)*pow(nz,6) + 8*nx*tz*pow(nz,7) + tx*pow(nz,8),
   -28*nz*(nz*ty + 2*ny*tz)*pow(nx,6) + 8*ny*tx*pow(nx,7) + ty*pow(nx,8) - 168*ny*tx*pow(nx,5)*pow(nz,2) + 
    70*(nz*ty + 4*ny*tz)*pow(nx,4)*pow(nz,3) + 280*ny*tx*pow(nx,3)*pow(nz,4) - 28*(nz*ty + 6*ny*tz)*pow(nx,2)*pow(nz,5) - 56*nx*ny*tx*pow(nz,6) + 
    (nz*ty + 8*ny*tz)*pow(nz,7),8*nz*tx*pow(nx,7) + tz*pow(nx,8) - 28*tz*pow(nx,6)*pow(nz,2) - 56*tx*pow(nx,5)*pow(nz,3) + 
    70*tz*pow(nx,4)*pow(nz,4) + 56*tx*pow(nx,3)*pow(nz,5) - 28*tz*pow(nx,2)*pow(nz,6) - 8*nx*tx*pow(nz,7) + tz*pow(nz,8),
   4*(2*(ny*ty - nz*tz)*pow(nx,7) + 7*tx*pow(nx,6)*(pow(ny,2) - pow(nz,2)) + 35*tx*pow(nx,4)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
      14*nz*pow(nx,5)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 14*pow(nx,3)*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) - 
      21*tx*pow(nx,2)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4) + 2*nx*(-7*ny*nz*ty - 21*tz*pow(ny,2) + 4*tz*pow(nz,2))*pow(nz,5) + 
      tx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,6)),8*(7*ny*nz*tx*pow(nx,6) + (nz*ty + ny*tz)*pow(nx,7) - 7*(nz*ty + 3*ny*tz)*pow(nx,5)*pow(nz,2) - 
      35*ny*tx*pow(nx,4)*pow(nz,3) + 7*(nz*ty + 5*ny*tz)*pow(nx,3)*pow(nz,4) + 21*ny*tx*pow(nx,2)*pow(nz,5) - nx*(nz*ty + 7*ny*tz)*pow(nz,6) - 
      ny*tx*pow(nz,7)),4*(14*ny*tx*pow(nx,5)*(pow(ny,2) - 3*pow(nz,2)) - 140*ny*tx*pow(nx,3)*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
      7*pow(nx,6)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
      7*pow(nx,2)*pow(nz,3)*(15*nz*ty*pow(ny,2) + 20*tz*pow(ny,3) - 18*ny*tz*pow(nz,2) - 3*ty*pow(nz,3)) + 
      35*nz*pow(nx,4)*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 14*nx*tx*(5*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,4) + 
      (-7*nz*ty*pow(ny,2) - 14*tz*pow(ny,3) + 8*ny*tz*pow(nz,2) + ty*pow(nz,3))*pow(nz,5)),
   4*(-14*nz*tx*pow(nx,5)*(-3*pow(ny,2) + pow(nz,2)) + 7*pow(nx,6)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 
      35*pow(nx,4)*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 28*tx*pow(nx,3)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
      21*pow(nx,2)*(2*ny*nz*ty + 5*tz*pow(ny,2) - tz*pow(nz,2))*pow(nz,4) - 6*nx*tx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,5) + 
      (-2*ny*nz*ty - 7*tz*pow(ny,2) + tz*pow(nz,2))*pow(nz,6)),
   2*(28*pow(nx,5)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
      42*tx*pow(nx,2)*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
      35*tx*pow(nx,4)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) + tx*pow(nz,4)*(35*pow(ny,4) - 42*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) - 
      28*nz*pow(nx,3)*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
      4*nx*pow(nz,3)*(35*nz*ty*pow(ny,3) + 35*tz*pow(ny,4) - 63*tz*pow(ny,2)*pow(nz,2) - 21*ny*ty*pow(nz,3) + 6*tz*pow(nz,4))),
   8*(35*ny*nz*tx*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 7*pow(nx,5)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      14*pow(nx,3)*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      nx*(21*nz*ty*pow(ny,2) + 35*tz*pow(ny,3) - 21*ny*tz*pow(nz,2) - 3*ty*pow(nz,3))*pow(nz,4) + tx*(7*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,5) + 
      pow(nx,2)*(-70*tx*pow(ny,3)*pow(nz,3) + 42*ny*tx*pow(nz,5))),
   2*(-28*nx*ny*tx*pow(nz,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
      28*ny*tx*pow(nx,3)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
      35*pow(nx,4)*(-4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4)) - 
      14*nz*pow(nx,2)*(15*nz*ty*pow(ny,4) + 6*tz*pow(ny,5) - 40*tz*pow(ny,3)*pow(nz,2) - 30*ty*pow(ny,2)*pow(nz,3) + 18*ny*tz*pow(nz,4) + 
         3*ty*pow(nz,5)) + pow(nz,3)*(35*nz*ty*pow(ny,4) + 28*tz*pow(ny,5) - 84*tz*pow(ny,3)*pow(nz,2) - 42*ty*pow(ny,2)*pow(nz,3) + 
         24*ny*tz*pow(nz,4) + 3*ty*pow(nz,5))),56*nz*tx*pow(nx,3)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    8*nx*tx*pow(nz,3)*(35*pow(ny,4) - 42*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
    70*pow(nx,4)*(4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4)) + 
    2*pow(nz,4)*(28*nz*ty*pow(ny,3) + 35*tz*pow(ny,4) - 42*tz*pow(ny,2)*pow(nz,2) - 12*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) - 
    28*pow(nx,2)*pow(nz,2)*(20*nz*ty*pow(ny,3) + 15*tz*pow(ny,4) - 30*tz*pow(ny,2)*pow(nz,2) - 12*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)),
   4*(14*pow(nx,3)*(-5*nz*tz*pow(ny,4) + ty*pow(ny,5) - 10*ty*pow(ny,3)*pow(nz,2) + 10*tz*pow(ny,2)*pow(nz,3) + 5*ny*ty*pow(nz,4) - tz*pow(nz,5)) + 
      7*tx*pow(nx,2)*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) + 
      tx*pow(nz,2)*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 
      2*nx*nz*(-21*nz*ty*pow(ny,5) - 7*tz*pow(ny,6) + 70*tz*pow(ny,4)*pow(nz,2) + 70*ty*pow(ny,3)*pow(nz,3) - 63*tz*pow(ny,2)*pow(nz,4) - 
         21*ny*ty*pow(nz,5) + 4*tz*pow(nz,6))),8*(ny*tx*pow(nz,3)*(-7*pow(ny,4) + 14*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 
      7*ny*nz*tx*pow(nx,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
      7*pow(nx,3)*(5*nz*ty*pow(ny,4) + tz*pow(ny,5) - 10*tz*pow(ny,3)*pow(nz,2) - 10*ty*pow(ny,2)*pow(nz,3) + 5*ny*tz*pow(nz,4) + ty*pow(nz,5)) - 
      nx*pow(nz,2)*(35*nz*ty*pow(ny,4) + 21*tz*pow(ny,5) - 70*tz*pow(ny,3)*pow(nz,2) - 42*ty*pow(ny,2)*pow(nz,3) + 21*ny*tz*pow(nz,4) + 
         3*ty*pow(nz,5))),4*(2*nx*ny*tx*(pow(ny,6) - 21*pow(ny,4)*pow(nz,2) + 35*pow(ny,2)*pow(nz,4) - 7*pow(nz,6)) + 
      7*pow(nx,2)*(-6*nz*tz*pow(ny,5) + ty*pow(ny,6) - 15*ty*pow(ny,4)*pow(nz,2) + 20*tz*pow(ny,3)*pow(nz,3) + 15*ty*pow(ny,2)*pow(nz,4) - 
         6*ny*tz*pow(nz,5) - ty*pow(nz,6)) + nz*(-7*nz*ty*pow(ny,6) - 2*tz*pow(ny,7) + 28*tz*pow(ny,5)*pow(nz,2) + 35*ty*pow(ny,4)*pow(nz,3) - 
         42*tz*pow(ny,3)*pow(nz,4) - 21*ty*pow(ny,2)*pow(nz,5) + 8*ny*tz*pow(nz,6) + ty*pow(nz,7))),
   4*(-2*nx*nz*tx*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 
      7*pow(nx,2)*(6*nz*ty*pow(ny,5) + tz*pow(ny,6) - 15*tz*pow(ny,4)*pow(nz,2) - 20*ty*pow(ny,3)*pow(nz,3) + 15*tz*pow(ny,2)*pow(nz,4) + 
         6*ny*ty*pow(nz,5) - tz*pow(nz,6)) + pow(nz,2)*(-14*nz*ty*pow(ny,5) - 7*tz*pow(ny,6) + 35*tz*pow(ny,4)*pow(nz,2) + 
         28*ty*pow(ny,3)*pow(nz,3) - 21*tz*pow(ny,2)*pow(nz,4) - 6*ny*ty*pow(nz,5) + tz*pow(nz,6))),
   -28*nz*(nz*tx + 2*nx*tz)*pow(ny,6) + 8*nx*ty*pow(ny,7) + tx*pow(ny,8) - 168*nx*ty*pow(ny,5)*pow(nz,2) + 
    70*(nz*tx + 4*nx*tz)*pow(ny,4)*pow(nz,3) + 280*nx*ty*pow(ny,3)*pow(nz,4) - 28*(nz*tx + 6*nx*tz)*pow(ny,2)*pow(nz,5) - 56*nx*ny*ty*pow(nz,6) + 
    (nz*tx + 8*nx*tz)*pow(nz,7),8*(7*nx*nz*ty*pow(ny,6) + (nz*tx + nx*tz)*pow(ny,7) - 7*(nz*tx + 3*nx*tz)*pow(ny,5)*pow(nz,2) - 
      35*nx*ty*pow(ny,4)*pow(nz,3) + 7*(nz*tx + 5*nx*tz)*pow(ny,3)*pow(nz,4) + 21*nx*ty*pow(ny,2)*pow(nz,5) - ny*(nz*tx + 7*nx*tz)*pow(nz,6) - 
      nx*ty*pow(nz,7)),-8*nz*tz*pow(ny,7) + ty*pow(ny,8) - 28*ty*pow(ny,6)*pow(nz,2) + 56*tz*pow(ny,5)*pow(nz,3) + 70*ty*pow(ny,4)*pow(nz,4) - 
    56*tz*pow(ny,3)*pow(nz,5) - 28*ty*pow(ny,2)*pow(nz,6) + 8*ny*tz*pow(nz,7) + ty*pow(nz,8),
   8*nz*ty*pow(ny,7) + tz*pow(ny,8) - 28*tz*pow(ny,6)*pow(nz,2) - 56*ty*pow(ny,5)*pow(nz,3) + 70*tz*pow(ny,4)*pow(nz,4) + 
    56*ty*pow(ny,3)*pow(nz,5) - 28*tz*pow(ny,2)*pow(nz,6) - 8*ny*ty*pow(nz,7) + tz*pow(nz,8),
   -8*nz*rz*pow(nx,7) + rx*pow(nx,8) - 28*rx*pow(nx,6)*pow(nz,2) + 56*rz*pow(nx,5)*pow(nz,3) + 70*rx*pow(nx,4)*pow(nz,4) - 
    56*rz*pow(nx,3)*pow(nz,5) - 28*rx*pow(nx,2)*pow(nz,6) + 8*nx*rz*pow(nz,7) + rx*pow(nz,8),
   -28*nz*(nz*ry + 2*ny*rz)*pow(nx,6) + 8*ny*rx*pow(nx,7) + ry*pow(nx,8) - 168*ny*rx*pow(nx,5)*pow(nz,2) + 
    70*(nz*ry + 4*ny*rz)*pow(nx,4)*pow(nz,3) + 280*ny*rx*pow(nx,3)*pow(nz,4) - 28*(nz*ry + 6*ny*rz)*pow(nx,2)*pow(nz,5) - 56*nx*ny*rx*pow(nz,6) + 
    (nz*ry + 8*ny*rz)*pow(nz,7),8*nz*rx*pow(nx,7) + rz*pow(nx,8) - 28*rz*pow(nx,6)*pow(nz,2) - 56*rx*pow(nx,5)*pow(nz,3) + 
    70*rz*pow(nx,4)*pow(nz,4) + 56*rx*pow(nx,3)*pow(nz,5) - 28*rz*pow(nx,2)*pow(nz,6) - 8*nx*rx*pow(nz,7) + rz*pow(nz,8),
   4*(2*(ny*ry - nz*rz)*pow(nx,7) + 7*rx*pow(nx,6)*(pow(ny,2) - pow(nz,2)) + 35*rx*pow(nx,4)*pow(nz,2)*(-3*pow(ny,2) + pow(nz,2)) + 
      14*nz*pow(nx,5)*(-3*ny*nz*ry - 3*rz*pow(ny,2) + 2*rz*pow(nz,2)) + 14*pow(nx,3)*(5*ny*nz*ry + 10*rz*pow(ny,2) - 3*rz*pow(nz,2))*pow(nz,3) - 
      21*rx*pow(nx,2)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,4) + 2*nx*(-7*ny*nz*ry - 21*rz*pow(ny,2) + 4*rz*pow(nz,2))*pow(nz,5) + 
      rx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,6)),8*(7*ny*nz*rx*pow(nx,6) + (nz*ry + ny*rz)*pow(nx,7) - 7*(nz*ry + 3*ny*rz)*pow(nx,5)*pow(nz,2) - 
      35*ny*rx*pow(nx,4)*pow(nz,3) + 7*(nz*ry + 5*ny*rz)*pow(nx,3)*pow(nz,4) + 21*ny*rx*pow(nx,2)*pow(nz,5) - nx*(nz*ry + 7*ny*rz)*pow(nz,6) - 
      ny*rx*pow(nz,7)),4*(14*ny*rx*pow(nx,5)*(pow(ny,2) - 3*pow(nz,2)) - 140*ny*rx*pow(nx,3)*(pow(ny,2) - pow(nz,2))*pow(nz,2) + 
      7*pow(nx,6)*(-2*ny*nz*rz + ry*pow(ny,2) - ry*pow(nz,2)) + 
      7*pow(nx,2)*pow(nz,3)*(15*nz*ry*pow(ny,2) + 20*rz*pow(ny,3) - 18*ny*rz*pow(nz,2) - 3*ry*pow(nz,3)) + 
      35*nz*pow(nx,4)*(-3*nz*ry*pow(ny,2) - 2*rz*pow(ny,3) + 4*ny*rz*pow(nz,2) + ry*pow(nz,3)) + 14*nx*rx*(5*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,4) + 
      (-7*nz*ry*pow(ny,2) - 14*rz*pow(ny,3) + 8*ny*rz*pow(nz,2) + ry*pow(nz,3))*pow(nz,5)),
   4*(-14*nz*rx*pow(nx,5)*(-3*pow(ny,2) + pow(nz,2)) + 7*pow(nx,6)*(2*ny*nz*ry + rz*pow(ny,2) - rz*pow(nz,2)) + 
      35*pow(nx,4)*pow(nz,2)*(-2*ny*nz*ry - 3*rz*pow(ny,2) + rz*pow(nz,2)) + 28*rx*pow(nx,3)*(-5*pow(ny,2) + pow(nz,2))*pow(nz,3) + 
      21*pow(nx,2)*(2*ny*nz*ry + 5*rz*pow(ny,2) - rz*pow(nz,2))*pow(nz,4) - 6*nx*rx*(-7*pow(ny,2) + pow(nz,2))*pow(nz,5) + 
      (-2*ny*nz*ry - 7*rz*pow(ny,2) + rz*pow(nz,2))*pow(nz,6)),
   2*(28*pow(nx,5)*(-3*nz*rz*pow(ny,2) + ry*pow(ny,3) - 3*ny*ry*pow(nz,2) + rz*pow(nz,3)) - 
      42*rx*pow(nx,2)*pow(nz,2)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) + 
      35*rx*pow(nx,4)*(pow(ny,4) - 6*pow(ny,2)*pow(nz,2) + pow(nz,4)) + rx*pow(nz,4)*(35*pow(ny,4) - 42*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) - 
      28*nz*pow(nx,3)*(10*nz*ry*pow(ny,3) + 5*rz*pow(ny,4) - 20*rz*pow(ny,2)*pow(nz,2) - 10*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)) + 
      4*nx*pow(nz,3)*(35*nz*ry*pow(ny,3) + 35*rz*pow(ny,4) - 63*rz*pow(ny,2)*pow(nz,2) - 21*ny*ry*pow(nz,3) + 6*rz*pow(nz,4))),
   8*(35*ny*nz*rx*pow(nx,4)*(pow(ny,2) - pow(nz,2)) + 7*pow(nx,5)*(3*nz*ry*pow(ny,2) + rz*pow(ny,3) - 3*ny*rz*pow(nz,2) - ry*pow(nz,3)) + 
      14*pow(nx,3)*pow(nz,2)*(-5*nz*ry*pow(ny,2) - 5*rz*pow(ny,3) + 5*ny*rz*pow(nz,2) + ry*pow(nz,3)) + 
      nx*(21*nz*ry*pow(ny,2) + 35*rz*pow(ny,3) - 21*ny*rz*pow(nz,2) - 3*ry*pow(nz,3))*pow(nz,4) + rx*(7*pow(ny,3) - 3*ny*pow(nz,2))*pow(nz,5) + 
      pow(nx,2)*(-70*rx*pow(ny,3)*pow(nz,3) + 42*ny*rx*pow(nz,5))),
   2*(-28*nx*ny*rx*pow(nz,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
      28*ny*rx*pow(nx,3)*(pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 5*pow(nz,4)) + 
      35*pow(nx,4)*(-4*nz*rz*pow(ny,3) + ry*pow(ny,4) - 6*ry*pow(ny,2)*pow(nz,2) + 4*ny*rz*pow(nz,3) + ry*pow(nz,4)) - 
      14*nz*pow(nx,2)*(15*nz*ry*pow(ny,4) + 6*rz*pow(ny,5) - 40*rz*pow(ny,3)*pow(nz,2) - 30*ry*pow(ny,2)*pow(nz,3) + 18*ny*rz*pow(nz,4) + 
         3*ry*pow(nz,5)) + pow(nz,3)*(35*nz*ry*pow(ny,4) + 28*rz*pow(ny,5) - 84*rz*pow(ny,3)*pow(nz,2) - 42*ry*pow(ny,2)*pow(nz,3) + 
         24*ny*rz*pow(nz,4) + 3*ry*pow(nz,5))),56*nz*rx*pow(nx,3)*(5*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + pow(nz,4)) - 
    8*nx*rx*pow(nz,3)*(35*pow(ny,4) - 42*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
    70*pow(nx,4)*(4*nz*ry*pow(ny,3) + rz*pow(ny,4) - 6*rz*pow(ny,2)*pow(nz,2) - 4*ny*ry*pow(nz,3) + rz*pow(nz,4)) + 
    2*pow(nz,4)*(28*nz*ry*pow(ny,3) + 35*rz*pow(ny,4) - 42*rz*pow(ny,2)*pow(nz,2) - 12*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)) - 
    28*pow(nx,2)*pow(nz,2)*(20*nz*ry*pow(ny,3) + 15*rz*pow(ny,4) - 30*rz*pow(ny,2)*pow(nz,2) - 12*ny*ry*pow(nz,3) + 3*rz*pow(nz,4)),
   4*(14*pow(nx,3)*(-5*nz*rz*pow(ny,4) + ry*pow(ny,5) - 10*ry*pow(ny,3)*pow(nz,2) + 10*rz*pow(ny,2)*pow(nz,3) + 5*ny*ry*pow(nz,4) - rz*pow(nz,5)) + 
      7*rx*pow(nx,2)*(pow(ny,6) - 15*pow(ny,4)*pow(nz,2) + 15*pow(ny,2)*pow(nz,4) - pow(nz,6)) + 
      rx*pow(nz,2)*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 
      2*nx*nz*(-21*nz*ry*pow(ny,5) - 7*rz*pow(ny,6) + 70*rz*pow(ny,4)*pow(nz,2) + 70*ry*pow(ny,3)*pow(nz,3) - 63*rz*pow(ny,2)*pow(nz,4) - 
         21*ny*ry*pow(nz,5) + 4*rz*pow(nz,6))),8*(ny*rx*pow(nz,3)*(-7*pow(ny,4) + 14*pow(ny,2)*pow(nz,2) - 3*pow(nz,4)) + 
      7*ny*nz*rx*pow(nx,2)*(3*pow(ny,4) - 10*pow(ny,2)*pow(nz,2) + 3*pow(nz,4)) + 
      7*pow(nx,3)*(5*nz*ry*pow(ny,4) + rz*pow(ny,5) - 10*rz*pow(ny,3)*pow(nz,2) - 10*ry*pow(ny,2)*pow(nz,3) + 5*ny*rz*pow(nz,4) + ry*pow(nz,5)) - 
      nx*pow(nz,2)*(35*nz*ry*pow(ny,4) + 21*rz*pow(ny,5) - 70*rz*pow(ny,3)*pow(nz,2) - 42*ry*pow(ny,2)*pow(nz,3) + 21*ny*rz*pow(nz,4) + 
         3*ry*pow(nz,5))),4*(2*nx*ny*rx*(pow(ny,6) - 21*pow(ny,4)*pow(nz,2) + 35*pow(ny,2)*pow(nz,4) - 7*pow(nz,6)) + 
      7*pow(nx,2)*(-6*nz*rz*pow(ny,5) + ry*pow(ny,6) - 15*ry*pow(ny,4)*pow(nz,2) + 20*rz*pow(ny,3)*pow(nz,3) + 15*ry*pow(ny,2)*pow(nz,4) - 
         6*ny*rz*pow(nz,5) - ry*pow(nz,6)) + nz*(-7*nz*ry*pow(ny,6) - 2*rz*pow(ny,7) + 28*rz*pow(ny,5)*pow(nz,2) + 35*ry*pow(ny,4)*pow(nz,3) - 
         42*rz*pow(ny,3)*pow(nz,4) - 21*ry*pow(ny,2)*pow(nz,5) + 8*ny*rz*pow(nz,6) + ry*pow(nz,7))),
   4*(-2*nx*nz*rx*(-7*pow(ny,6) + 35*pow(ny,4)*pow(nz,2) - 21*pow(ny,2)*pow(nz,4) + pow(nz,6)) + 
      7*pow(nx,2)*(6*nz*ry*pow(ny,5) + rz*pow(ny,6) - 15*rz*pow(ny,4)*pow(nz,2) - 20*ry*pow(ny,3)*pow(nz,3) + 15*rz*pow(ny,2)*pow(nz,4) + 
         6*ny*ry*pow(nz,5) - rz*pow(nz,6)) + pow(nz,2)*(-14*nz*ry*pow(ny,5) - 7*rz*pow(ny,6) + 35*rz*pow(ny,4)*pow(nz,2) + 
         28*ry*pow(ny,3)*pow(nz,3) - 21*rz*pow(ny,2)*pow(nz,4) - 6*ny*ry*pow(nz,5) + rz*pow(nz,6))),
   -28*nz*(nz*rx + 2*nx*rz)*pow(ny,6) + 8*nx*ry*pow(ny,7) + rx*pow(ny,8) - 168*nx*ry*pow(ny,5)*pow(nz,2) + 
    70*(nz*rx + 4*nx*rz)*pow(ny,4)*pow(nz,3) + 280*nx*ry*pow(ny,3)*pow(nz,4) - 28*(nz*rx + 6*nx*rz)*pow(ny,2)*pow(nz,5) - 56*nx*ny*ry*pow(nz,6) + 
    (nz*rx + 8*nx*rz)*pow(nz,7),8*(7*nx*nz*ry*pow(ny,6) + (nz*rx + nx*rz)*pow(ny,7) - 7*(nz*rx + 3*nx*rz)*pow(ny,5)*pow(nz,2) - 
      35*nx*ry*pow(ny,4)*pow(nz,3) + 7*(nz*rx + 5*nx*rz)*pow(ny,3)*pow(nz,4) + 21*nx*ry*pow(ny,2)*pow(nz,5) - ny*(nz*rx + 7*nx*rz)*pow(nz,6) - 
      nx*ry*pow(nz,7)),-8*nz*rz*pow(ny,7) + ry*pow(ny,8) - 28*ry*pow(ny,6)*pow(nz,2) + 56*rz*pow(ny,5)*pow(nz,3) + 70*ry*pow(ny,4)*pow(nz,4) - 
    56*rz*pow(ny,3)*pow(nz,5) - 28*ry*pow(ny,2)*pow(nz,6) + 8*ny*rz*pow(nz,7) + ry*pow(nz,8),
   8*nz*ry*pow(ny,7) + rz*pow(ny,8) - 28*rz*pow(ny,6)*pow(nz,2) - 56*ry*pow(ny,5)*pow(nz,3) + 70*rz*pow(ny,4)*pow(nz,4) + 
    56*ry*pow(ny,3)*pow(nz,5) - 28*rz*pow(ny,2)*pow(nz,6) - 8*ny*ry*pow(nz,7) + rz*pow(nz,8),
   -14*nz*tx*tz*pow(nx,6) + 70*tx*tz*pow(nx,4)*pow(nz,3) - 42*tx*tz*pow(nx,2)*pow(nz,5) + 2*tx*tz*pow(nz,7) + pow(nx,7)*(pow(tx,2) - pow(tz,2)) - 
    21*pow(nx,5)*pow(nz,2)*(pow(tx,2) - pow(tz,2)) + 35*pow(nx,3)*pow(nz,4)*(pow(tx,2) - pow(tz,2)) + 7*nx*pow(nz,6)*(-pow(tx,2) + pow(tz,2)),
   -42*nz*tx*(nz*ty + 2*ny*tz)*pow(nx,5) + 2*tx*ty*pow(nx,7) + 70*tx*(nz*ty + 4*ny*tz)*pow(nx,3)*pow(nz,3) - 14*nx*tx*(nz*ty + 6*ny*tz)*pow(nz,5) + 
    35*pow(nx,4)*pow(nz,2)*(2*nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 21*pow(nx,2)*pow(nz,4)*(2*nz*ty*tz - 5*ny*pow(tx,2) + 5*ny*pow(tz,2)) + 
    pow(nz,6)*(2*nz*ty*tz - 7*ny*pow(tx,2) + 7*ny*pow(tz,2)) - 7*pow(nx,6)*(2*nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2))),
   2*tx*tz*pow(nx,7) - 42*tx*tz*pow(nx,5)*pow(nz,2) + 70*tx*tz*pow(nx,3)*pow(nz,4) - 14*nx*tx*tz*pow(nz,6) + 
    7*nz*pow(nx,6)*(pow(tx,2) - pow(tz,2)) - 35*pow(nx,4)*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + 21*pow(nx,2)*pow(nz,5)*(pow(tx,2) - pow(tz,2)) + 
    pow(nz,7)*(-pow(tx,2) + pow(tz,2)),14*tx*(ny*ty - nz*tz)*pow(nx,6) + 70*nz*tx*pow(nx,4)*(-3*ny*nz*ty - 3*tz*pow(ny,2) + 2*tz*pow(nz,2)) + 
    42*tx*pow(nx,2)*(5*ny*nz*ty + 10*tz*pow(ny,2) - 3*tz*pow(nz,2))*pow(nz,3) + 2*tx*(-7*ny*nz*ty - 21*tz*pow(ny,2) + 4*tz*pow(nz,2))*pow(nz,5) - 
    7*nx*pow(nz,4)*(12*ny*nz*ty*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 15*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    35*pow(nx,3)*pow(nz,2)*(8*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    21*pow(nx,5)*(-4*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    pow(nx,7)*(pow(ty,2) - pow(tz,2)),2*(7*tx*(nz*ty + ny*tz)*pow(nx,6) + ty*tz*pow(nx,7) - 35*tx*(nz*ty + 3*ny*tz)*pow(nx,4)*pow(nz,2) + 
      21*tx*(nz*ty + 5*ny*tz)*pow(nx,2)*pow(nz,4) - tx*(nz*ty + 7*ny*tz)*pow(nz,6) + 
      35*pow(nx,3)*pow(nz,3)*(nz*ty*tz - 2*ny*pow(tx,2) + 2*ny*pow(tz,2)) - 7*nx*pow(nz,5)*(nz*ty*tz - 3*ny*pow(tx,2) + 3*ny*pow(tz,2)) - 
      21*nz*pow(nx,5)*(nz*ty*tz + ny*(-pow(tx,2) + pow(tz,2)))),
   42*tx*pow(nx,5)*(-2*ny*nz*tz + ty*pow(ny,2) - ty*pow(nz,2)) + 
    14*nx*tx*pow(nz,3)*(15*nz*ty*pow(ny,2) + 20*tz*pow(ny,3) - 18*ny*tz*pow(nz,2) - 3*ty*pow(nz,3)) + 
    140*nz*tx*pow(nx,3)*(-3*nz*ty*pow(ny,2) - 2*tz*pow(ny,3) + 4*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
    35*pow(nx,4)*(-6*nz*ty*tz*pow(ny,2) + 4*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       pow(ny,3)*(pow(tx,2) - pow(tz,2))) - 21*pow(nx,2)*pow(nz,2)*
     (-20*nz*ty*tz*pow(ny,2) + 6*ty*tz*pow(nz,3) - 5*ny*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 10*pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 
    pow(nz,4)*(-42*nz*ty*tz*pow(ny,2) + 8*ty*tz*pow(nz,3) - 7*ny*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) + 
       35*pow(ny,3)*(pow(tx,2) - pow(tz,2))) - 7*pow(nx,6)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2))),
   42*tx*pow(nx,5)*(2*ny*nz*ty + tz*pow(ny,2) - tz*pow(nz,2)) + 140*tx*pow(nx,3)*pow(nz,2)*(-2*ny*nz*ty - 3*tz*pow(ny,2) + tz*pow(nz,2)) + 
    42*nx*tx*(2*ny*nz*ty + 5*tz*pow(ny,2) - tz*pow(nz,2))*pow(nz,4) - 
    pow(nz,5)*(14*ny*nz*ty*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 21*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    21*pow(nx,2)*pow(nz,3)*(10*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 10*pow(ny,2)*(pow(tx,2) - pow(tz,2))) - 
    35*nz*pow(nx,4)*(6*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - pow(tz,2))) + 
    7*pow(nx,6)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2))),70*tx*pow(nx,4)*(-3*nz*tz*pow(ny,2) + ty*pow(ny,3) - 3*ny*ty*pow(nz,2) + tz*pow(nz,3)) - 
    42*nz*tx*pow(nx,2)*(10*nz*ty*pow(ny,3) + 5*tz*pow(ny,4) - 20*tz*pow(ny,2)*pow(nz,2) - 10*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
    2*tx*pow(nz,3)*(35*nz*ty*pow(ny,3) + 35*tz*pow(ny,4) - 63*tz*pow(ny,2)*pow(nz,2) - 21*ny*ty*pow(nz,3) + 6*tz*pow(nz,4)) + 
    35*pow(nx,3)*(-8*nz*ty*tz*pow(ny,3) + 16*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       6*pow(ny,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,4)*(pow(tx,2) - pow(tz,2))) - 
    7*nx*pow(nz,2)*(-40*nz*ty*tz*pow(ny,3) + 36*ny*ty*tz*pow(nz,3) - 15*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
       3*pow(nz,4)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 15*pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    21*pow(nx,5)*(-4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2))),
   2*(35*tx*pow(nx,4)*(3*nz*ty*pow(ny,2) + tz*pow(ny,3) - 3*ny*tz*pow(nz,2) - ty*pow(nz,3)) + 
      42*tx*pow(nx,2)*pow(nz,2)*(-5*nz*ty*pow(ny,2) - 5*tz*pow(ny,3) + 5*ny*tz*pow(nz,2) + ty*pow(nz,3)) + 
      tx*(21*nz*ty*pow(ny,2) + 35*tz*pow(ny,3) - 21*ny*tz*pow(nz,2) - 3*ty*pow(nz,3))*pow(nz,4) + 
      70*nz*pow(nx,3)*(-3*nz*ty*tz*pow(ny,2) + ty*tz*pow(nz,3) - ny*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         pow(ny,3)*(pow(tx,2) - pow(tz,2))) - 7*nx*pow(nz,3)*
       (-15*nz*ty*tz*pow(ny,2) + 3*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 10*pow(ny,3)*(pow(tx,2) - pow(tz,2)))\
       + 21*pow(nx,5)*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2)))),
   70*tx*pow(nx,3)*(-4*nz*tz*pow(ny,3) + ty*pow(ny,4) - 6*ty*pow(ny,2)*pow(nz,2) + 4*ny*tz*pow(nz,3) + ty*pow(nz,4)) - 
    14*nx*nz*tx*(15*nz*ty*pow(ny,4) + 6*tz*pow(ny,5) - 40*tz*pow(ny,3)*pow(nz,2) - 30*ty*pow(ny,2)*pow(nz,3) + 18*ny*tz*pow(nz,4) + 
       3*ty*pow(nz,5)) + pow(nz,2)*(70*nz*ty*tz*pow(ny,4) - 126*ty*tz*pow(ny,2)*pow(nz,3) + 12*ty*tz*pow(nz,5) + 
       35*pow(ny,3)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 21*ny*pow(nz,4)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 
       21*pow(ny,5)*(pow(tx,2) - pow(tz,2))) + 21*pow(nx,2)*
     (-10*nz*ty*tz*pow(ny,4) + 40*ty*tz*pow(ny,2)*pow(nz,3) - 6*ty*tz*pow(nz,5) + 5*ny*pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       10*pow(ny,3)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,5)*(pow(tx,2) - pow(tz,2))) + 
    35*pow(nx,4)*(-6*nz*ty*tz*pow(ny,2) + 2*ty*tz*pow(nz,3) + pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 3*ny*pow(nz,2)*(-pow(ty,2) + pow(tz,2))),
   70*tx*pow(nx,3)*(4*nz*ty*pow(ny,3) + tz*pow(ny,4) - 6*tz*pow(ny,2)*pow(nz,2) - 4*ny*ty*pow(nz,3) + tz*pow(nz,4)) - 
    14*nx*tx*pow(nz,2)*(20*nz*ty*pow(ny,3) + 15*tz*pow(ny,4) - 30*tz*pow(ny,2)*pow(nz,2) - 12*ny*ty*pow(nz,3) + 3*tz*pow(nz,4)) + 
    pow(nz,3)*(70*nz*ty*tz*pow(ny,3) - 42*ny*ty*tz*pow(nz,3) + 21*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 
       3*pow(nz,4)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 35*pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    21*nz*pow(nx,2)*(-20*nz*ty*tz*pow(ny,3) + 20*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       10*pow(ny,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 5*pow(ny,4)*(pow(tx,2) - pow(tz,2))) + 
    35*pow(nx,4)*(2*ty*tz*pow(ny,3) - 6*ny*ty*tz*pow(nz,2) + 3*nz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,3)*(-pow(ty,2) + pow(tz,2))),
   42*tx*pow(nx,2)*(-5*nz*tz*pow(ny,4) + ty*pow(ny,5) - 10*ty*pow(ny,3)*pow(nz,2) + 10*tz*pow(ny,2)*pow(nz,3) + 5*ny*ty*pow(nz,4) - tz*pow(nz,5)) + 
    2*nz*tx*(-21*nz*ty*pow(ny,5) - 7*tz*pow(ny,6) + 70*tz*pow(ny,4)*pow(nz,2) + 70*ty*pow(ny,3)*pow(nz,3) - 63*tz*pow(ny,2)*pow(nz,4) - 
       21*ny*ty*pow(nz,5) + 4*tz*pow(nz,6)) + 7*nx*(-12*nz*ty*tz*pow(ny,5) + 80*ty*tz*pow(ny,3)*pow(nz,3) - 36*ny*ty*tz*pow(nz,5) - 
       pow(nz,6)*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) + 15*pow(ny,2)*pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 
       15*pow(ny,4)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(ny,6)*(pow(tx,2) - pow(tz,2))) + 
    35*pow(nx,3)*(-8*nz*ty*tz*pow(ny,3) + 8*ny*ty*tz*pow(nz,3) + pow(ny,4)*(pow(ty,2) - pow(tz,2)) + pow(nz,4)*(pow(ty,2) - pow(tz,2)) + 
       6*pow(ny,2)*pow(nz,2)*(-pow(ty,2) + pow(tz,2))),2*(21*tx*pow(nx,2)*
       (5*nz*ty*pow(ny,4) + tz*pow(ny,5) - 10*tz*pow(ny,3)*pow(nz,2) - 10*ty*pow(ny,2)*pow(nz,3) + 5*ny*tz*pow(nz,4) + ty*pow(nz,5)) - 
      tx*pow(nz,2)*(35*nz*ty*pow(ny,4) + 21*tz*pow(ny,5) - 70*tz*pow(ny,3)*pow(nz,2) - 42*ty*pow(ny,2)*pow(nz,3) + 21*ny*tz*pow(nz,4) + 
         3*ty*pow(nz,5)) + 7*nx*nz*(-15*nz*ty*tz*pow(ny,4) + 30*ty*tz*pow(ny,2)*pow(nz,3) - 3*ty*tz*pow(nz,5) + 
         3*ny*pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) - 10*pow(ny,3)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         3*pow(ny,5)*(pow(tx,2) - pow(tz,2))) + 35*pow(nx,3)*
       (ty*tz*pow(ny,4) - 6*ty*tz*pow(ny,2)*pow(nz,2) + ty*tz*pow(nz,4) + 2*nz*pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 
         2*ny*pow(nz,3)*(-pow(ty,2) + pow(tz,2)))),14*ty*(nx*tx - nz*tz)*pow(ny,6) + 
    70*nz*ty*pow(ny,4)*(-3*nx*nz*tx - 3*tz*pow(nx,2) + 2*tz*pow(nz,2)) + 42*ty*pow(ny,2)*(5*nx*nz*tx + 10*tz*pow(nx,2) - 3*tz*pow(nz,2))*pow(nz,3) + 
    2*ty*(-7*nx*nz*tx - 21*tz*pow(nx,2) + 4*tz*pow(nz,2))*pow(nz,5) + pow(ny,7)*(pow(tx,2) - pow(tz,2)) - 
    21*pow(ny,5)*(4*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    35*pow(ny,3)*pow(nz,2)*(8*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 6*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) - 
    7*ny*pow(nz,4)*(12*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) + 15*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   42*ty*pow(ny,5)*(2*nx*nz*tx + tz*pow(nx,2) - tz*pow(nz,2)) + 140*ty*pow(ny,3)*pow(nz,2)*(-2*nx*nz*tx - 3*tz*pow(nx,2) + tz*pow(nz,2)) + 
    42*ny*ty*(2*nx*nz*tx + 5*tz*pow(nx,2) - tz*pow(nz,2))*pow(nz,4) + 7*pow(ny,6)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
    35*nz*pow(ny,4)*(6*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 3*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
    21*pow(ny,2)*pow(nz,3)*(10*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 10*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) - 
    pow(nz,5)*(14*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) + 21*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),
   -42*nz*ty*(nz*tx + 2*nx*tz)*pow(ny,5) + 2*tx*ty*pow(ny,7) + 70*ty*(nz*tx + 4*nx*tz)*pow(ny,3)*pow(nz,3) - 14*ny*ty*(nz*tx + 6*nx*tz)*pow(nz,5) + 
    35*pow(ny,4)*pow(nz,2)*(2*nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 21*pow(ny,2)*pow(nz,4)*(2*nz*tx*tz - 5*nx*pow(ty,2) + 5*nx*pow(tz,2)) + 
    pow(nz,6)*(2*nz*tx*tz - 7*nx*pow(ty,2) + 7*nx*pow(tz,2)) - 7*pow(ny,6)*(2*nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2))),
   2*(7*ty*(nz*tx + nx*tz)*pow(ny,6) + tx*tz*pow(ny,7) - 35*ty*(nz*tx + 3*nx*tz)*pow(ny,4)*pow(nz,2) + 21*ty*(nz*tx + 5*nx*tz)*pow(ny,2)*pow(nz,4) - 
      ty*(nz*tx + 7*nx*tz)*pow(nz,6) + 35*pow(ny,3)*pow(nz,3)*(nz*tx*tz - 2*nx*pow(ty,2) + 2*nx*pow(tz,2)) - 
      7*ny*pow(nz,5)*(nz*tx*tz - 3*nx*pow(ty,2) + 3*nx*pow(tz,2)) - 21*nz*pow(ny,5)*(nz*tx*tz + nx*(-pow(ty,2) + pow(tz,2)))),
   -14*nz*ty*tz*pow(ny,6) + 70*ty*tz*pow(ny,4)*pow(nz,3) - 42*ty*tz*pow(ny,2)*pow(nz,5) + 2*ty*tz*pow(nz,7) + pow(ny,7)*(pow(ty,2) - pow(tz,2)) - 
    21*pow(ny,5)*pow(nz,2)*(pow(ty,2) - pow(tz,2)) + 35*pow(ny,3)*pow(nz,4)*(pow(ty,2) - pow(tz,2)) + 7*ny*pow(nz,6)*(-pow(ty,2) + pow(tz,2)),
   2*ty*tz*pow(ny,7) - 42*ty*tz*pow(ny,5)*pow(nz,2) + 70*ty*tz*pow(ny,3)*pow(nz,4) - 14*ny*ty*tz*pow(nz,6) + 
    7*nz*pow(ny,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(ny,4)*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + 21*pow(ny,2)*pow(nz,5)*(pow(ty,2) - pow(tz,2)) + 
    pow(nz,7)*(-pow(ty,2) + pow(tz,2)),-7*nz*(rz*tx + rx*tz)*pow(nx,6) + (rx*tx - rz*tz)*pow(nx,7) - 21*(rx*tx - rz*tz)*pow(nx,5)*pow(nz,2) + 
    35*(rz*tx + rx*tz)*pow(nx,4)*pow(nz,3) + 35*(rx*tx - rz*tz)*pow(nx,3)*pow(nz,4) - 21*(rz*tx + rx*tz)*pow(nx,2)*pow(nz,5) + 
    7*nx*(-(rx*tx) + rz*tz)*pow(nz,6) + (rz*tx + rx*tz)*pow(nz,7),
   -21*nz*(nz*ry*tx + 2*ny*rz*tx + nz*rx*ty + 2*ny*rx*tz)*pow(nx,5) - 7*(-(ny*rx*tx) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,6) + 
    (ry*tx + rx*ty)*pow(nx,7) + 35*(-3*ny*rx*tx + nz*rz*ty + nz*ry*tz + 3*ny*rz*tz)*pow(nx,4)*pow(nz,2) + 
    35*(nz*ry*tx + 4*ny*rz*tx + nz*rx*ty + 4*ny*rx*tz)*pow(nx,3)*pow(nz,3) - 
    21*(-5*ny*rx*tx + nz*rz*ty + nz*ry*tz + 5*ny*rz*tz)*pow(nx,2)*pow(nz,4) - 7*nx*(nz*ry*tx + 6*ny*rz*tx + nz*rx*ty + 6*ny*rx*tz)*pow(nz,5) + 
    (-7*ny*rx*tx + nz*rz*ty + nz*ry*tz + 7*ny*rz*tz)*pow(nz,6),
   7*nz*(rx*tx - rz*tz)*pow(nx,6) + (rz*tx + rx*tz)*pow(nx,7) - 21*(rz*tx + rx*tz)*pow(nx,5)*pow(nz,2) - 35*(rx*tx - rz*tz)*pow(nx,4)*pow(nz,3) + 
    35*(rz*tx + rx*tz)*pow(nx,3)*pow(nz,4) + 21*(rx*tx - rz*tz)*pow(nx,2)*pow(nz,5) - 7*nx*(rz*tx + rx*tz)*pow(nz,6) + (-(rx*tx) + rz*tz)*pow(nz,7),
   7*(ny*(ry*tx + rx*ty) - nz*(rz*tx + rx*tz))*pow(nx,6) + (ry*ty - rz*tz)*pow(nx,7) - 
    35*nz*pow(nx,4)*(3*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(nz,2)) + 
    35*pow(nx,3)*pow(nz,2)*(4*ny*nz*(rz*ty + ry*tz) + (-6*rx*tx + 6*rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2)) + 
    21*pow(nx,5)*(-2*ny*nz*(rz*ty + ry*tz) + (rx*tx - rz*tz)*pow(ny,2) - (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    21*pow(nx,2)*(5*ny*nz*(ry*tx + rx*ty) + 10*(rz*tx + rx*tz)*pow(ny,2) - 3*(rz*tx + rx*tz)*pow(nz,2))*pow(nz,3) - 
    7*nx*(6*ny*nz*(rz*ty + ry*tz) - 15*(rx*tx - rz*tz)*pow(ny,2) + (3*rx*tx + ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,4) + 
    (-7*ny*nz*(ry*tx + rx*ty) - 21*(rz*tx + rx*tz)*pow(ny,2) + 4*(rz*tx + rx*tz)*pow(nz,2))*pow(nz,5),
   -21*nz*(-2*ny*rx*tx + nz*rz*ty + nz*ry*tz + 2*ny*rz*tz)*pow(nx,5) + 7*(nz*ry*tx + ny*rz*tx + nz*rx*ty + ny*rx*tz)*pow(nx,6) + 
    (rz*ty + ry*tz)*pow(nx,7) - 35*(nz*ry*tx + 3*ny*rz*tx + nz*rx*ty + 3*ny*rx*tz)*pow(nx,4)*pow(nz,2) + 
    35*(-4*ny*rx*tx + nz*rz*ty + nz*ry*tz + 4*ny*rz*tz)*pow(nx,3)*pow(nz,3) + 
    21*(nz*ry*tx + 5*ny*rz*tx + nz*rx*ty + 5*ny*rx*tz)*pow(nx,2)*pow(nz,4) - 7*nx*(-6*ny*rx*tx + nz*rz*ty + nz*ry*tz + 6*ny*rz*tz)*pow(nz,5) - 
    (nz*ry*tx + 7*ny*rz*tx + nz*rx*ty + 7*ny*rx*tz)*pow(nz,6),
   -7*(-(ny*ry*ty) + nz*rz*ty + nz*ry*tz + ny*rz*tz)*pow(nx,6) + 
    21*pow(nx,5)*(-2*ny*nz*(rz*tx + rx*tz) + (ry*tx + rx*ty)*pow(ny,2) - (ry*tx + rx*ty)*pow(nz,2)) + 
    7*nx*pow(nz,3)*(15*nz*(ry*tx + rx*ty)*pow(ny,2) + 20*(rz*tx + rx*tz)*pow(ny,3) - 18*ny*(rz*tx + rx*tz)*pow(nz,2) - 
       3*(ry*tx + rx*ty)*pow(nz,3)) + 70*nz*pow(nx,3)*(-3*nz*(ry*tx + rx*ty)*pow(ny,2) - 2*(rz*tx + rx*tz)*pow(ny,3) + 
       4*ny*(rz*tx + rx*tz)*pow(nz,2) + (ry*tx + rx*ty)*pow(nz,3)) + 
    35*pow(nx,4)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (rx*tx - rz*tz)*pow(ny,3) - 3*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + 
       2*(rz*ty + ry*tz)*pow(nz,3)) - 21*pow(nx,2)*pow(nz,2)*
     (-10*nz*(rz*ty + ry*tz)*pow(ny,2) + 10*(rx*tx - rz*tz)*pow(ny,3) - 5*ny*(2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2) + 3*(rz*ty + ry*tz)*pow(nz,3)) + 
    (-21*nz*(rz*ty + ry*tz)*pow(ny,2) + 35*(rx*tx - rz*tz)*pow(ny,3) - 7*ny*(3*rx*tx + ry*ty - 4*rz*tz)*pow(nz,2) + 4*(rz*ty + ry*tz)*pow(nz,3))*
     pow(nz,4),7*(nz*ry*ty + ny*rz*ty + ny*ry*tz - nz*rz*tz)*pow(nx,6) + 
    21*pow(nx,5)*(2*ny*nz*(ry*tx + rx*ty) + (rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    70*pow(nx,3)*pow(nz,2)*(2*ny*nz*(ry*tx + rx*ty) + 3*(rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2)) - 
    35*nz*pow(nx,4)*(3*ny*nz*(rz*ty + ry*tz) + (-3*rx*tx + 3*rz*tz)*pow(ny,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    21*pow(nx,2)*(5*ny*nz*(rz*ty + ry*tz) - 10*(rx*tx - rz*tz)*pow(ny,2) + (2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3) + 
    21*nx*(2*ny*nz*(ry*tx + rx*ty) + 5*(rz*tx + rx*tz)*pow(ny,2) - (rz*tx + rx*tz)*pow(nz,2))*pow(nz,4) - 
    (7*ny*nz*(rz*ty + ry*tz) - 21*(rx*tx - rz*tz)*pow(ny,2) + (3*rx*tx + ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,5),
   21*pow(nx,5)*(-2*ny*nz*(rz*ty + ry*tz) + (ry*ty - rz*tz)*pow(ny,2) + (-(ry*ty) + rz*tz)*pow(nz,2)) + 
    35*pow(nx,4)*(-3*nz*(rz*tx + rx*tz)*pow(ny,2) + (ry*tx + rx*ty)*pow(ny,3) - 3*ny*(ry*tx + rx*ty)*pow(nz,2) + (rz*tx + rx*tz)*pow(nz,3)) - 
    21*nz*pow(nx,2)*(10*nz*(ry*tx + rx*ty)*pow(ny,3) + 5*(rz*tx + rx*tz)*pow(ny,4) - 20*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       10*ny*(ry*tx + rx*ty)*pow(nz,3) + 3*(rz*tx + rx*tz)*pow(nz,4)) + 
    pow(nz,3)*(35*nz*(ry*tx + rx*ty)*pow(ny,3) + 35*(rz*tx + rx*tz)*pow(ny,4) - 63*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       21*ny*(ry*tx + rx*ty)*pow(nz,3) + 6*(rz*tx + rx*tz)*pow(nz,4)) + 
    35*pow(nx,3)*(-4*nz*(rz*ty + ry*tz)*pow(ny,3) + (rx*tx - rz*tz)*pow(ny,4) - 6*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,2)*pow(nz,2) + 
       8*ny*(rz*ty + ry*tz)*pow(nz,3) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4)) - 
    7*nx*pow(nz,2)*(-20*nz*(rz*ty + ry*tz)*pow(ny,3) + 15*(rx*tx - rz*tz)*pow(ny,4) - 15*(2*rx*tx + ry*ty - 3*rz*tz)*pow(ny,2)*pow(nz,2) + 
       18*ny*(rz*ty + ry*tz)*pow(nz,3) + 3*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,4)),
   21*pow(nx,5)*(2*ny*nz*(ry*ty - rz*tz) + (rz*ty + ry*tz)*pow(ny,2) - (rz*ty + ry*tz)*pow(nz,2)) + 
    35*pow(nx,4)*(3*nz*(ry*tx + rx*ty)*pow(ny,2) + (rz*tx + rx*tz)*pow(ny,3) - 3*ny*(rz*tx + rx*tz)*pow(nz,2) - (ry*tx + rx*ty)*pow(nz,3)) + 
    42*pow(nx,2)*pow(nz,2)*(-5*nz*(ry*tx + rx*ty)*pow(ny,2) - 5*(rz*tx + rx*tz)*pow(ny,3) + 5*ny*(rz*tx + rx*tz)*pow(nz,2) + 
       (ry*tx + rx*ty)*pow(nz,3)) + 70*nz*pow(nx,3)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + 2*(rx*tx - rz*tz)*pow(ny,3) - 
       2*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3)) - 
    7*nx*pow(nz,3)*(-15*nz*(rz*ty + ry*tz)*pow(ny,2) + 20*(rx*tx - rz*tz)*pow(ny,3) - 6*ny*(2*rx*tx + ry*ty - 3*rz*tz)*pow(nz,2) + 
       3*(rz*ty + ry*tz)*pow(nz,3)) + (21*nz*(ry*tx + rx*ty)*pow(ny,2) + 35*(rz*tx + rx*tz)*pow(ny,3) - 21*ny*(rz*tx + rx*tz)*pow(nz,2) - 
       3*(ry*tx + rx*ty)*pow(nz,3))*pow(nz,4),35*pow(nx,4)*(-3*nz*(rz*ty + ry*tz)*pow(ny,2) + (ry*ty - rz*tz)*pow(ny,3) + 
       3*ny*(-(ry*ty) + rz*tz)*pow(nz,2) + (rz*ty + ry*tz)*pow(nz,3)) + 
    35*pow(nx,3)*(-4*nz*(rz*tx + rx*tz)*pow(ny,3) + (ry*tx + rx*ty)*pow(ny,4) - 6*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,2) + 
       4*ny*(rz*tx + rx*tz)*pow(nz,3) + (ry*tx + rx*ty)*pow(nz,4)) - 
    7*nx*nz*(15*nz*(ry*tx + rx*ty)*pow(ny,4) + 6*(rz*tx + rx*tz)*pow(ny,5) - 40*(rz*tx + rx*tz)*pow(ny,3)*pow(nz,2) - 
       30*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,3) + 18*ny*(rz*tx + rx*tz)*pow(nz,4) + 3*(ry*tx + rx*ty)*pow(nz,5)) + 
    21*pow(nx,2)*(-5*nz*(rz*ty + ry*tz)*pow(ny,4) + (rx*tx - rz*tz)*pow(ny,5) - 10*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,3)*pow(nz,2) + 
       20*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,3) + 5*ny*(rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4) - 3*(rz*ty + ry*tz)*pow(nz,5)) + 
    pow(nz,2)*(35*nz*(rz*ty + ry*tz)*pow(ny,4) - 21*(rx*tx - rz*tz)*pow(ny,5) + 35*(2*rx*tx + ry*ty - 3*rz*tz)*pow(ny,3)*pow(nz,2) - 
       63*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,3) - 21*ny*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,4) + 6*(rz*ty + ry*tz)*pow(nz,5)),
   35*pow(nx,4)*(3*nz*(ry*ty - rz*tz)*pow(ny,2) + (rz*ty + ry*tz)*pow(ny,3) - 3*ny*(rz*ty + ry*tz)*pow(nz,2) + (-(ry*ty) + rz*tz)*pow(nz,3)) + 
    35*pow(nx,3)*(4*nz*(ry*tx + rx*ty)*pow(ny,3) + (rz*tx + rx*tz)*pow(ny,4) - 6*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       4*ny*(ry*tx + rx*ty)*pow(nz,3) + (rz*tx + rx*tz)*pow(nz,4)) - 
    7*nx*pow(nz,2)*(20*nz*(ry*tx + rx*ty)*pow(ny,3) + 15*(rz*tx + rx*tz)*pow(ny,4) - 30*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,2) - 
       12*ny*(ry*tx + rx*ty)*pow(nz,3) + 3*(rz*tx + rx*tz)*pow(nz,4)) + 
    21*nz*pow(nx,2)*(-10*nz*(rz*ty + ry*tz)*pow(ny,3) + 5*(rx*tx - rz*tz)*pow(ny,4) - 10*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,2)*pow(nz,2) + 
       10*ny*(rz*ty + ry*tz)*pow(nz,3) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4)) + 
    pow(nz,3)*(35*nz*(rz*ty + ry*tz)*pow(ny,3) - 35*(rx*tx - rz*tz)*pow(ny,4) + 21*(2*rx*tx + ry*ty - 3*rz*tz)*pow(ny,2)*pow(nz,2) - 
       21*ny*(rz*ty + ry*tz)*pow(nz,3) - 3*(rx*tx + ry*ty - 2*rz*tz)*pow(nz,4)),
   35*pow(nx,3)*(-4*nz*(rz*ty + ry*tz)*pow(ny,3) + (ry*ty - rz*tz)*pow(ny,4) + 6*(-(ry*ty) + rz*tz)*pow(ny,2)*pow(nz,2) + 
       4*ny*(rz*ty + ry*tz)*pow(nz,3) + (ry*ty - rz*tz)*pow(nz,4)) + 
    21*pow(nx,2)*(-5*nz*(rz*tx + rx*tz)*pow(ny,4) + (ry*tx + rx*ty)*pow(ny,5) - 10*(ry*tx + rx*ty)*pow(ny,3)*pow(nz,2) + 
       10*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,3) + 5*ny*(ry*tx + rx*ty)*pow(nz,4) - (rz*tx + rx*tz)*pow(nz,5)) + 
    nz*(-21*nz*(ry*tx + rx*ty)*pow(ny,5) - 7*(rz*tx + rx*tz)*pow(ny,6) + 70*(rz*tx + rx*tz)*pow(ny,4)*pow(nz,2) + 
       70*(ry*tx + rx*ty)*pow(ny,3)*pow(nz,3) - 63*(rz*tx + rx*tz)*pow(ny,2)*pow(nz,4) - 21*ny*(ry*tx + rx*ty)*pow(nz,5) + 
       4*(rz*tx + rx*tz)*pow(nz,6)) + 7*nx*(-6*nz*(rz*ty + ry*tz)*pow(ny,5) + (rx*tx - rz*tz)*pow(ny,6) - 
       15*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,4)*pow(nz,2) + 40*(rz*ty + ry*tz)*pow(ny,3)*pow(nz,3) + 
       15*(rx*tx + 2*ry*ty - 3*rz*tz)*pow(ny,2)*pow(nz,4) - 18*ny*(rz*ty + ry*tz)*pow(nz,5) - (rx*tx + 3*ry*ty - 4*rz*tz)*pow(nz,6)),
   35*pow(nx,3)*(4*nz*(ry*ty - rz*tz)*pow(ny,3) + (rz*ty + ry*tz)*pow(ny,4) - 6*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,2) + 
       4*ny*(-(ry*ty) + rz*tz)*pow(nz,3) + (rz*ty + ry*tz)*pow(nz,4)) + 
    21*pow(nx,2)*(5*nz*(ry*tx + rx*ty)*pow(ny,4) + (rz*tx + rx*tz)*pow(ny,5) - 10*(rz*tx + rx*tz)*pow(ny,3)*pow(nz,2) - 
       10*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,3) + 5*ny*(rz*tx + rx*tz)*pow(nz,4) + (ry*tx + rx*ty)*pow(nz,5)) - 
    pow(nz,2)*(35*nz*(ry*tx + rx*ty)*pow(ny,4) + 21*(rz*tx + rx*tz)*pow(ny,5) - 70*(rz*tx + rx*tz)*pow(ny,3)*pow(nz,2) - 
       42*(ry*tx + rx*ty)*pow(ny,2)*pow(nz,3) + 21*ny*(rz*tx + rx*tz)*pow(nz,4) + 3*(ry*tx + rx*ty)*pow(nz,5)) - 
    7*nx*nz*(15*nz*(rz*ty + ry*tz)*pow(ny,4) + (-6*rx*tx + 6*rz*tz)*pow(ny,5) + 20*(rx*tx + ry*ty - 2*rz*tz)*pow(ny,3)*pow(nz,2) - 
       30*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,3) - 6*ny*(rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,4) + 3*(rz*ty + ry*tz)*pow(nz,5)),
   7*(nx*(ry*tx + rx*ty) - nz*(rz*ty + ry*tz))*pow(ny,6) + (rx*tx - rz*tz)*pow(ny,7) - 
    35*nz*pow(ny,4)*(3*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - 2*(rz*ty + ry*tz)*pow(nz,2)) + 
    35*pow(ny,3)*pow(nz,2)*(4*nx*nz*(rz*tx + rx*tz) + 6*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2)) - 
    21*pow(ny,5)*(2*nx*nz*(rz*tx + rx*tz) + (-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    21*pow(ny,2)*(5*nx*nz*(ry*tx + rx*ty) + 10*(rz*ty + ry*tz)*pow(nx,2) - 3*(rz*ty + ry*tz)*pow(nz,2))*pow(nz,3) - 
    7*ny*(6*nx*nz*(rz*tx + rx*tz) + 15*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 3*ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,4) + 
    (-7*nx*nz*(ry*tx + rx*ty) - 21*(rz*ty + ry*tz)*pow(nx,2) + 4*(rz*ty + ry*tz)*pow(nz,2))*pow(nz,5),
   7*(nz*rx*tx + nx*rz*tx + nx*rx*tz - nz*rz*tz)*pow(ny,6) + 
    21*pow(ny,5)*(2*nx*nz*(ry*tx + rx*ty) + (rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    70*pow(ny,3)*pow(nz,2)*(2*nx*nz*(ry*tx + rx*ty) + 3*(rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2)) - 
    35*nz*pow(ny,4)*(3*nx*nz*(rz*tx + rx*tz) + 3*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + ry*ty - 2*rz*tz)*pow(nz,2)) + 
    21*pow(ny,2)*(5*nx*nz*(rz*tx + rx*tz) + 10*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 2*ry*ty - 3*rz*tz)*pow(nz,2))*pow(nz,3) + 
    21*ny*(2*nx*nz*(ry*tx + rx*ty) + 5*(rz*ty + ry*tz)*pow(nx,2) - (rz*ty + ry*tz)*pow(nz,2))*pow(nz,4) - 
    (7*nx*nz*(rz*tx + rx*tz) + 21*(-(ry*ty) + rz*tz)*pow(nx,2) + (rx*tx + 3*ry*ty - 4*rz*tz)*pow(nz,2))*pow(nz,5),
   -21*nz*(nz*ry*tx + nz*rx*ty + 2*nx*rz*ty + 2*nx*ry*tz)*pow(ny,5) - 7*(nz*rz*tx - nx*ry*ty + nz*rx*tz + nx*rz*tz)*pow(ny,6) + 
    (ry*tx + rx*ty)*pow(ny,7) + 35*(nz*rz*tx - 3*nx*ry*ty + nz*rx*tz + 3*nx*rz*tz)*pow(ny,4)*pow(nz,2) + 
    35*(nz*ry*tx + nz*rx*ty + 4*nx*rz*ty + 4*nx*ry*tz)*pow(ny,3)*pow(nz,3) - 
    21*(nz*rz*tx - 5*nx*ry*ty + nz*rx*tz + 5*nx*rz*tz)*pow(ny,2)*pow(nz,4) - 7*ny*(nz*ry*tx + nz*rx*ty + 6*nx*rz*ty + 6*nx*ry*tz)*pow(nz,5) + 
    (nz*rz*tx - 7*nx*ry*ty + nz*rx*tz + 7*nx*rz*tz)*pow(nz,6),
   -21*nz*(nz*rz*tx - 2*nx*ry*ty + nz*rx*tz + 2*nx*rz*tz)*pow(ny,5) + 7*(nz*ry*tx + nz*rx*ty + nx*rz*ty + nx*ry*tz)*pow(ny,6) + 
    (rz*tx + rx*tz)*pow(ny,7) - 35*(nz*ry*tx + nz*rx*ty + 3*nx*rz*ty + 3*nx*ry*tz)*pow(ny,4)*pow(nz,2) + 
    35*(nz*rz*tx - 4*nx*ry*ty + nz*rx*tz + 4*nx*rz*tz)*pow(ny,3)*pow(nz,3) + 
    21*(nz*ry*tx + nz*rx*ty + 5*nx*rz*ty + 5*nx*ry*tz)*pow(ny,2)*pow(nz,4) - 7*ny*(nz*rz*tx - 6*nx*ry*ty + nz*rx*tz + 6*nx*rz*tz)*pow(nz,5) - 
    (nz*ry*tx + nz*rx*ty + 7*nx*rz*ty + 7*nx*ry*tz)*pow(nz,6),
   -7*nz*(rz*ty + ry*tz)*pow(ny,6) + (ry*ty - rz*tz)*pow(ny,7) - 21*(ry*ty - rz*tz)*pow(ny,5)*pow(nz,2) + 35*(rz*ty + ry*tz)*pow(ny,4)*pow(nz,3) + 
    35*(ry*ty - rz*tz)*pow(ny,3)*pow(nz,4) - 21*(rz*ty + ry*tz)*pow(ny,2)*pow(nz,5) + 7*ny*(-(ry*ty) + rz*tz)*pow(nz,6) + (rz*ty + ry*tz)*pow(nz,7),
   7*nz*(ry*ty - rz*tz)*pow(ny,6) + (rz*ty + ry*tz)*pow(ny,7) - 21*(rz*ty + ry*tz)*pow(ny,5)*pow(nz,2) - 35*(ry*ty - rz*tz)*pow(ny,4)*pow(nz,3) + 
    35*(rz*ty + ry*tz)*pow(ny,3)*pow(nz,4) + 21*(ry*ty - rz*tz)*pow(ny,2)*pow(nz,5) - 7*ny*(rz*ty + ry*tz)*pow(nz,6) + (-(ry*ty) + rz*tz)*pow(nz,7),
   -15*tx*pow(nx,4)*pow(nz,2)*(pow(tx,2) - 3*pow(tz,2)) + 15*tx*pow(nx,2)*pow(nz,4)*(pow(tx,2) - 3*pow(tz,2)) - 
    tx*pow(nz,6)*(pow(tx,2) - 3*pow(tz,2)) + 6*nz*tz*pow(nx,5)*(-3*pow(tx,2) + pow(tz,2)) - 20*tz*pow(nx,3)*pow(nz,3)*(-3*pow(tx,2) + pow(tz,2)) + 
    6*nx*tz*pow(nz,5)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,6)*(pow(tx,3) - 3*tx*pow(tz,2)),
   3*(2*tx*pow(nx,5)*(-6*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 2*nx*tx*pow(nz,4)*(-6*nz*ty*tz + 5*ny*(pow(tx,2) - 3*pow(tz,2))) + 
      ty*pow(nx,6)*(pow(tx,2) - pow(tz,2)) + 20*tx*pow(nx,3)*pow(nz,2)*(2*nz*ty*tz - ny*pow(tx,2) + 3*ny*pow(tz,2)) + 
      5*pow(nx,2)*pow(nz,3)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - 4*ny*tz*(-3*pow(tx,2) + pow(tz,2))) + 
      pow(nz,5)*(2*ny*tz*(-3*pow(tx,2) + pow(tz,2)) + nz*ty*(-pow(tx,2) + pow(tz,2))) - 
      5*nz*pow(nx,4)*(6*ny*tz*pow(tx,2) + 3*nz*ty*(pow(tx,2) - pow(tz,2)) - 2*ny*pow(tz,3))),
   6*nz*tx*pow(nx,5)*(pow(tx,2) - 3*pow(tz,2)) - 20*tx*pow(nx,3)*pow(nz,3)*(pow(tx,2) - 3*pow(tz,2)) + 6*nx*tx*pow(nz,5)*(pow(tx,2) - 3*pow(tz,2)) + 
    15*tz*pow(nx,4)*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) - 15*tz*pow(nx,2)*pow(nz,4)*(-3*pow(tx,2) + pow(tz,2)) + 
    tz*pow(nz,6)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,6)*(3*tz*pow(tx,2) - pow(tz,3)),
   3*(5*tx*pow(nx,2)*pow(nz,2)*(24*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 6*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) - 
      tx*pow(nz,4)*(12*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 5*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
      5*tx*pow(nx,4)*(-12*ny*nz*ty*tz - pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
      tx*pow(nx,6)*(pow(ty,2) - pow(tz,2)) + pow(nx,5)*(6*ny*ty*(pow(tx,2) - pow(tz,2)) + 2*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      2*nx*pow(nz,3)*(15*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2)) + 
         10*pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 20*nz*pow(nx,3)*
       (tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 3*ny*nz*ty*(-pow(tx,2) + pow(tz,2)) + pow(ny,2)*(-3*tz*pow(tx,2) + pow(tz,3)))),
   6*(tx*ty*tz*pow(nx,6) + 5*nz*tx*pow(nx,4)*(-3*nz*ty*tz + ny*(pow(tx,2) - 3*pow(tz,2))) + 
      tx*pow(nz,5)*(-(nz*ty*tz) + ny*(pow(tx,2) - 3*pow(tz,2))) + 5*tx*pow(nx,2)*pow(nz,3)*(3*nz*ty*tz - 2*ny*pow(tx,2) + 6*ny*pow(tz,2)) + 
      nx*pow(nz,4)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - 5*ny*tz*(-3*pow(tx,2) + pow(tz,2))) - 
      10*pow(nx,3)*pow(nz,2)*(nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2))) + 
      pow(nx,5)*(3*nz*ty*(pow(tx,2) - pow(tz,2)) - ny*tz*(-3*pow(tx,2) + pow(tz,2)))),
   20*tx*pow(nx,3)*(-18*nz*ty*tz*pow(ny,2) + 12*ty*tz*pow(nz,3) - 3*ny*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
       pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) - 6*nx*tx*pow(nz,2)*
     (-60*nz*ty*tz*pow(ny,2) + 18*ty*tz*pow(nz,3) - 5*ny*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) + 
       10*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + pow(nx,6)*(pow(ty,3) - 3*ty*pow(tz,2)) + 
    15*nz*pow(nx,2)*(ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 18*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       12*ny*tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 4*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) - 
    pow(nz,3)*(ty*pow(nz,3)*(9*pow(tx,2) + pow(ty,2) - 12*pow(tz,2)) + 6*ny*tz*pow(nz,2)*(9*pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) - 
       45*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 20*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) - 
    18*tx*pow(nx,5)*(2*nz*ty*tz + ny*(-pow(ty,2) + pow(tz,2))) + 
    15*pow(nx,4)*(-(ty*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
       2*ny*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))),
   6*nx*tx*pow(nz,3)*(30*ny*nz*ty*tz + pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 10*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) - 
    20*nz*tx*pow(nx,3)*(18*ny*nz*ty*tz + pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) - 3*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2))) + 
    18*tx*pow(nx,5)*(2*ny*ty*tz + nz*(pow(ty,2) - pow(tz,2))) - 
    45*pow(nx,2)*pow(nz,2)*(4*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)) + 
       pow(ny,2)*(6*tz*pow(tx,2) - 2*pow(tz,3))) + 15*pow(nx,4)*
     (6*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
    pow(nz,4)*(18*ny*nz*ty*(pow(tx,2) - pow(tz,2)) + tz*pow(nz,2)*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2)) + 
       15*pow(ny,2)*(3*tz*pow(tx,2) - pow(tz,3))) + pow(nx,6)*(3*tz*pow(ty,2) - pow(tz,3)),
   3*(tx*pow(nz,2)*(40*nz*ty*tz*pow(ny,3) - 36*ny*ty*tz*pow(nz,3) + 5*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 
         pow(nz,4)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) - 5*pow(ny,4)*(pow(tx,2) - 3*pow(tz,2))) + 
      5*tx*pow(nx,2)*(-24*nz*ty*tz*pow(ny,3) + 48*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 
         6*pow(ny,2)*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,4)*(pow(tx,2) - 3*pow(tz,2))) + 
      15*tx*pow(nx,4)*(-4*ny*nz*ty*tz + pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,2)*(-pow(ty,2) + pow(tz,2))) + 
      20*pow(nx,3)*(-(ny*ty*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 
         tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + nz*tz*pow(ny,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      2*pow(nx,5)*(nz*tz*(-3*pow(ty,2) + pow(tz,2)) + ny*(pow(ty,3) - 3*ty*pow(tz,2))) - 
      2*nx*nz*(-5*ny*ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + 3*tz*pow(nz,4)*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         30*nz*ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 30*tz*pow(ny,2)*pow(nz,2)*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)) + 
         5*pow(ny,4)*(3*tz*pow(tx,2) - pow(tz,3)))),2*(tx*pow(nz,3)*
       (45*nz*ty*tz*pow(ny,2) - 9*ty*tz*pow(nz,3) + 3*ny*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) - 
         10*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + 30*nz*tx*pow(nx,2)*
       (-9*nz*ty*tz*pow(ny,2) + 3*ty*tz*pow(nz,3) - ny*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + pow(ny,3)*(pow(tx,2) - 3*pow(tz,2))) + 
      45*tx*pow(nx,4)*(ty*tz*pow(ny,2) - ty*tz*pow(nz,2) + ny*nz*(pow(ty,2) - pow(tz,2))) + 
      3*nx*pow(nz,2)*(ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 30*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         15*ny*tz*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 10*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2))) + 
      10*pow(nx,3)*(-(ty*pow(nz,3)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 9*nz*ty*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         3*ny*tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,3)*(3*tz*pow(tx,2) - pow(tz,3))) + 
      3*pow(nx,5)*(3*ny*tz*pow(ty,2) + nz*pow(ty,3) - 3*nz*ty*pow(tz,2) - ny*pow(tz,3))),
   3*(2*nx*tx*(-30*nz*ty*tz*pow(ny,4) + 120*ty*tz*pow(ny,2)*pow(nz,3) - 18*ty*tz*pow(nz,5) + 
         5*ny*pow(nz,4)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 10*pow(ny,3)*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
         pow(ny,5)*(pow(tx,2) - 3*pow(tz,2))) + 20*tx*pow(nx,3)*
       (-6*nz*ty*tz*pow(ny,2) + 2*ty*tz*pow(nz,3) + pow(ny,3)*(pow(ty,2) - pow(tz,2)) + 3*ny*pow(nz,2)*(-pow(ty,2) + pow(tz,2))) + 
      5*pow(nx,2)*(ty*pow(nz,4)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*ty*pow(ny,2)*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
         3*ty*pow(ny,4)*(pow(tx,2) - pow(tz,2)) + 12*ny*tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 
         4*nz*tz*pow(ny,3)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      5*pow(nx,4)*(-(ty*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2))) + 2*ny*nz*tz*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,2)*(pow(ty,3) - 3*ty*pow(tz,2))) + 
      nz*(5*ty*pow(ny,2)*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - ty*pow(nz,5)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 
         15*nz*ty*pow(ny,4)*(pow(tx,2) - pow(tz,2)) + 20*tz*pow(ny,3)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
         6*ny*tz*pow(nz,4)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(ny,5)*(-6*tz*pow(tx,2) + 2*pow(tz,3)))),
   3*(2*nx*nz*tx*(-60*nz*ty*tz*pow(ny,3) + 60*ny*ty*tz*pow(nz,3) + pow(nz,4)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 
         10*pow(ny,2)*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 5*pow(ny,4)*(pow(tx,2) - 3*pow(tz,2))) + 
      20*tx*pow(nx,3)*(2*ty*tz*pow(ny,3) - 6*ny*ty*tz*pow(nz,2) + 3*nz*pow(ny,2)*(pow(ty,2) - pow(tz,2)) + pow(nz,3)*(-pow(ty,2) + pow(tz,2))) + 
      pow(nz,2)*(2*ny*ty*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 20*nz*ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 
         15*tz*pow(ny,2)*pow(nz,2)*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 5*tz*pow(ny,4)*(-3*pow(tx,2) + pow(tz,2)) + 
         tz*pow(nz,4)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      5*pow(nx,2)*(-4*ny*ty*pow(nz,3)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 12*nz*ty*pow(ny,3)*(pow(tx,2) - pow(tz,2)) + 
         3*tz*pow(nz,4)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 6*tz*pow(ny,2)*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + 
         pow(ny,4)*(3*tz*pow(tx,2) - pow(tz,3))) + 5*pow(nx,4)*
       (2*ny*nz*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,2)*(3*tz*pow(ty,2) - pow(tz,3)))),
   20*ty*pow(ny,3)*(-18*nz*tx*tz*pow(nx,2) + 12*tx*tz*pow(nz,3) - 3*nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
       pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) - 6*ny*ty*pow(nz,2)*
     (-60*nz*tx*tz*pow(nx,2) + 18*tx*tz*pow(nz,3) - 5*nx*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 
       10*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + pow(ny,6)*(pow(tx,3) - 3*tx*pow(tz,2)) - 
    18*ty*pow(ny,5)*(2*nz*tx*tz + nx*(-pow(tx,2) + pow(tz,2))) + 
    15*nz*pow(ny,2)*(tx*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 18*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
       12*nx*tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 4*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
    pow(nz,3)*(tx*pow(nz,3)*(pow(tx,2) + 9*pow(ty,2) - 12*pow(tz,2)) + 6*nx*tz*pow(nz,2)*(3*pow(tx,2) + 9*pow(ty,2) - 4*pow(tz,2)) - 
       45*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 20*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
    15*pow(ny,4)*(tx*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 2*nx*nz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
       3*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))),2*(ty*pow(nz,3)*
       (45*nz*tx*tz*pow(nx,2) - 9*tx*tz*pow(nz,3) + 3*nx*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 
         10*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 30*nz*ty*pow(ny,2)*
       (-9*nz*tx*tz*pow(nx,2) + 3*tx*tz*pow(nz,3) - nx*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + pow(nx,3)*(pow(ty,2) - 3*pow(tz,2))) + 
      45*ty*pow(ny,4)*(tx*tz*pow(nx,2) - tx*tz*pow(nz,2) + nx*nz*(pow(tx,2) - pow(tz,2))) + 
      3*ny*pow(nz,2)*(tx*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) - 30*nz*tx*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         15*nx*tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 10*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2))) - 
      10*pow(ny,3)*(tx*pow(nz,3)*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*nx*tz*pow(nz,2)*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 9*nz*tx*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 
      3*pow(ny,5)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3))),
   3*(5*ty*pow(ny,2)*pow(nz,2)*(24*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
      ty*pow(nz,4)*(12*nx*nz*tx*tz + pow(nz,2)*(pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 5*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
      5*ty*pow(ny,4)*(12*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
      ty*pow(ny,6)*(pow(tx,2) - pow(tz,2)) + pow(ny,5)*(6*nx*tx*(pow(ty,2) - pow(tz,2)) + 2*nz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
      2*ny*pow(nz,3)*(15*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2)) + 
         10*pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + 20*nz*pow(ny,3)*
       (tz*pow(nz,2)*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 3*nx*nz*tx*(-pow(ty,2) + pow(tz,2)) + pow(nx,2)*(-3*tz*pow(ty,2) + pow(tz,3)))),
   6*ny*ty*pow(nz,3)*(30*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 10*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) - 
    20*nz*ty*pow(ny,3)*(18*nx*nz*tx*tz + pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 3*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2))) + 
    18*ty*pow(ny,5)*(2*nx*tx*tz + nz*(pow(tx,2) - pow(tz,2))) - 
    45*pow(ny,2)*pow(nz,2)*(4*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)) + 
       pow(nx,2)*(6*tz*pow(ty,2) - 2*pow(tz,3))) + 15*pow(ny,4)*
     (6*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)) + pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + 
    pow(nz,4)*(18*nx*nz*tx*(pow(ty,2) - pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2)) + 
       15*pow(nx,2)*(3*tz*pow(ty,2) - pow(tz,3))) + pow(ny,6)*(3*tz*pow(tx,2) - pow(tz,3)),
   3*(2*ty*pow(ny,5)*(-6*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 2*ny*ty*pow(nz,4)*(-6*nz*tx*tz + 5*nx*(pow(ty,2) - 3*pow(tz,2))) + 
      tx*pow(ny,6)*(pow(ty,2) - pow(tz,2)) + 20*ty*pow(ny,3)*pow(nz,2)*(2*nz*tx*tz - nx*pow(ty,2) + 3*nx*pow(tz,2)) + 
      5*pow(ny,2)*pow(nz,3)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - 4*nx*tz*(-3*pow(ty,2) + pow(tz,2))) + 
      pow(nz,5)*(2*nx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*tx*(-pow(ty,2) + pow(tz,2))) - 
      5*nz*pow(ny,4)*(6*nx*tz*pow(ty,2) + 3*nz*tx*(pow(ty,2) - pow(tz,2)) - 2*nx*pow(tz,3))),
   6*(tx*ty*tz*pow(ny,6) + 5*nz*ty*pow(ny,4)*(-3*nz*tx*tz + nx*(pow(ty,2) - 3*pow(tz,2))) + 
      ty*pow(nz,5)*(-(nz*tx*tz) + nx*(pow(ty,2) - 3*pow(tz,2))) + 5*ty*pow(ny,2)*pow(nz,3)*(3*nz*tx*tz - 2*nx*pow(ty,2) + 6*nx*pow(tz,2)) + 
      ny*pow(nz,4)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - 5*nx*tz*(-3*pow(ty,2) + pow(tz,2))) - 
      10*pow(ny,3)*pow(nz,2)*(nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2))) + 
      pow(ny,5)*(3*nz*tx*(pow(ty,2) - pow(tz,2)) - nx*tz*(-3*pow(ty,2) + pow(tz,2)))),
   -15*ty*pow(ny,4)*pow(nz,2)*(pow(ty,2) - 3*pow(tz,2)) + 15*ty*pow(ny,2)*pow(nz,4)*(pow(ty,2) - 3*pow(tz,2)) - 
    ty*pow(nz,6)*(pow(ty,2) - 3*pow(tz,2)) + 6*nz*tz*pow(ny,5)*(-3*pow(ty,2) + pow(tz,2)) - 20*tz*pow(ny,3)*pow(nz,3)*(-3*pow(ty,2) + pow(tz,2)) + 
    6*ny*tz*pow(nz,5)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,6)*(pow(ty,3) - 3*ty*pow(tz,2)),
   6*nz*ty*pow(ny,5)*(pow(ty,2) - 3*pow(tz,2)) - 20*ty*pow(ny,3)*pow(nz,3)*(pow(ty,2) - 3*pow(tz,2)) + 6*ny*ty*pow(nz,5)*(pow(ty,2) - 3*pow(tz,2)) + 
    15*tz*pow(ny,4)*pow(nz,2)*(-3*pow(ty,2) + pow(tz,2)) - 15*tz*pow(ny,2)*pow(nz,4)*(-3*pow(ty,2) + pow(tz,2)) + 
    tz*pow(nz,6)*(-3*pow(ty,2) + pow(tz,2)) + pow(ny,6)*(3*tz*pow(ty,2) - pow(tz,3)),
   pow(nx,6)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - 6*nz*pow(nx,5)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    20*pow(nx,3)*pow(nz,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*nx*pow(nz,5)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,4)*pow(nz,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) - 15*pow(nx,2)*pow(nz,4)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    pow(nz,6)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))),pow(nx,6)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
    6*pow(nx,5)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    20*pow(nx,3)*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) - 
    6*nx*pow(nz,4)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(10*rz*tx*tz - 5*rx*pow(tx,2) + 5*rx*pow(tz,2))) + 
    pow(nz,5)*(nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) - 6*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
    15*nz*pow(nx,4)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 2*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
    15*pow(nx,2)*pow(nz,3)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 4*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))),
   pow(nx,6)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 15*pow(nx,4)*pow(nz,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
    15*pow(nx,2)*pow(nz,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*nz*pow(nx,5)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    20*pow(nx,3)*pow(nz,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) - 6*nx*pow(nz,5)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))) + 
    pow(nz,6)*(-2*rx*tx*tz + rz*(-pow(tx,2) + pow(tz,2))),15*pow(nx,4)*
     (-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) - pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + pow(nx,6)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
    15*pow(nx,2)*pow(nz,2)*(8*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*ry*tx*ty - 6*rz*tx*tz + 2*rx*pow(tx,2) + rx*pow(ty,2) - 3*rx*pow(tz,2)) + 
       pow(ny,2)*(12*rz*tx*tz - 6*rx*pow(tx,2) + 6*rx*pow(tz,2))) + 
    6*pow(nx,5)*(-(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       ny*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    6*nx*pow(nz,3)*(-5*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
       pow(nz,2)*(6*rx*tx*tz + 2*ry*ty*tz + 3*rz*pow(tx,2) + rz*pow(ty,2) - 4*rz*pow(tz,2)) - 10*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))
       ) + 20*nz*pow(nx,3)*(3*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
      - pow(nz,4)*(12*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*tx*(ry*ty - 4*rz*tz) + rx*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2))) + 
       15*pow(ny,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))),
   2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,6) - 15*nz*pow(nx,4)*
       (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*rz*tx*tz - rx*pow(tx,2) + rx*pow(tz,2))) + 
      15*pow(nx,2)*pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(4*rz*tx*tz - 2*rx*pow(tx,2) + 2*rx*pow(tz,2))) - 
      pow(nz,5)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) + 
      3*pow(nx,5)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) - 
      10*pow(nx,3)*pow(nz,2)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 3*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
      3*nx*pow(nz,4)*(nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 5*ny*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))),
   20*pow(nx,3)*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*ny*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    pow(nx,6)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 6*nx*pow(nz,2)*
     (-20*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + 6*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       5*ny*pow(nz,2)*(2*ry*tx*ty - 6*rz*tx*tz + 2*rx*pow(tx,2) + rx*pow(ty,2) - 3*rx*pow(tz,2)) + 
       10*pow(ny,3)*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2))) + 
    6*pow(nx,5)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + ny*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    15*pow(nx,4)*(-(pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
       2*ny*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)))
      - pow(nz,3)*(pow(nz,3)*(6*rx*tx*ty - 8*rz*ty*tz + ry*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2))) + 
       6*ny*pow(nz,2)*(2*(3*rx*tx + ry*ty)*tz + rz*(3*pow(tx,2) + pow(ty,2) - 4*pow(tz,2))) - 
       20*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 15*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    15*nz*pow(nx,2)*(pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       4*ny*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 6*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))),
   pow(nx,6)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 6*pow(nx,5)*
     (2*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    20*nz*pow(nx,3)*(6*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,2)*(6*rz*tx*tz - 3*rx*pow(tx,2) + 3*rx*pow(tz,2))) + 
    15*pow(nx,4)*(-(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) + 
       pow(ny,2)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 2*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    pow(nz,4)*(-6*ny*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) + 
       pow(nz,2)*(6*rx*tx*tz + 2*ry*ty*tz + 3*rz*pow(tx,2) + rz*pow(ty,2) - 4*rz*pow(tz,2)) - 15*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))
       ) + 15*pow(nx,2)*pow(nz,2)*(4*ny*nz*(-2*rx*tx*ty + 2*rz*ty*tz - ry*pow(tx,2) + ry*pow(tz,2)) + 
       pow(nz,2)*(4*rx*tx*tz + 2*ry*ty*tz + 2*rz*pow(tx,2) + rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(ny,2)*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2)))
      + 6*nx*pow(nz,3)*(10*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*tx*(ry*ty - 3*rz*tz) + rx*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 10*pow(ny,2)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))),
   15*pow(nx,2)*(-8*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) + 16*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       pow(nz,4)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       6*pow(ny,2)*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,4)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    15*pow(nx,4)*(-4*ny*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(ny,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
       pow(nz,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    20*pow(nx,3)*(pow(nz,3)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       3*ny*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       3*nz*pow(ny,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    6*nx*nz*(-5*ny*pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       10*pow(ny,2)*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       3*pow(nz,4)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 5*pow(ny,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
       10*nz*pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    pow(nz,2)*(40*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) - 36*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       15*pow(ny,2)*pow(nz,2)*(2*tx*(ry*ty - 3*rz*tz) + rx*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       3*pow(nz,4)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 15*pow(ny,4)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))) - 
    6*pow(nx,5)*(nz*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + ny*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   2*(30*nz*pow(nx,2)*(-3*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
         ny*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
      15*pow(nx,4)*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
         ny*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
      3*nx*pow(nz,2)*(pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
         5*ny*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
         10*pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 10*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
      10*pow(nx,3)*(-(pow(nz,3)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)))) - 
         3*ny*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(ny,3)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 
         3*nz*pow(ny,2)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
      3*pow(nx,5)*(nz*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2)) + ny*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
      pow(nz,3)*(15*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) - 3*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
         3*ny*pow(nz,2)*(2*tx*(ry*ty - 3*rz*tz) + rx*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
         10*pow(ny,3)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2))))),
   6*nx*(-10*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,4) + 40*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2)*pow(nz,3) - 
       6*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,5) + 5*ny*pow(nz,4)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       10*pow(ny,3)*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,5)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2)))) + 
    20*pow(nx,3)*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,2) + 2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       pow(ny,3)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 3*ny*pow(nz,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)))\
     + 15*pow(nx,2)*(pow(nz,4)*(2*ty*(rx*tx - 3*rz*tz) + ry*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
       4*ny*pow(nz,3)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       6*pow(ny,2)*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       4*nz*pow(ny,3)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,4)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    nz*(15*pow(ny,2)*pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       20*pow(ny,3)*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       3*pow(nz,5)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       18*ny*pow(nz,4)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 6*pow(ny,5)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
       15*nz*pow(ny,4)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    15*pow(nx,4)*(pow(ny,2)*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2)) + pow(nz,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2)) - 
       2*ny*nz*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),20*pow(nx,3)*
     (2*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) - 6*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
       3*nz*pow(ny,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + pow(nz,3)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)))\
     + pow(nz,2)*(6*ny*pow(nz,3)*(4*rx*tx*ty - 6*rz*ty*tz + ry*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) + 
       15*pow(ny,2)*pow(nz,2)*(2*(2*rx*tx + ry*ty)*tz + rz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2))) - 
       3*pow(nz,4)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 15*pow(ny,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) - 
       20*nz*pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    15*pow(nx,2)*(pow(nz,4)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       4*ny*pow(nz,3)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       6*pow(ny,2)*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(ny,4)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 4*nz*pow(ny,3)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    15*pow(nx,4)*(2*ny*nz*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2)) + pow(ny,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)) + 
       pow(nz,2)*(-2*ry*ty*tz - rz*pow(ty,2) + rz*pow(tz,2))) + 
    6*nx*nz*(-20*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,3) + 20*ny*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
       pow(nz,4)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       10*pow(ny,2)*pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
       5*pow(ny,4)*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))),
   20*pow(ny,3)*(-6*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + 4*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       3*nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
    pow(ny,6)*(-2*rz*tx*tz + rx*(pow(tx,2) - pow(tz,2))) - pow(nz,3)*
     (pow(nz,3)*(6*ry*tx*ty - 8*rz*tx*tz + rx*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2))) + 
       6*nx*pow(nz,2)*(2*(rx*tx + 3*ry*ty)*tz + rz*(pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2))) - 
       20*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 15*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
    15*nz*pow(ny,2)*(pow(nz,3)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
       4*nx*pow(nz,2)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
       4*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 6*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) - 
    15*pow(ny,4)*(pow(nz,2)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))
       ) + 6*pow(ny,5)*(-2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) - 
    6*ny*pow(nz,2)*(-20*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + 6*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
       5*nx*pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 
       10*pow(nx,3)*(-2*rz*ty*tz + ry*pow(ty,2) - ry*pow(tz,2))),
   2*(30*nz*pow(ny,2)*(-3*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) + (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) - 
         nx*pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + pow(nx,3)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2)))) + 
      3*ny*pow(nz,2)*(pow(nz,3)*(4*ry*tx*ty - 6*rz*tx*tz + rx*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
         5*nx*pow(nz,2)*(2*(rx*tx + 2*ry*ty)*tz + rz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
         10*pow(nx,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 10*nz*pow(nx,2)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2))) + 
      15*pow(ny,4)*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - (rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,2) + 
         nx*nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
      3*pow(ny,5)*(nz*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + nx*(2*rx*tx*tz + rz*pow(tx,2) - rz*pow(tz,2))) + 
      pow(nz,3)*(15*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nx,2) - 3*(rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(nz,3) + 
         3*nx*pow(nz,2)*(2*ty*(rx*tx - 3*rz*tz) + ry*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) + 
         10*pow(nx,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))) - 
      10*pow(ny,3)*(pow(nz,3)*(2*tx*(ry*ty - 2*rz*tz) + rx*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
         3*nx*pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
         3*nz*pow(nx,2)*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,3)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))))),
   pow(ny,6)*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2)) - 
    6*pow(ny,5)*(nz*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       nx*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2))) + 
    15*pow(ny,2)*pow(nz,2)*(8*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 6*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2)))
      - pow(nz,4)*(12*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*rx*tx*ty - 8*rz*ty*tz + ry*pow(tx,2) + 3*ry*pow(ty,2) - 4*ry*pow(tz,2)) + 15*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))
       ) - 6*ny*pow(nz,3)*(-5*nx*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 6*ry*ty*tz + rz*pow(tx,2) + 3*rz*pow(ty,2) - 4*rz*pow(tz,2)) - 10*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))
       ) + 20*nz*pow(ny,3)*(3*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 3*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 15*pow(ny,4)*(4*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,6)*(2*rx*tx*tz + rz*(pow(tx,2) - pow(tz,2))) + 6*pow(ny,5)*
     (2*nx*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nz*(2*rx*tx*ty - 2*rz*ty*tz + ry*pow(tx,2) - ry*pow(tz,2))) + 
    6*ny*pow(nz,3)*(10*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 
       pow(nz,2)*(2*rx*tx*ty - 6*rz*ty*tz + ry*pow(tx,2) + 2*ry*pow(ty,2) - 3*ry*pow(tz,2)) + 10*pow(nx,2)*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))
       ) - pow(nz,4)*(-6*nx*nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 6*ry*ty*tz + rz*pow(tx,2) + 3*rz*pow(ty,2) - 4*rz*pow(tz,2)) - 15*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))
       ) + 15*pow(ny,2)*pow(nz,2)*(4*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + 
       pow(nz,2)*(2*rx*tx*tz + 4*ry*ty*tz + rz*pow(tx,2) + 2*rz*pow(ty,2) - 3*rz*pow(tz,2)) - 6*pow(nx,2)*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))
      - 20*nz*pow(ny,3)*(6*nx*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + pow(nz,2)*(2*ty*(rx*tx - 2*rz*tz) + ry*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       3*pow(nx,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2)))) - 
    15*pow(ny,4)*(pow(nz,2)*(2*(rx*tx + ry*ty)*tz + rz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 
       2*nx*nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) + pow(nx,2)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2)))),
   pow(ny,6)*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) - 
    6*pow(ny,5)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    20*pow(ny,3)*pow(nz,2)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 3*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
    6*ny*pow(nz,4)*(2*nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 5*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
    pow(nz,5)*(nz*(-2*ry*tx*ty + 2*rz*tx*tz - rx*pow(ty,2) + rx*pow(tz,2)) - 6*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) - 
    15*nz*pow(ny,4)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 2*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
    15*pow(ny,2)*pow(nz,3)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 4*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))),
   2*((rz*tx*ty + ry*tx*tz + rx*ty*tz)*pow(ny,6) - 15*nz*pow(ny,4)*
       (nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
      15*pow(ny,2)*pow(nz,3)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 2*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) - 
      pow(nz,5)*(nz*(rz*tx*ty + ry*tx*tz + rx*ty*tz) + 3*nx*(2*rz*ty*tz - ry*pow(ty,2) + ry*pow(tz,2))) + 
      3*pow(ny,5)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) - 
      10*pow(ny,3)*pow(nz,2)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 3*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2))) + 
      3*ny*pow(nz,4)*(nz*(2*ry*tx*ty - 2*rz*tx*tz + rx*pow(ty,2) - rx*pow(tz,2)) + 5*nx*(2*ry*ty*tz + rz*pow(ty,2) - rz*pow(tz,2)))),
   pow(ny,6)*(-2*rz*ty*tz + ry*(pow(ty,2) - pow(tz,2))) - 6*nz*pow(ny,5)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    20*pow(ny,3)*pow(nz,3)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 6*ny*pow(nz,5)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 
    15*pow(ny,4)*pow(nz,2)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) - 15*pow(ny,2)*pow(nz,4)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + 
    pow(nz,6)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))),pow(ny,6)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    15*pow(ny,4)*pow(nz,2)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) + 15*pow(ny,2)*pow(nz,4)*(2*ry*ty*tz + rz*(pow(ty,2) - pow(tz,2))) - 
    6*nz*pow(ny,5)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + 20*pow(ny,3)*pow(nz,3)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) - 
    6*ny*pow(nz,5)*(2*rz*ty*tz + ry*(-pow(ty,2) + pow(tz,2))) + pow(nz,6)*(-2*ry*ty*tz + rz*(-pow(ty,2) + pow(tz,2))),
   -20*nz*tx*tz*pow(nx,4)*(pow(tx,2) - pow(tz,2)) + 40*tx*tz*pow(nx,2)*pow(nz,3)*(pow(tx,2) - pow(tz,2)) + 
    4*tx*tz*pow(nz,5)*(-pow(tx,2) + pow(tz,2)) + pow(nx,5)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    10*pow(nx,3)*pow(nz,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nx*pow(nz,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   -40*nz*tx*pow(nx,3)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*ny*tz*(pow(tx,2) - pow(tz,2))) + 
    20*nx*tx*pow(nz,3)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 4*ny*tz*(pow(tx,2) - pow(tz,2))) + 4*tx*ty*pow(nx,5)*(pow(tx,2) - 3*pow(tz,2)) + 
    5*pow(nx,4)*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    10*pow(nx,2)*pow(nz,2)*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + 5*ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   4*tx*tz*pow(nx,5)*(pow(tx,2) - pow(tz,2)) - 40*tx*tz*pow(nx,3)*pow(nz,2)*(pow(tx,2) - pow(tz,2)) + 
    20*nx*tx*tz*pow(nz,4)*(pow(tx,2) - pow(tz,2)) + 5*nz*pow(nx,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    10*pow(nx,2)*pow(nz,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + pow(nz,5)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)),
   2*(-2*tx*pow(nz,3)*(tz*pow(nz,2)*(3*pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) - 5*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 
         10*tz*pow(ny,2)*(-pow(tx,2) + pow(tz,2))) + 10*tx*pow(nx,4)*
       (ny*ty*(pow(tx,2) - 3*pow(tz,2)) + nz*tz*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
      20*nz*tx*pow(nx,2)*(3*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
      pow(nx,5)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
      5*nx*pow(nz,2)*(8*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
      5*pow(nx,3)*(8*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(5*tx*pow(nx,4)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + ny*tz*(pow(tx,2) - pow(tz,2))) - 
      10*tx*pow(nx,2)*pow(nz,2)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 3*ny*tz*(pow(tx,2) - pow(tz,2))) + 
      tx*pow(nz,4)*(nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 5*ny*tz*(pow(tx,2) - pow(tz,2))) - ty*tz*pow(nx,5)*(-3*pow(tx,2) + pow(tz,2)) - 
      5*nx*pow(nz,3)*(nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      5*nz*pow(nx,3)*(2*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   2*(20*tx*pow(nx,3)*(-(ty*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) - 
         2*ny*nz*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 2*tx*ty*pow(nx,5)*(pow(ty,2) - 3*pow(tz,2)) + 
      10*nx*nz*tx*(ty*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 6*nz*ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) + 
         4*ny*tz*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + pow(ny,3)*(-4*tz*pow(tx,2) + 4*pow(tz,3))) + 
      5*pow(nx,4)*(-2*nz*ty*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         ny*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,2)*(-2*ty*tz*pow(nz,3)*(9*pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) - 20*nz*ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2)) - 
         5*pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         5*ny*pow(nz,2)*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
      5*pow(nx,2)*(4*ty*tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 12*nz*ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2)) + 
         pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         3*ny*pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(20*tx*pow(nx,3)*(2*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
      10*nx*tx*pow(nz,2)*(4*ny*nz*ty*(pow(tx,2) - 3*pow(tz,2)) + 6*tz*pow(ny,2)*(pow(tx,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + tx*pow(nx,5)*(6*tz*pow(ty,2) - 2*pow(tz,3)) + 
      5*pow(nx,4)*(-2*ny*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,3)*(-10*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) - 5*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
      5*nz*pow(nx,2)*(12*ny*nz*ty*tz*(-3*pow(tx,2) + pow(tz,2)) + 3*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   40*tx*pow(nx,2)*(-3*ny*ty*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + ty*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2)) + 
       tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) - 3*nz*tz*pow(ny,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) - 
    4*nz*tx*(-5*ny*ty*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + 10*nz*ty*pow(ny,3)*(pow(tx,2) - 3*pow(tz,2)) + 
       3*tz*pow(nz,4)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 5*tz*pow(ny,4)*(pow(tx,2) - pow(tz,2)) + 
       10*tz*pow(ny,2)*pow(nz,2)*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
    20*tx*pow(nx,4)*(nz*tz*(-3*pow(ty,2) + pow(tz,2)) + ny*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(nx,5)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    10*pow(nx,3)*(-8*ny*nz*ty*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       2*pow(ny,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    5*nx*(16*ny*ty*tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 16*nz*ty*tz*pow(ny,3)*(-3*pow(tx,2) + pow(tz,2)) + 
       pow(ny,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
       6*pow(ny,2)*pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       pow(nz,4)*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))),
   4*(10*tx*pow(nx,2)*(-(ty*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*nz*ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) - 
         3*ny*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + tz*pow(ny,3)*(pow(tx,2) - pow(tz,2))) + 
      ty*tz*pow(nx,5)*(pow(ty,2) - pow(tz,2)) + tx*pow(nz,2)*
       (ty*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) - 10*nz*ty*pow(ny,2)*(pow(tx,2) - 3*pow(tz,2)) + 
         5*ny*tz*pow(nz,2)*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + 10*tz*pow(ny,3)*(-pow(tx,2) + pow(tz,2))) + 
      5*tx*pow(nx,4)*(3*ny*tz*pow(ty,2) + nz*pow(ty,3) - 3*nz*ty*pow(tz,2) - ny*pow(tz,3)) + 
      10*pow(nx,3)*(-(ty*tz*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2)) + 
         ny*nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      5*nx*nz*(ty*tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 6*nz*ty*tz*pow(ny,2)*(-3*pow(tx,2) + pow(tz,2)) + 
         pow(ny,3)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         ny*pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   40*ty*pow(ny,2)*(-3*nx*tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
       tz*pow(nz,3)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 3*nz*tz*pow(nx,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
    4*nz*ty*(-5*nx*tx*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 10*nz*tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
       3*tz*pow(nz,4)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 5*tz*pow(nx,4)*(pow(ty,2) - pow(tz,2)) + 
       10*tz*pow(nx,2)*pow(nz,2)*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
    20*ty*pow(ny,4)*(nz*tz*(-3*pow(tx,2) + pow(tz,2)) + nx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(ny,5)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    10*pow(ny,3)*(8*nx*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    5*ny*(16*nx*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 16*nz*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 
       pow(nx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       6*pow(nx,2)*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       pow(nz,4)*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))),
   20*ny*ty*(-4*nx*tx*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 4*nz*tx*pow(nx,3)*(pow(ty,2) - 3*pow(tz,2)) + 
       tz*pow(nz,4)*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*tz*pow(nx,2)*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       tz*pow(nx,4)*(pow(ty,2) - pow(tz,2))) + 40*ty*pow(ny,3)*
     (2*nx*nz*tx*(pow(tx,2) - 3*pow(tz,2)) + tz*pow(nz,2)*(-3*pow(tx,2) + pow(tz,2)) + pow(nx,2)*(3*tz*pow(tx,2) - pow(tz,3))) + 
    5*pow(ny,4)*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    10*pow(ny,2)*(12*nx*tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 4*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) - 
       6*nz*pow(nx,2)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,3)*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    nz*(20*nx*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 40*nz*tx*tz*pow(nx,3)*(-3*pow(ty,2) + pow(tz,2)) + 
       5*pow(nx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       10*pow(nx,2)*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       pow(nz,4)*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))),
   2*(-20*ty*pow(ny,3)*(tx*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         2*nx*nz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + 2*tx*ty*pow(ny,5)*(pow(tx,2) - 3*pow(tz,2)) + 
      10*ny*nz*ty*(tx*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 6*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         4*nx*tz*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 4*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 
      5*pow(ny,4)*(-2*nz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         nx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,2)*(-2*tx*tz*pow(nz,3)*(pow(tx,2) + 9*pow(ty,2) - 4*pow(tz,2)) - 20*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) - 
         5*pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         5*nx*pow(nz,2)*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
      5*pow(ny,2)*(4*tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 12*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         3*nx*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(tx*tz*pow(ny,5)*(pow(tx,2) - pow(tz,2)) - 10*ty*pow(ny,2)*
       (tx*pow(nz,3)*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 3*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         3*nx*tz*pow(nz,2)*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 
      ty*pow(nz,2)*(tx*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 10*nz*tx*pow(nx,2)*(pow(ty,2) - 3*pow(tz,2)) + 
         5*nx*tz*pow(nz,2)*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 10*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2))) + 
      5*ty*pow(ny,4)*(3*nx*tz*pow(tx,2) + nz*pow(tx,3) - 3*nz*tx*pow(tz,2) - nx*pow(tz,3)) + 
      10*pow(ny,3)*(-(tx*tz*pow(nz,2)*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) - tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         nx*nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      5*ny*nz*(tx*tz*pow(nz,3)*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 6*nz*tx*tz*pow(nx,2)*(-3*pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         nx*pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(-2*ty*pow(nz,3)*(tz*pow(nz,2)*(3*pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) - 5*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 
         10*tz*pow(nx,2)*(-pow(ty,2) + pow(tz,2))) + 10*ty*pow(ny,4)*
       (nx*tx*(pow(ty,2) - 3*pow(tz,2)) + nz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) - 
      20*nz*ty*pow(ny,2)*(3*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
      pow(ny,5)*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      5*ny*pow(nz,2)*(-8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) - 
      5*pow(ny,3)*(-8*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(20*ty*pow(ny,3)*(2*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) - 
      10*ny*ty*pow(nz,2)*(4*nx*nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 6*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         tz*pow(nz,2)*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + pow(ny,5)*(6*ty*tz*pow(tx,2) - 2*ty*pow(tz,3)) + 
      5*pow(ny,4)*(-2*nx*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,3)*(-10*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 5*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) - 
      5*nz*pow(ny,2)*(-12*nx*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   -40*nz*ty*pow(ny,3)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 2*nx*tz*(pow(ty,2) - pow(tz,2))) + 
    20*ny*ty*pow(nz,3)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 4*nx*tz*(pow(ty,2) - pow(tz,2))) + 4*tx*ty*pow(ny,5)*(pow(ty,2) - 3*pow(tz,2)) + 
    5*pow(ny,4)*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    10*pow(ny,2)*pow(nz,2)*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 5*nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*(5*ty*pow(ny,4)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + nx*tz*(pow(ty,2) - pow(tz,2))) - 
      10*ty*pow(ny,2)*pow(nz,2)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 3*nx*tz*(pow(ty,2) - pow(tz,2))) + 
      ty*pow(nz,4)*(nz*tx*(pow(ty,2) - 3*pow(tz,2)) + 5*nx*tz*(pow(ty,2) - pow(tz,2))) - tx*tz*pow(ny,5)*(-3*pow(ty,2) + pow(tz,2)) - 
      5*ny*pow(nz,3)*(nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      5*nz*pow(ny,3)*(2*nz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   -20*nz*ty*tz*pow(ny,4)*(pow(ty,2) - pow(tz,2)) + 40*ty*tz*pow(ny,2)*pow(nz,3)*(pow(ty,2) - pow(tz,2)) + 
    4*ty*tz*pow(nz,5)*(-pow(ty,2) + pow(tz,2)) + pow(ny,5)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    10*pow(ny,3)*pow(nz,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*ny*pow(nz,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   4*ty*tz*pow(ny,5)*(pow(ty,2) - pow(tz,2)) - 40*ty*tz*pow(ny,3)*pow(nz,2)*(pow(ty,2) - pow(tz,2)) + 
    20*ny*ty*tz*pow(nz,4)*(pow(ty,2) - pow(tz,2)) + 5*nz*pow(ny,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    10*pow(ny,2)*pow(nz,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + pow(nz,5)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)),
   pow(nx,5)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    10*pow(nx,3)*pow(nz,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    5*nx*pow(nz,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    5*nz*pow(nx,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    10*pow(nx,2)*pow(nz,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    pow(nz,5)*(-3*rx*tz*pow(tx,2) - rz*pow(tx,3) + 3*rz*tx*pow(tz,2) + rx*pow(tz,3)),
   pow(nx,5)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    10*nz*pow(nx,3)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    5*nx*pow(nz,3)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       4*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    5*pow(nx,4)*(ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) - 
    10*pow(nx,2)*pow(nz,2)*(3*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) + 
    pow(nz,4)*(5*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   5*nz*pow(nx,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
    10*pow(nx,2)*pow(nz,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nz,5)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    pow(nx,5)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
    10*pow(nx,3)*pow(nz,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    5*nx*pow(nz,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)),
   pow(nx,5)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
    10*pow(nx,3)*(pow(nz,2)*(-(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(2*pow(tx,2) + pow(ty,2) - pow(tz,2)) + 
          3*ry*ty*(-pow(tx,2) + pow(tz,2))) - 2*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    5*pow(nx,4)*(-(nz*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2)))) + 
       ny*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) - 
    5*nx*pow(nz,2)*(-4*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       6*pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(rz*tz*(9*pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 
          rx*(-2*pow(tx,3) - 3*tx*pow(ty,2) + 9*tx*pow(tz,2)))) - 
    pow(nz,3)*(pow(nz,2)*(3*rz*tx*(pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) + tz*(6*ry*tx*ty + 9*rx*pow(tx,2) + 3*rx*pow(ty,2) - 4*rx*pow(tz,2))) - 
       5*ny*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       10*pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    10*nz*pow(nx,2)*(pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       3*ny*nz*(6*rz*tx*ty*tz - 3*rx*ty*pow(tx,2) - ry*pow(tx,3) + 3*ry*tx*pow(tz,2) + 3*rx*ty*pow(tz,2)) + 
       pow(ny,2)*(-9*rx*tz*pow(tx,2) - 3*rz*pow(tx,3) + 9*rz*tx*pow(tz,2) + 3*rx*pow(tz,3))),
   pow(nx,5)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    5*pow(nx,4)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) - 
    10*pow(nx,2)*pow(nz,2)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       3*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    pow(nz,4)*(nz*(-6*rz*tx*ty*tz + 3*rx*ty*pow(tx,2) + ry*pow(tx,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       5*ny*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    10*nz*pow(nx,3)*(2*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))) - 
    5*nx*pow(nz,3)*(4*ny*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(nx,5)*(ty*(-6*rz*tx*tz + rx*(pow(ty,2) - 3*pow(tz,2))) + 3*ry*tx*(pow(ty,2) - pow(tz,2))) + 
    5*pow(nx,4)*(-(nz*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
       ny*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)))) - 
    pow(nz,2)*(pow(nz,3)*(tz*(18*rx*tx*ty + ry*(9*pow(tx,2) + 3*pow(ty,2) - 4*pow(tz,2))) + rz*ty*(9*pow(tx,2) + pow(ty,2) - 12*pow(tz,2))) - 
       5*ny*pow(nz,2)*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) + 3*ry*ty*(pow(tx,2) - pow(tz,2)) + 
          rz*tz*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2))) - 
       10*nz*pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       10*pow(ny,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    10*pow(nx,3)*(-(pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)))) - 
       2*ny*nz*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2))) + 
       pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    10*pow(nx,2)*(pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       3*nz*pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       3*ny*pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
          3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2))))) + 
    5*nx*nz*(pow(nz,3)*(ty*(-18*rz*tx*tz + rx*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       4*ny*pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       6*nz*pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       4*pow(ny,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))),
   pow(nx,5)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    5*pow(nx,4)*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       ny*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2)))) + 
    pow(nz,3)*(pow(nz,2)*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2)) + 3*ry*ty*(pow(tx,2) - pow(tz,2)) + 
          rz*tz*(-9*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2))) + 
       5*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) - 
       10*pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2)))) - 
    10*nz*pow(nx,2)*(3*ny*nz*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) - 
       3*pow(ny,2)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + 5*nx*pow(nz,2)*(pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       4*ny*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       6*pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    10*pow(nx,3)*(-(pow(nz,2)*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
            tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2)))) + 
       2*ny*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(ny,2)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))),
   pow(nx,5)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    10*pow(nx,3)*(-2*ny*nz*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       pow(nz,2)*(-(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 3*rz*tz*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2)) + 
          3*rx*tx*(-pow(ty,2) + pow(tz,2))) + pow(ny,2)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
          rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2)))) + 
    10*pow(nx,2)*(pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
       3*ny*pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) - 
       3*nz*pow(ny,2)*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2))) + 
       pow(ny,3)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    5*nx*(4*ny*pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,4)*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + 
          3*rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
       4*nz*pow(ny,3)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       pow(ny,4)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       6*pow(ny,2)*pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
          3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2))))) + 
    nz*(5*ny*pow(nz,3)*(ty*(-18*rz*tx*tz + rx*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       10*pow(ny,2)*pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       3*pow(nz,4)*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2))) - 
       10*nz*pow(ny,3)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       5*pow(ny,4)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    5*pow(nx,4)*(ny*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       nz*(-6*ry*tx*ty*tz - 3*rx*tz*pow(ty,2) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)) + rx*pow(tz,3))),
   5*pow(nx,4)*(nz*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       ny*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)))) + 
    10*pow(nx,3)*(-(pow(nz,2)*(tz*(6*rx*tx*ty + ry*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)))) + 
       2*ny*nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2)))) + 
    5*nx*nz*(pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       6*nz*pow(ny,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       4*pow(ny,3)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       4*ny*pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 
          3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2))))) + 
    pow(nx,5)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    pow(nz,2)*(pow(nz,3)*(ty*(-18*rz*tx*tz + rx*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) + 
       5*ny*pow(nz,2)*(3*tz*(2*ry*tx*ty + rx*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*tx*(2*pow(tx,2) + 3*pow(ty,2) - 9*pow(tz,2))) - 
       10*nz*pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) - 
       10*pow(ny,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    10*pow(nx,2)*(-(pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)))) - 
       3*ny*pow(nz,2)*(rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(tx,2) + 3*rx*pow(ty,2) - 2*rx*pow(tz,2))) + 
       3*nz*pow(ny,2)*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(ny,3)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))),
   pow(ny,5)*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
    10*pow(ny,2)*(pow(nz,3)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       3*nx*pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,3)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
       3*nz*pow(nx,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) + 
    5*ny*(4*nx*pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,4)*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + 
          3*rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
       4*nz*pow(nx,3)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       pow(nx,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       6*pow(nx,2)*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) - 
    10*pow(ny,3)*(2*nx*nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + nz*(10*pow(nx,2)*pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       5*nx*pow(nz,3)*(ry*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + ty*(-18*rz*tx*tz + 3*rx*pow(tx,2) + 2*rx*pow(ty,2) - 9*rx*pow(tz,2))) - 
       10*nz*pow(nx,3)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
       3*pow(nz,4)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) - 
       5*pow(nx,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    5*pow(ny,4)*(nx*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nz*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   nz*(5*nx*pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,4)*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + 
          3*rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
       10*nz*pow(nx,3)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       5*pow(nx,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       10*pow(nx,2)*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) - 
    10*pow(ny,2)*(3*nx*pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       3*nz*pow(nx,2)*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nx,3)*(-6*ry*tx*ty*tz + rx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*rz*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,3)*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2)) + 3*(ry*ty*(pow(tx,2) - pow(tz,2)) + rz*tz*(-2*pow(tx,2) - pow(ty,2) + pow(tz,2)))))\
     + 5*pow(ny,4)*(nz*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       nx*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3))) + 
    5*ny*(pow(nz,4)*(3*tz*(4*rx*tx*ty + ry*(2*pow(tx,2) + pow(ty,2) - pow(tz,2))) + rz*ty*(6*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       4*nx*pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       4*nz*pow(nx,3)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
       6*pow(nx,2)*pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    10*pow(ny,3)*(pow(nx,2)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       2*nx*nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2))) + 
       pow(nz,2)*(-6*rx*tx*ty*tz - 3*ry*tz*pow(tx,2) + 3*rz*ty*(-pow(tx,2) + pow(tz,2)) + ry*pow(tz,3))),
   pow(ny,5)*(ry*(pow(tx,3) - 3*tx*pow(tz,2)) - 3*ty*(2*rz*tx*tz + rx*(-pow(tx,2) + pow(tz,2)))) - 
    5*pow(ny,4)*(nz*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       nx*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2)))) - 
    10*pow(ny,3)*(pow(nz,2)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       pow(nx,2)*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       2*nx*nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) - 
    pow(nz,2)*(pow(nz,3)*(tz*(18*ry*tx*ty + rx*(3*pow(tx,2) + 9*pow(ty,2) - 4*pow(tz,2))) + rz*tx*(pow(tx,2) + 9*pow(ty,2) - 12*pow(tz,2))) - 
       5*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
          rz*tz*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2))) - 
       10*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       10*pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2)))) + 
    10*pow(ny,2)*(pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
       3*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       3*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) + 
    5*ny*nz*(4*nx*pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,3)*(ry*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + ty*(-18*rz*tx*tz + 3*rx*pow(tx,2) + 2*rx*pow(ty,2) - 9*rx*pow(tz,2))) - 
       6*nz*pow(nx,2)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       4*pow(nx,3)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   -10*pow(ny,3)*(pow(nz,2)*(tz*(6*ry*tx*ty + rx*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + rz*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       2*nx*nz*(rz*tz*(3*pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 3*ry*ty*(-pow(tx,2) + pow(tz,2)) + 3*rx*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nx,2)*(-6*ry*tx*ty*tz + rx*tz*(-3*pow(ty,2) + pow(tz,2)) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)))) + 
    5*pow(ny,4)*(nx*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
       nz*(3*ty*(-2*rz*tx*tz + rx*pow(tx,2) - rx*pow(tz,2)) + ry*(pow(tx,3) - 3*tx*pow(tz,2)))) + 
    5*ny*nz*(pow(nz,3)*(3*tz*(4*ry*tx*ty + rx*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2))) - 
       6*nz*pow(nx,2)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
       4*pow(nx,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
       4*nx*pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 
          3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2))))) + 
    pow(ny,5)*(3*rx*tz*pow(tx,2) + rz*pow(tx,3) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) + 
    pow(nz,2)*(5*nx*pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       pow(nz,3)*(ry*tx*(pow(tx,2) + 6*pow(ty,2) - 9*pow(tz,2)) + ty*(-18*rz*tx*tz + 3*rx*pow(tx,2) + 2*rx*pow(ty,2) - 9*rx*pow(tz,2))) - 
       10*nz*pow(nx,2)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) - 
       10*pow(nx,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) - 
    10*pow(ny,2)*(pow(nz,3)*(ty*(-12*rz*tx*tz + rx*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 6*pow(tz,2))) + 
       3*nz*pow(nx,2)*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*nx*pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,3)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,5)*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) - 
    5*pow(ny,4)*(nx*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       nz*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2)))) + 
    5*ny*pow(nz,2)*(pow(nz,2)*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
          rz*tz*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2))) + 
       4*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       6*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2)))) - 
    10*pow(ny,3)*(2*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     - pow(nz,3)*(-5*nx*nz*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       pow(nz,2)*(3*rz*ty*(pow(tx,2) + pow(ty,2) - 4*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 9*ry*pow(ty,2) - 4*ry*pow(tz,2))) - 
       10*pow(nx,2)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    10*nz*pow(ny,2)*(pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) + 
       3*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       3*pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,5)*(3*rz*ty*(pow(tx,2) - pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) - ry*pow(tz,2))) + 
    5*pow(ny,4)*(nz*(3*ry*ty*(pow(tx,2) - pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       nx*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2)))) + 
    pow(nz,3)*(pow(nz,2)*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 3*rx*tx*(pow(ty,2) - pow(tz,2)) + 
          rz*tz*(-3*pow(tx,2) - 9*pow(ty,2) + 4*pow(tz,2))) + 
       5*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       10*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2)))) - 
    10*nz*pow(ny,2)*(3*nx*nz*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) - 
       3*pow(nx,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       pow(nz,2)*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + 3*(rx*tx*(pow(ty,2) - pow(tz,2)) + rz*tz*(-pow(tx,2) - 2*pow(ty,2) + pow(tz,2)))))\
     - 10*pow(ny,3)*(2*nx*nz*(6*rz*tx*ty*tz - rx*pow(ty,3) + 3*rx*ty*pow(tz,2) + 3*ry*tx*(-pow(ty,2) + pow(tz,2))) + 
       pow(nz,2)*(rz*ty*(3*pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + tz*(6*rx*tx*ty + 3*ry*pow(tx,2) + 3*ry*pow(ty,2) - 2*ry*pow(tz,2))) + 
       pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))) + 
    5*ny*pow(nz,2)*(pow(nz,2)*(3*tz*(2*rx*tx*ty + ry*(pow(tx,2) + 2*pow(ty,2) - pow(tz,2))) + rz*ty*(3*pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 
       4*nx*nz*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
       6*pow(nx,2)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3))),
   pow(ny,5)*(3*ry*tx*(pow(ty,2) - pow(tz,2)) + ty*(-6*rz*tx*tz + rx*pow(ty,2) - 3*rx*pow(tz,2))) + 
    pow(nz,4)*(5*nx*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
       nz*(-6*ry*tx*ty*tz - 3*rx*tz*pow(ty,2) + 3*rz*tx*(-pow(ty,2) + pow(tz,2)) + rx*pow(tz,3))) - 
    10*nz*pow(ny,3)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       2*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    5*ny*pow(nz,3)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       4*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    10*pow(ny,2)*pow(nz,2)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       3*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))) - 
    5*pow(ny,4)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,5)*(3*rz*tx*(pow(ty,2) - pow(tz,2)) + tz*(6*ry*tx*ty + 3*rx*pow(ty,2) - rx*pow(tz,2))) + 
    5*pow(ny,4)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) - 
    10*pow(ny,2)*pow(nz,2)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       3*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    pow(nz,4)*(nz*(-6*rz*tx*ty*tz + 3*ry*tx*pow(ty,2) + rx*pow(ty,3) - 3*ry*tx*pow(tz,2) - 3*rx*ty*pow(tz,2)) + 
       5*nx*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3))) + 
    5*ny*pow(nz,3)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       4*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))) - 
    10*nz*pow(ny,3)*(nz*(6*ry*tx*ty*tz + 3*rz*tx*pow(ty,2) + 3*rx*tz*pow(ty,2) - 3*rz*tx*pow(tz,2) - rx*pow(tz,3)) - 
       2*nx*(-3*rz*tz*pow(ty,2) + ry*pow(ty,3) - 3*ry*ty*pow(tz,2) + rz*pow(tz,3))),
   pow(ny,5)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    10*pow(ny,3)*pow(nz,2)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    5*ny*pow(nz,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    5*nz*pow(ny,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    10*pow(ny,2)*pow(nz,3)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    pow(nz,5)*(-3*ry*tz*pow(ty,2) - rz*pow(ty,3) + 3*rz*ty*pow(tz,2) + ry*pow(tz,3)),
   5*nz*pow(ny,4)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) - 
    10*pow(ny,2)*pow(nz,3)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(nz,5)*(rz*tz*(-3*pow(ty,2) + pow(tz,2)) + ry*(pow(ty,3) - 3*ty*pow(tz,2))) + 
    pow(ny,5)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) - 
    10*pow(ny,3)*pow(nz,2)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)) + 
    5*ny*pow(nz,4)*(3*ry*tz*pow(ty,2) + rz*pow(ty,3) - 3*rz*ty*pow(tz,2) - ry*pow(tz,3)),
   -4*nz*tz*pow(nx,3)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 4*nx*tz*pow(nz,3)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    6*tx*pow(nx,2)*pow(nz,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + tx*pow(nz,4)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(nx,4)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)),
   5*ty*pow(nx,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
    6*nz*pow(nx,2)*(2*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(4*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    4*tx*pow(nx,3)*(20*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) - 
    4*nx*tx*pow(nz,2)*(20*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + 3*ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))),
   -6*tz*pow(nx,2)*pow(nz,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + tz*pow(nz,4)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    4*nz*tx*pow(nx,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 4*nx*tx*pow(nz,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(nx,4)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)),
   2*(5*tx*pow(nx,4)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      2*pow(nx,3)*(nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
         5*ny*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      3*tx*pow(nx,2)*(40*ny*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) - 
         pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
         pow(ny,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      tx*pow(nz,2)*(40*ny*nz*ty*tz*(pow(tx,2) - pow(tz,2)) - 3*pow(ny,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
         pow(nz,2)*(pow(tx,4) + 5*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 10*pow(tz,4))) - 
      2*nx*nz*(15*ny*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         2*tz*pow(nz,2)*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
         3*pow(ny,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)))),
   4*(5*tx*ty*tz*pow(nx,4)*(pow(tx,2) - pow(tz,2)) + pow(nx,3)*
       (ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
      nx*pow(nz,2)*(3*ny*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 5*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
      tx*pow(nz,3)*(5*nz*ty*tz*(pow(tx,2) - pow(tz,2)) - ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      3*nz*tx*pow(nx,2)*(10*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) + ny*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)))),
   2*(5*ty*pow(nx,4)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      20*tx*pow(nx,3)*(-2*nz*ty*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         ny*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      3*pow(nx,2)*(5*ty*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         2*ny*nz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
         5*ty*pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      2*nx*tx*(20*ty*tz*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 60*nz*ty*tz*pow(ny,2)*(-pow(tx,2) + pow(tz,2)) - 
         3*ny*pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
         pow(ny,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      nz*(-15*nz*ty*pow(ny,2)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         4*ny*tz*pow(nz,2)*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
         5*ty*pow(nz,3)*(pow(tx,4) + pow(tx,2)*(pow(ty,2) - 9*pow(tz,2)) - pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
         2*pow(ny,3)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)))),
   2*(20*tx*pow(nx,3)*(2*ny*ty*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      2*nx*nz*tx*(60*ny*nz*ty*tz*(-pow(tx,2) + pow(tz,2)) - 
         pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
         3*pow(ny,2)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      pow(nx,4)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)) + 
      pow(nz,2)*(-10*ny*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         tz*pow(nz,2)*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
         3*pow(ny,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5))) + 
      3*pow(nx,2)*(tz*pow(nz,2)*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
         10*ny*nz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + pow(ny,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)))),
   5*tx*pow(nx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    4*nz*tz*pow(nx,3)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
    30*tx*pow(nx,2)*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    4*nx*tz*pow(nz,3)*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
    pow(ny,4)*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
    20*ty*pow(ny,3)*(4*nz*tx*tz*(-pow(tx,2) + pow(tz,2)) + nx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(ny,2)*(tx*pow(nz,2)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) - 
       10*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*nx*nz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    20*ny*ty*(4*tx*tz*pow(nz,3)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 12*nz*tx*tz*pow(nx,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
       2*pow(nx,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       3*nx*pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    tx*pow(nz,4)*(pow(tx,4) + 10*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) + 5*(pow(ty,4) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))),
   4*(15*ty*pow(ny,2)*(2*tx*tz*pow(nx,2)*(pow(tx,2) - pow(tz,2)) + 2*tx*tz*pow(nz,2)*(-pow(tx,2) + pow(tz,2)) + 
         nx*nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
      ny*(tx*pow(nz,3)*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) - 
         2*tz*pow(nx,3)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         30*nz*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         3*nx*tz*pow(nz,2)*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      5*ty*(tx*tz*pow(nz,4)*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) - 6*tx*tz*pow(nx,2)*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         tx*tz*pow(nx,4)*(pow(ty,2) - pow(tz,2)) + 2*nz*pow(nx,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         nx*pow(nz,3)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      pow(ny,3)*(nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)))),
   5*ty*pow(ny,4)*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
    4*pow(ny,3)*(nz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       10*nx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    30*ty*pow(ny,2)*(8*nx*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) - 
       2*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    4*ny*(5*tx*pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       3*nz*tz*pow(nx,2)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
       15*nx*tx*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*pow(nz,3)*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    ty*(80*nx*tx*tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 80*nz*tx*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2)) + 
       pow(nx,4)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
       6*pow(nx,2)*pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       pow(nz,4)*(5*pow(tx,4) + pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - 9*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 30*pow(tz,4))),
   20*nz*tx*pow(nx,3)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    6*tz*pow(nx,2)*pow(nz,2)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
    20*nx*tx*pow(nz,3)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    tz*pow(nz,4)*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
    20*ty*pow(ny,3)*(4*nx*tx*tz*(pow(tx,2) - pow(tz,2)) + nz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    20*ny*ty*(12*nx*tx*tz*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 4*tx*tz*pow(nx,3)*(-pow(ty,2) + pow(tz,2)) - 
       6*nz*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,3)*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(ny,4)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5)) + pow(nx,4)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) + 
    6*pow(ny,2)*(tz*pow(nz,2)*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       20*nx*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       2*pow(nx,2)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5))),
   2*(5*tx*pow(ny,4)*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      20*ty*pow(ny,3)*(-2*nz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
         nx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      nz*(-2*tz*pow(nx,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 15*nz*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         4*nx*tz*pow(nz,2)*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
         5*tx*pow(nz,3)*(pow(ty,4) + pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) - 
      3*pow(ny,2)*(-5*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         2*nx*nz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         5*tx*pow(nz,2)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      2*ny*ty*(20*tx*tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 60*nz*tx*tz*pow(nx,2)*(-pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
         3*nx*pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))),
   4*(5*tx*ty*tz*pow(ny,4)*(pow(tx,2) - pow(tz,2)) + 2*pow(ny,3)*
       (nx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         5*nz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      30*ty*pow(ny,2)*(-(tx*tz*pow(nz,2)*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) + tx*tz*pow(nx,2)*(pow(ty,2) - pow(tz,2)) + 
         nx*nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ny*(tz*pow(nx,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 15*nz*tx*pow(nx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         3*nx*tz*pow(nz,2)*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
         5*tx*pow(nz,3)*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      nz*ty*(5*tx*tz*pow(nz,3)*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 30*nz*tx*tz*pow(nx,2)*(-pow(ty,2) + pow(tz,2)) + 
         pow(nx,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
         nx*pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))),
   2*(5*ty*pow(ny,4)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      2*pow(ny,3)*(nz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
         5*nx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ty*pow(nz,2)*(40*nx*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - 3*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 5*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 10*pow(tz,4))) - 
      3*ty*pow(ny,2)*(40*nx*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) - 
      2*ny*nz*(15*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         2*tz*pow(nz,2)*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
         3*pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)))),
   2*(20*ty*pow(ny,3)*(2*nx*tx*tz*(pow(ty,2) - pow(tz,2)) + nz*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
      2*ny*nz*ty*(60*nx*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - 3*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         pow(nz,2)*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
      pow(ny,4)*(5*pow(tx,2)*(3*tz*pow(ty,2) - pow(tz,3)) - 5*pow(ty,2)*pow(tz,3) + pow(tz,5)) + 
      pow(nz,2)*(-10*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         tz*pow(nz,2)*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
         3*pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5))) + 
      3*pow(ny,2)*(tz*pow(nz,2)*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
         10*nx*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + pow(nx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)))),
   5*tx*pow(ny,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    6*nz*pow(ny,2)*(2*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(4*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    4*ty*pow(ny,3)*(20*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) - 
    4*ny*ty*pow(nz,2)*(20*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + 3*nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))),
   4*(5*tx*ty*tz*pow(ny,4)*(pow(ty,2) - pow(tz,2)) + pow(ny,3)*
       (nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
      ny*pow(nz,2)*(3*nx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 5*nz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      ty*pow(nz,3)*(5*nz*tx*tz*(pow(ty,2) - pow(tz,2)) - nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      3*nz*ty*pow(ny,2)*(10*nz*tx*tz*(-pow(ty,2) + pow(tz,2)) + nx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))),
   -4*nz*tz*pow(ny,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 4*ny*tz*pow(nz,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
    6*ty*pow(ny,2)*pow(nz,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + ty*pow(nz,4)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(ny,4)*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)),
   -6*tz*pow(ny,2)*pow(nz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + tz*pow(nz,4)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    4*nz*ty*pow(ny,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 4*ny*ty*pow(nz,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
    pow(ny,4)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)),
   pow(nx,4)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(nx,2)*pow(nz,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    4*nz*pow(nx,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    4*nx*pow(nz,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   pow(nx,4)*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    4*pow(nx,3)*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    4*nx*pow(nz,2)*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       3*ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    6*nz*pow(nx,2)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       2*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    pow(nz,3)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
          ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       4*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))),
   4*nz*pow(nx,3)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    4*nx*pow(nz,3)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,4)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(nx,2)*pow(nz,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))),
   2*(pow(nx,4)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,2)*(8*ny*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         3*pow(ny,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(2*tx*(ry*ty*(pow(tx,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2))) + 
            rx*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
      2*nx*nz*(3*ny*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         3*pow(ny,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         2*pow(nz,2)*(2*tz*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
      3*pow(nx,2)*(-8*ny*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         pow(ny,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(4*tx*(-(ry*ty*(pow(tx,2) - 3*pow(tz,2))) + rz*tz*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) - 
            rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      2*pow(nx,3)*(ny*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   4*(pow(nx,4)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
      pow(nz,3)*(nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
      3*nz*pow(nx,2)*(2*nz*(-3*rx*ty*tz*pow(tx,2) - ry*tz*pow(tx,3) - rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + ry*tx*pow(tz,3) + rx*ty*pow(tz,3)) + 
         ny*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
      pow(nx,3)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
      nx*pow(nz,2)*(nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         3*ny*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   2*(pow(nx,4)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      4*pow(nx,3)*(-2*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
            rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
         ny*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
            rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      nz*(-3*nz*pow(ny,2)*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         2*pow(ny,3)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,3)*(2*ty*(rx*tx*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + rz*tz*(-9*pow(tx,2) - pow(ty,2) + 4*pow(tz,2))) + 
            ry*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
         4*ny*pow(nz,2)*(2*tz*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
      2*nx*(4*pow(nz,3)*(tz*(rx*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) + 
            rz*tx*ty*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
         12*nz*pow(ny,2)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         pow(ny,3)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         3*ny*pow(nz,2)*(4*tx*(ry*ty*(pow(tx,2) - 3*pow(tz,2)) + rz*tz*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
            rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      3*pow(nx,2)*(pow(ny,2)*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
            ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
         2*ny*nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   2*(pow(nx,4)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
         rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      4*pow(nx,3)*(2*ny*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))) + 
         nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
            rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      pow(nz,2)*(-2*ny*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         3*pow(ny,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(2*tz*(rx*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 3*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
      2*nx*nz*(-12*ny*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         3*pow(ny,2)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(4*tx*(-(ry*ty*(pow(tx,2) - 3*pow(tz,2))) + rz*tz*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) - 
            rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      3*pow(nx,2)*(2*ny*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + 
            ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         pow(ny,2)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) - 
         pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   pow(ny,4)*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,4)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(nx,2)*pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
       rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    4*nz*pow(nx,3)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(nz,4)*(4*tx*(ry*ty*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2)) + 3*rz*tz*(-pow(tx,2) - 3*pow(ty,2) + 2*pow(tz,2))) + 
       rx*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    4*nx*pow(nz,3)*(4*tz*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
       rz*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    4*pow(ny,3)*(-4*nz*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nx*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    4*ny*(4*pow(nz,3)*(tz*(rx*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) + 
          rz*tx*ty*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
       12*nz*pow(nx,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       2*pow(nx,3)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       3*nx*pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
          ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    6*pow(ny,2)*(2*pow(nx,2)*(-2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tx*(ry*ty*(pow(tx,2) - 3*pow(tz,2)) + rz*tz*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
          rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       2*nx*nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   4*(pow(nz,4)*(tz*(rx*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + ry*tx*(2*pow(tx,2) + 3*pow(ty,2) - 3*pow(tz,2))) + 
         rz*tx*ty*(2*pow(tx,2) + pow(ty,2) - 9*pow(tz,2))) - 
      6*pow(nx,2)*pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
         rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + pow(nx,4)*
       (tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
      2*nz*pow(nx,3)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
      nx*pow(nz,3)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
         ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      3*pow(ny,2)*(2*pow(nx,2)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         2*pow(nz,2)*(-3*rx*ty*tz*pow(tx,2) - ry*tz*pow(tx,3) - rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + ry*tx*pow(tz,3) + rx*ty*pow(tz,3)) + 
         nx*nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
      pow(ny,3)*(nz*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
      ny*(-6*nz*pow(nx,2)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
            rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         2*pow(nx,3)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
            rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,3)*(4*tx*(ry*ty*(pow(tx,2) - 3*pow(tz,2)) + rz*tz*(-2*pow(tx,2) - 3*pow(ty,2) + 3*pow(tz,2))) + 
            rx*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
         3*nx*pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
            rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   16*nx*pow(nz,3)*(tz*(rx*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + ry*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
       rz*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 16*nz*pow(nx,3)*
     (tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
    pow(ny,4)*(4*ty*(rz*tz*(-3*pow(tx,2) + pow(tz,2)) + rx*(pow(tx,3) - 3*tx*pow(tz,2))) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,4)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(nx,2)*pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
       ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(nz,4)*(4*ty*(rx*tx*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) + 3*rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
       ry*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) - 
    6*pow(ny,2)*(8*nx*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
       2*pow(nx,2)*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
          ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    4*pow(ny,3)*(2*nx*(-2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) - 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       nz*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    4*ny*(pow(nx,3)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       3*nx*pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
       3*nz*pow(nx,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
          rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       pow(nz,3)*(4*tz*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
          rz*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)))),
   pow(ny,4)*(4*rx*tx*tz*(pow(tx,2) - pow(tz,2)) + rz*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
    4*nz*pow(nx,3)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
       rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nx,4)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    4*nx*pow(nz,3)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
       rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
    6*pow(nx,2)*pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
       rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
    pow(nz,4)*(4*tz*(ry*ty*(6*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + rx*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
       rz*(pow(tx,4) + pow(ty,4) + 6*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    4*pow(ny,3)*(4*nx*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(3*rx*ty*pow(tx,2) + ry*pow(tx,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
       nz*(4*ty*(-3*rz*tz*pow(tx,2) + rx*pow(tx,3) - 3*rx*tx*pow(tz,2) + rz*pow(tz,3)) + ry*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    4*ny*(12*nx*pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
          rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) - 
       4*pow(nx,3)*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       6*nz*pow(nx,2)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
          ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,3)*(4*ty*(rx*tx*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) + rz*tz*(-6*pow(tx,2) - pow(ty,2) + 3*pow(tz,2))) + 
          ry*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    6*pow(ny,2)*(-4*nx*nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
          rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
       2*pow(nx,2)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
          rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       pow(nz,2)*(4*tz*(rx*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + ry*ty*(3*pow(tx,2) - pow(tz,2))) + 
          rz*(pow(tx,4) + 6*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 6*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))),
   2*(pow(ny,4)*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
         rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
      4*pow(ny,3)*(2*nz*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
            rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
         nx*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
            ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      nz*(-3*nz*pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         2*pow(nx,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,3)*(2*ry*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 9*pow(ty,2) - 4*pow(tz,2)) + 
            rx*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
         4*nx*pow(nz,2)*(2*tz*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
      2*ny*(4*pow(nz,3)*(tz*(rx*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + ry*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
            rz*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 
         12*nz*pow(nx,2)*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         pow(nx,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         3*nx*pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
            ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
      3*pow(ny,2)*(-(pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
              rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
         pow(nz,2)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
            rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
         2*nx*nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   4*(pow(ny,4)*(rz*tx*ty*(pow(tx,2) - 3*pow(tz,2)) + tz*(rx*ty*(3*pow(tx,2) - pow(tz,2)) + ry*(pow(tx,3) - tx*pow(tz,2)))) - 
      6*pow(ny,2)*(pow(nz,2)*(tz*(rx*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + ry*tx*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2))) + 
            rz*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2))) + 
         pow(nx,2)*(-(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + ry*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*ty*tz*(-pow(ty,2) + pow(tz,2))) + 
         nx*nz*(2*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2))) - 
            ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      2*pow(ny,3)*(nz*(2*ry*tx*ty*(pow(tx,2) - 3*pow(tz,2)) - 2*rz*tx*tz*(pow(tx,2) + 3*pow(ty,2) - 2*pow(tz,2)) + 
            rx*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
            rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      nz*(pow(nz,3)*(tz*(rx*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + ry*tx*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2))) + 
            rz*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 9*pow(tz,2))) - 
         6*nz*pow(nx,2)*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         pow(nx,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         nx*pow(nz,2)*(4*ty*(-(rx*tx*(pow(ty,2) - 3*pow(tz,2))) + rz*tz*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2))) - 
            ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      ny*(3*nz*pow(nx,2)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nx,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         pow(nz,3)*(4*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 6*pow(tz,2)) - 4*rz*tx*tz*(pow(tx,2) + 6*pow(ty,2) - 3*pow(tz,2)) + 
            rx*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) - 
         3*nx*pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   2*(pow(ny,4)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
         ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      pow(nz,2)*(8*nx*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         3*pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 3*pow(ty,2) + 4*pow(tz,2))) + 
            ry*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
      2*ny*nz*(-3*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         3*pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         2*pow(nz,2)*(2*tz*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
      3*pow(ny,2)*(8*nx*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
            ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
      2*pow(ny,3)*(-(nx*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
              rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
         nz*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   2*(pow(ny,4)*(6*ry*ty*tz*pow(tx,2) + 6*rx*tx*tz*pow(ty,2) - 2*rx*tx*pow(tz,3) - 2*ry*ty*pow(tz,3) + 
         rz*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      4*pow(ny,3)*(2*nx*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) + 
         nz*(2*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - pow(ty,2) + 2*pow(tz,2))) + 
            ry*(3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      pow(nz,2)*(-2*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         3*pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(2*tz*(ry*ty*(3*pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 3*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
      2*ny*nz*(12*nx*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         3*pow(nx,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(4*ty*(rx*tx*(pow(ty,2) - 3*pow(tz,2)) + rz*tz*(-3*pow(tx,2) - 2*pow(ty,2) + 3*pow(tz,2))) + 
            ry*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
      3*pow(ny,2)*(-2*nx*nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         pow(nx,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         pow(nz,2)*(4*tz*(ry*ty*(3*pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*tx*(3*pow(ty,2) - pow(tz,2))) + 
            rz*(pow(ty,4) + 6*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   pow(ny,4)*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    4*ny*pow(nz,2)*(4*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       3*nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    4*pow(ny,3)*(4*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
       nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    6*nz*pow(ny,2)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    pow(nz,3)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
          rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       4*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))),
   4*(pow(ny,4)*(tz*(rx*ty*(pow(ty,2) - pow(tz,2)) + ry*tx*(3*pow(ty,2) - pow(tz,2))) + rz*tx*ty*(pow(ty,2) - 3*pow(tz,2))) + 
      pow(nz,3)*(nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
      3*nz*pow(ny,2)*(2*nz*(rz*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + tz*(3*ry*tx*pow(ty,2) + rx*pow(ty,3) - ry*tx*pow(tz,2) - rx*ty*pow(tz,2))) - 
         nx*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      pow(ny,3)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
      ny*pow(nz,2)*(nz*(4*ry*tx*ty*(pow(ty,2) - 3*pow(tz,2)) + 4*rz*tx*tz*(-3*pow(ty,2) + pow(tz,2)) + 
            rx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         3*nx*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,4)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(ny,2)*pow(nz,2)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    4*nz*pow(ny,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    4*ny*pow(nz,3)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   4*nz*pow(ny,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    4*ny*pow(nz,3)*(4*rz*ty*tz*(-pow(ty,2) + pow(tz,2)) + ry*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(ny,4)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    6*pow(ny,2)*pow(nz,2)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,4)*(4*ry*ty*tz*(pow(ty,2) - pow(tz,2)) + rz*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))),
   -6*nz*tx*tz*pow(nx,2)*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
    2*tx*tz*pow(nz,3)*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(nx,3)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    3*nx*pow(nz,2)*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6)),
   3*(2*tx*ty*pow(nx,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
      2*nx*nz*tx*(2*ny*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      pow(nx,2)*(-6*nz*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         ny*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
      pow(nz,2)*(2*nz*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         ny*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6)))),
   -6*nx*tx*tz*pow(nz,2)*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(nx,3)*(6*tz*pow(tx,5) - 20*pow(tx,3)*pow(tz,3) + 6*tx*pow(tz,5)) + 
    3*nz*pow(nx,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    pow(nz,3)*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6)),
   3*(2*tx*pow(nx,2)*(nz*tz*(-3*pow(tx,4) + 30*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 9*pow(tz,4)) + 
         3*ny*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      2*nz*tx*(2*tz*pow(nz,2)*(pow(tx,4) + 5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
         3*ny*nz*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + pow(ny,2)*(-3*tz*pow(tx,4) + 10*pow(tx,2)*pow(tz,3) - 3*pow(tz,5))) + 
      nx*(-12*ny*nz*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
         pow(nz,2)*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 
            15*pow(ty,2)*pow(tz,4) - 4*pow(tz,6)) + pow(ny,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
      pow(nx,3)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   6*(ty*tz*pow(nx,3)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
      tx*pow(nz,2)*(ny*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + nz*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      tx*pow(nx,2)*(ny*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
      nx*nz*(-3*nz*ty*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
         ny*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))),
   2*ty*(10*tx*pow(nx,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       2*tz*pow(nz,3)*(15*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) - 
       3*nz*tz*pow(nx,2)*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) - 
       3*nx*tx*pow(nz,2)*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4))) + 
    18*ty*pow(ny,2)*(-(nz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + nx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    pow(ny,3)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) - 
    3*ny*(4*nx*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       pow(nz,2)*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
          4*pow(tz,6)) + 3*pow(nx,2)*(-5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) - 
          5*pow(ty,2)*pow(tz,4) + pow(tz,6))),4*tx*tz*pow(nx,3)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) - 
    6*nx*tx*tz*pow(nz,2)*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
    18*ny*ty*(-(tz*pow(nz,2)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 2*nx*nz*tx*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
       pow(nx,2)*(5*tz*pow(tx,4) - 10*pow(tx,2)*pow(tz,3) + pow(tz,5))) + 
    3*pow(ny,2)*(2*nx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) - 
    pow(nz,3)*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
       4*pow(tz,6)) + 9*nz*pow(nx,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 
       5*pow(ty,2)*pow(tz,4) - pow(tz,6)),3*(2*nz*tx*tz*pow(nx,2)*
       (-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
      2*tx*ty*pow(ny,3)*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
      2*tx*tz*pow(nz,3)*(pow(tx,4) + 5*pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
      pow(nx,3)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
      2*ny*ty*(-10*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         2*nx*nz*tz*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         tx*pow(nz,2)*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4))) + 
      pow(ny,2)*(-2*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         3*nx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) - 
      3*nx*pow(nz,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
         5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))),
   2*(-(ty*(2*tz*pow(nx,3)*(-15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 5*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
           30*nz*tx*pow(nx,2)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
           3*nx*tz*pow(nz,2)*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
           tx*pow(nz,3)*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)))) + 
      9*ty*pow(ny,2)*(nx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + nz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
      pow(ny,3)*(3*tz*pow(tx,5) - 10*pow(tx,3)*pow(tz,3) + 3*tx*pow(tz,5)) + 
      3*ny*(2*tx*tz*pow(nx,2)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
         tx*tz*pow(nz,2)*(-3*pow(tx,4) + 30*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 9*pow(tz,4)) + 
         3*nx*nz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   3*(2*ty*pow(ny,2)*(nz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         10*nx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      2*ty*(nz*tz*pow(nx,2)*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         tx*pow(nx,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*pow(nz,3)*(5*pow(tx,4) + pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
         nx*tx*pow(nz,2)*(-10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) + 
      ny*(4*nx*nz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
         3*pow(nx,2)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
            5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
         3*pow(nz,2)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
            5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))) + 
      pow(ny,3)*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   3*(2*nx*tx*tz*pow(nz,2)*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
      2*ty*tz*pow(ny,3)*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 2*tx*tz*pow(nx,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
      3*nz*pow(nx,2)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      2*ny*ty*(tz*pow(nz,2)*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         20*nx*nz*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         2*pow(nx,2)*(15*tz*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,3) + 3*pow(tz,5))) + 
      pow(ny,2)*(4*nx*tx*tz*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
         3*nz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      pow(nz,3)*(-5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 5*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) - 
         5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 2*pow(tz,6))),
   -18*nz*tx*tz*pow(nx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
    20*tx*ty*pow(ny,3)*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    4*tx*tz*pow(nz,3)*(15*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 45*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
    3*pow(ny,2)*(2*nz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       3*nx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    6*ny*ty*(-3*tx*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       2*nx*nz*tz*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       tx*pow(nz,2)*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) - 
    3*nx*pow(nz,2)*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       4*pow(tz,6)) + pow(nx,3)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   2*(2*tx*tz*pow(ny,3)*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
      6*ty*pow(ny,2)*(nx*tz*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         5*nz*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
      3*ny*(tx*tz*pow(nz,2)*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
         3*tx*tz*pow(nx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         3*nx*nz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
      ty*(tz*pow(nx,3)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 9*nz*tx*pow(nx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
         3*nx*tz*pow(nz,2)*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
         tx*pow(nz,3)*(-10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))))),
   3*(pow(ny,3)*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
      2*ty*pow(ny,2)*(nz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         3*nx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      2*nz*ty*(2*tz*pow(nz,2)*(pow(ty,4) + 5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
         3*nx*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + pow(nx,2)*(-3*tz*pow(ty,4) + 10*pow(ty,2)*pow(tz,3) - 3*pow(tz,5))) - 
      ny*(12*nx*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         pow(nz,2)*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
            4*pow(tz,6)) + pow(nx,2)*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   4*ty*tz*pow(ny,3)*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    9*pow(ny,2)*(2*nx*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       nz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    6*ny*ty*(tz*pow(nz,2)*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       6*nx*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + pow(nx,2)*(3*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + 3*pow(tz,5))) - 
    nz*(18*nx*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       pow(nz,2)*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          4*pow(tz,6)) + 3*pow(nx,2)*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   3*(2*tx*ty*pow(ny,3)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
      2*ny*nz*ty*(2*nx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      pow(ny,2)*(-6*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      pow(nz,2)*(2*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nx*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   6*(tx*tz*pow(ny,3)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
      ty*pow(nz,2)*(nx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      ty*pow(ny,2)*(nx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 3*nz*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
      ny*nz*(-3*nz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
         nx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   -6*nz*ty*tz*pow(ny,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    2*ty*tz*pow(nz,3)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(ny,3)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    3*ny*pow(nz,2)*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)),
   -6*ny*ty*tz*pow(nz,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
    pow(ny,3)*(6*tz*pow(ty,5) - 20*pow(ty,3)*pow(tz,3) + 6*ty*pow(tz,5)) + 
    3*nz*pow(ny,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    pow(nz,3)*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)),
   pow(nx,3)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    3*nx*pow(nz,2)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    3*nz*pow(nx,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    pow(nz,3)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   pow(nx,3)*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(nx,2)*(ny*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
       nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) + 
    pow(nz,2)*(-3*ny*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) - 
    3*nx*nz*(2*ny*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   pow(nz,3)*(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    3*nz*pow(nx,2)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    pow(nx,3)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
    3*nx*pow(nz,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))),
   pow(nx,3)*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    3*nx*(-(pow(nz,2)*(rx*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
            5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
            2*rz*tz*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) + 
       pow(ny,2)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
       2*ny*nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) + 
    nz*(-3*pow(ny,2)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) - 
       3*ny*nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
       2*pow(nz,2)*(rz*tx*(pow(tx,4) + 5*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 10*pow(tz,4)) + 
          tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))))\
     + 3*pow(nx,2)*(ny*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
       nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 
                3*pow(tz,4))))),pow(nx,3)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*nz*(-2*ny*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) + 
    3*pow(nx,2)*(ny*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))) - 
    pow(nz,2)*(3*ny*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   pow(ny,3)*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    10*pow(nx,3)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    2*pow(nz,3)*(5*rz*ty*(pow(tx,4) + pow(tx,2)*(pow(ty,2) - 9*pow(tz,2)) - pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
       tz*(10*rx*tx*ty*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
    3*nz*pow(nx,2)*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    3*nx*pow(nz,2)*(ry*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       5*ty*(-4*rz*tx*tz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    3*pow(ny,2)*(-(nz*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
            tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))))) + 
       nx*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*(pow(nz,2)*(rx*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
          2*rz*tz*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
       pow(nx,2)*(-5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       2*nx*nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 
                3*pow(tz,4))))),3*nz*pow(nx,2)*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 
          3*pow(tz,4)) + 5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    pow(nz,3)*(rx*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
       2*rz*tz*(5*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))) + 
    2*pow(nx,3)*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*pow(nz,2)*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))
      + 3*pow(ny,2)*(nz*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       nx*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)))) + 
    3*ny*(pow(nx,2)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
       pow(nz,2)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
       2*nx*nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))),
   pow(nx,3)*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nx*pow(nz,2)*(5*rx*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
       rz*tz*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    pow(ny,3)*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
       5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nz*pow(nx,2)*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    pow(nz,3)*(tz*(20*ry*tx*ty*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
          rx*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
       rz*tx*(pow(tx,4) + 10*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) + 5*(pow(ty,4) - 18*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)))) - 
    3*pow(ny,2)*(nx*(-5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) - 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       nz*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 
                3*pow(tz,4))))) - 3*ny*(-10*pow(nx,2)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       2*nx*nz*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
       pow(nz,2)*(ry*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          5*ty*(-4*rz*tx*tz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
             rx*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   pow(ny,3)*(rx*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
    2*pow(nx,3)*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    30*nz*pow(nx,2)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*pow(nz,2)*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
    pow(nz,3)*(ry*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
       5*ty*(-4*rz*tx*tz*(2*pow(tx,2) + pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    3*pow(ny,2)*(nx*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(ry*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4)) + 
          5*ty*(4*rz*tx*tz*(-pow(tx,2) + pow(tz,2)) + rx*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*(-2*nx*nz*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
          5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))
         - 2*pow(nx,2)*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       pow(nz,2)*(rz*tx*(pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 15*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2))) + 
          tz*(20*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 
                3*pow(tz,4))))),pow(ny,3)*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
       5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*nz*pow(nx,2)*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))\
     + pow(nx,3)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) + 
    pow(nz,3)*(rz*ty*(5*pow(tx,4) + pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - 9*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 30*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
          ry*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)))) - 
    3*nx*pow(nz,2)*(5*ry*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ty*(-20*rz*tx*tz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) - 
    3*ny*(pow(nx,2)*(-5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       pow(nz,2)*(5*rx*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          5*ry*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
          rz*tz*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
       2*nx*nz*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))) - 
    3*pow(ny,2)*(-10*nx*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   3*nz*pow(nx,2)*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
    pow(nz,3)*(-5*rx*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
       5*ry*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       rz*tz*(5*pow(tx,4) + 5*pow(ty,4) + 30*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4))) + 
    pow(ny,3)*(5*rz*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(tx,2) - pow(tz,2)) + ry*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)))) + 
    pow(nx,3)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*pow(nz,2)*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    3*pow(ny,2)*(nz*(rz*tz*(-5*pow(tx,4) + 10*pow(ty,2)*pow(tz,2) + pow(tx,2)*(-30*pow(ty,2) + 20*pow(tz,2)) - 3*pow(tz,4)) + 
          5*ry*ty*(pow(tx,4) - 6*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 10*rx*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)))
         + 2*nx*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*(-2*pow(nx,2)*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
       20*nx*nz*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       pow(nz,2)*(5*rz*ty*(pow(tx,4) + 2*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 2*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             ry*(5*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   -3*nx*pow(nz,2)*(5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       2*rz*tz*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
       ry*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
    pow(nx,3)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    3*nz*pow(nx,2)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    10*pow(ny,3)*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    2*pow(nz,3)*(5*rz*tx*(pow(ty,4) + pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 9*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)))) - 
    3*pow(ny,2)*(nx*(-5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          rz*tz*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))) - 
    3*ny*(2*nx*nz*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))
          ) - pow(nx,2)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) + 
       pow(nz,2)*(5*ry*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ty*(-20*rz*tx*tz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
             rx*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))))),
   pow(nx,3)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    2*pow(ny,3)*(5*rz*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*ry*tx*ty*(pow(tx,2) - pow(tz,2)) + rx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*nx*pow(nz,2)*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
       tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))\
     + 3*nz*pow(nx,2)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    pow(nz,3)*(5*ry*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ty*(-20*rz*tx*tz*(pow(tx,2) + 2*pow(ty,2) - 3*pow(tz,2)) + 
          rx*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) + 
    6*pow(ny,2)*(nx*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       5*nz*(ry*tx*(pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(-2*rz*tx*tz*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + rx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*(-2*nx*nz*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
          5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))
         - pow(nx,2)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       pow(nz,2)*(5*rz*tx*(pow(ty,4) + 2*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(tx,2) + pow(ty,2) - 2*pow(tz,2)) + 
             rx*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))))),
   pow(ny,3)*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
       5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4))) - 
    3*ny*(pow(nz,2)*(5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          2*rz*tz*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
          ry*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
       pow(nx,2)*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       2*nx*nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    3*pow(ny,2)*(nz*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))
          ) - nx*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))) + 
    nz*(-3*pow(nx,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       2*pow(nz,2)*(rz*ty*(pow(ty,4) + 5*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 10*pow(tz,4)) + 
          tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4))))
         - 3*nx*nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   2*pow(ny,3)*(5*rz*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(10*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    3*pow(ny,2)*(nz*(rz*tz*(-5*pow(ty,4) + 20*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 3*pow(tz,4)) + 
          5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 10*ry*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - pow(ty,2)*pow(tz,2) + pow(tz,4)))
         + nx*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    nz*(pow(nz,2)*(5*rx*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          2*rz*tz*(5*pow(ty,4) + 5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 15*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
          ry*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
       3*pow(nx,2)*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       3*nx*nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*(-(pow(nx,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4)))) + 
       pow(nz,2)*(rz*ty*(pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)) + 
          tz*(20*rx*tx*ty*(pow(ty,2) - pow(tz,2)) + ry*(5*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))
          ) - 2*nx*nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(ny,3)*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    3*pow(ny,2)*(nx*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) + 
    pow(nz,2)*(3*nx*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) - 
    3*ny*nz*(2*nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(ny,3)*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    3*ny*nz*(2*nx*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*rz*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          tz*(20*ry*tx*ty*(pow(ty,2) - pow(tz,2)) + rx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))))) + 
    3*pow(ny,2)*(nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))) - 
    pow(nz,2)*(3*nx*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
       nz*(5*ry*tx*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ty*(20*rz*tx*tz*(-pow(ty,2) + pow(tz,2)) + rx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))),
   pow(ny,3)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    3*ny*pow(nz,2)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    3*nz*pow(ny,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    pow(nz,3)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   pow(nz,3)*(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    3*nz*pow(ny,2)*(-(rz*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + ry*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) + 
    pow(ny,3)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))) - 
    3*ny*pow(nz,2)*(ry*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rz*(pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 5*ty*pow(tz,4))),
   -(tx*pow(nz,2)*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6))) + 
    2*nx*nz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    pow(nx,2)*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   2*nx*tx*(-14*nz*ty*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
       ny*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6))) + 
    7*ty*pow(nx,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    nz*(-7*nz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
       2*ny*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))),
   2*nx*nz*tx*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)) + 
    tz*pow(nz,2)*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    pow(nx,2)*(7*tz*pow(tx,6) - 35*pow(tx,4)*pow(tz,3) + 21*pow(tx,2)*pow(tz,5) - pow(tz,7)),
   -(tx*pow(nz,2)*(pow(tx,6) + 21*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) - 105*pow(tx,2)*(2*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         7*(15*pow(ty,2) - 4*pow(tz,2))*pow(tz,4))) + 14*ny*ty*
     (-2*nz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       nx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    7*tx*pow(nx,2)*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
       3*pow(tz,6)) + 2*nx*nz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
       21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6)) + pow(ny,2)*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)),
   2*(7*ty*(tx*tz*pow(nz,2)*(-3*pow(tx,4) + 10*pow(tx,2)*pow(tz,2) - 3*pow(tz,4)) + 
         pow(nx,2)*(3*tz*pow(tx,5) - 10*pow(tx,3)*pow(tz,3) + 3*tx*pow(tz,5)) + 
         nx*nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
      ny*(-(nx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))) + 
         nz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6)))),
   7*ty*pow(ny,2)*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) - 
    7*ty*(4*nx*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       pow(nz,2)*(pow(tx,6) + 5*pow(tx,4)*(pow(ty,2) - 6*pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) - 
          4*pow(tz,6)) + pow(nx,2)*(-5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2))*pow(tz,2) - 
          5*pow(ty,2)*pow(tz,4) + 3*pow(tz,6))) + 2*ny*(7*nx*tx*
        (3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       nz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
          21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6))),14*ny*ty*(2*nx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       nz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    14*nx*nz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    tz*pow(nz,2)*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
       21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6)) + pow(nx,2)*(35*pow(tx,4)*(3*tz*pow(ty,2) - pow(tz,3)) - 
       42*pow(tx,2)*(5*pow(ty,2)*pow(tz,3) - pow(tz,5)) + 21*pow(ty,2)*pow(tz,5) - 3*pow(tz,7)) + 
    pow(ny,2)*(7*tz*pow(tx,6) - 35*pow(tx,4)*pow(tz,3) + 21*pow(tx,2)*pow(tz,5) - pow(tz,7)),
   14*ny*ty*(-2*nz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       nx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) - 
    7*tx*pow(nz,2)*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 15*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 
       5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 6*pow(tz,6)) + 
    7*tx*pow(ny,2)*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
       3*pow(tz,6)) + 2*nx*nz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
       21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
    7*tx*pow(nx,2)*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   2*(7*tx*ty*tz*pow(ny,2)*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
      7*ty*(tx*tz*pow(nz,2)*(-3*pow(tx,4) - 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         2*tx*tz*pow(nx,2)*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         nx*nz*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
      ny*(7*nz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
         nx*tz*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)))),
   7*ty*(-4*nx*nz*tx*tz*(3*pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       pow(nz,2)*(-5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 3*pow(tz,2)*(pow(ty,4) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) - 
          3*pow(tx,2)*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
       pow(nx,2)*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
    7*ty*pow(ny,2)*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
    2*ny*(nz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
          21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
       7*nx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   14*ny*ty*(4*nx*tx*tz*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nz*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    tz*pow(nz,2)*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
       21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
    14*nx*nz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
    pow(ny,2)*(35*pow(tx,4)*(3*tz*pow(ty,2) - pow(tz,3)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,3) - pow(tz,5)) + 21*pow(ty,2)*pow(tz,5) - 
       3*pow(tz,7)) + pow(nx,2)*(-35*pow(ty,4)*pow(tz,3) + 42*pow(ty,2)*pow(tz,5) + 
       21*pow(tx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) - 3*pow(tz,7)),
   14*ny*ty*(2*nz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
       nx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) - 
    7*tx*pow(nz,2)*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       4*pow(tz,6)) + 7*tx*pow(nx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    2*nx*nz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       4*pow(tz,6)) + 7*tx*pow(ny,2)*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
       3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   2*(14*tx*ty*tz*pow(ny,2)*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      7*ty*(tx*tz*pow(nz,2)*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
         tx*tz*pow(nx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         nx*nz*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
      ny*(nx*tz*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)) + 
         7*nz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6))))),
   -(ty*(28*nx*nz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         pow(nz,2)*(pow(ty,6) - 42*pow(ty,4)*pow(tz,2) + 105*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
            28*pow(tz,6)) - pow(nx,2)*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)))) + 
    7*ty*pow(ny,2)*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 
       3*pow(tz,6)) + 2*ny*(7*nx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       nz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          4*pow(tz,6))),14*ny*ty*(2*nx*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nz*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
    14*nx*nz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    tz*pow(nz,2)*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
       4*pow(tz,6)) + pow(ny,2)*(-35*pow(ty,4)*pow(tz,3) + 42*pow(ty,2)*pow(tz,5) + 
       21*pow(tx,2)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) - 3*pow(tz,7)) + 
    pow(nx,2)*(7*tz*pow(ty,6) - 35*pow(ty,4)*pow(tz,3) + 21*pow(ty,2)*pow(tz,5) - pow(tz,7)),
   2*ny*ty*(-14*nz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       nx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
    7*tx*pow(ny,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    nz*(-7*nz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       2*nx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   2*(7*tx*ty*tz*pow(ny,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
      nz*ty*(-7*nz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         nx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
      ny*(7*nz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
         nx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   -(ty*pow(nz,2)*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
    2*ny*nz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    pow(ny,2)*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6)),
   2*ny*nz*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) + 
    tz*pow(nz,2)*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    pow(ny,2)*(7*tz*pow(ty,6) - 35*pow(ty,4)*pow(tz,3) + 21*pow(ty,2)*pow(tz,5) - pow(tz,7)),
   pow(nx,2)*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) - 
    2*nx*nz*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nz,2)*(2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6))),
   2*nx*(-2*nz*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
          tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
       ny*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
          rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))) + 
    pow(nx,2)*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
       ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(-2*ny*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
          rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       nz*(6*ty*(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) - rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
          ry*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6)))),
   2*nx*nz*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nx,2)*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nz,2)*(-2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(-pow(tx,6) + 15*pow(tx,4)*pow(tz,2) - 15*pow(tx,2)*pow(tz,4) + pow(tz,6))),
   2*ny*(-2*nz*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
          tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
       nx*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
          ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))) - 
    pow(nz,2)*(6*tx*(-2*rz*tz*(pow(tx,4) + 5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
          ry*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4))) + 
       rx*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
          4*pow(tz,6))) - 2*nx*nz*(2*tz*(3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          rx*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
          4*pow(tz,6))) + pow(ny,2)*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nx,2)*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*(pow(nx,2)*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) - 
      pow(nz,2)*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
      ny*(nz*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
            rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
         nx*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
            rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))) + 
      nx*nz*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
         ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))),
   -4*nx*nz*(rz*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) + 
       tz*(ry*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)))) - 
    2*ny*(nz*(2*tz*(3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
             rx*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
          rz*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
             4*pow(tz,6))) + nx*(-6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
          2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) - 
          3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))) - 
    pow(nz,2)*(2*ty*(-2*rz*tz*(15*pow(tx,4) + 15*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
          rx*tx*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4))) + 
       ry*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
          4*pow(tz,6))) + pow(ny,2)*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
          rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))\
     + pow(nx,2)*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*ny*(2*nx*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
          tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
       nz*(6*ty*(-(rz*tz*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + rx*(pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 5*tx*pow(tz,4))) + 
          ry*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)))) - 
    pow(nz,2)*(2*tz*(3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + 
          rx*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*(pow(tx,6) + 15*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 45*pow(tx,2)*pow(tz,2)*(-2*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 
          4*pow(tz,6))) + pow(ny,2)*(2*rx*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
       rz*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    2*nx*nz*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nx,2)*(2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   pow(nx,2)*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    pow(nz,2)*(-2*ry*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) - 
       3*(-2*rz*tx*tz*(pow(tx,4) + 5*pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
          rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
             5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6)))) - 
    2*ny*(2*nz*(rz*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) + 
          tz*(ry*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
             rx*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)))) + 
       nx*(2*ty*(-10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
             rz*tz*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) - 
          3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))) - 
    2*nx*nz*(2*tz*(rx*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          ry*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
          5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))) + 
    pow(ny,2)*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
       2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
       3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*(2*pow(nx,2)*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
            rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) - 
      pow(nz,2)*(rz*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 6*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) + 
         tz*(ry*tx*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
            rx*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)))) + 
      pow(ny,2)*(3*rz*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*(3*rx*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4)) + ry*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)))) + 
      ny*(nz*(6*ry*tx*ty*(pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 5*pow(tz,4)) - 
            2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(3*pow(ty,2) - 2*pow(tz,2)) - 30*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
            3*rx*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
         nx*(2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
               3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
            3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))) + 
      nx*nz*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
            10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
         3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))),
   pow(nx,2)*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    4*nx*nz*(tz*(ry*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*tx*ty*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) - 
    2*ny*(nx*(-20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          2*rz*tx*tz*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) - 
          3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
       nz*(2*tz*(rx*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
             ry*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
          3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
             5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6)))) + 
    pow(nz,2)*(2*ty*(3*rz*tz*(5*pow(tx,4) + pow(ty,4) + 10*pow(tx,2)*(2*pow(ty,2) - 3*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 6*pow(tz,4)) + 
          rx*tx*(-10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) - 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
          5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))) + 
    pow(ny,2)*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
       3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*nx*nz*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
    pow(nx,2)*(3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
          5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    2*ny*(4*nx*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
             rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
       nz*(2*ty*(rz*tz*(-15*pow(tx,4) - 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) + 10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
             10*rx*tx*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4))) + 
          3*ry*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6)))) - 
    pow(nz,2)*(2*tz*(rx*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          ry*ty*(15*pow(tx,4) + 30*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) + 
          5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 2*pow(tz,6))) + 
    pow(ny,2)*(2*tz*(2*rx*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
          3*ry*ty*(5*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + pow(tz,4))) + 
       3*rz*(5*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   pow(ny,2)*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
       20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
       3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) - 
    2*ny*(nx*(-3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
             5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
          2*ty*(-3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
             rz*tz*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)))) + 
       2*nz*(tz*(ry*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
             rx*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
          rz*tx*ty*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))))) - 
    pow(nz,2)*(4*rz*tx*tz*(-15*pow(ty,4) + 45*pow(ty,2)*pow(tz,2) + 5*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 6*pow(tz,4)) + 
       2*ry*tx*ty*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
       rx*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))
      - 2*nx*nz*(2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ry*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))
      + pow(nx,2)*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   2*(pow(nx,2)*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      2*pow(ny,2)*(5*rz*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
         tz*(ry*tx*(5*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) + 3*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2))) + 
            rx*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      nx*nz*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
         2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
            3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
      pow(nz,2)*(tz*(ry*tx*(15*pow(ty,4) + 10*pow(tx,2)*(3*pow(ty,2) - pow(tz,2)) - 60*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
            rx*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
         rz*tx*ty*(10*pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) + 3*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4)))) + 
      ny*(nz*(2*rz*tx*tz*(-15*pow(ty,4) + 60*pow(ty,2)*pow(tz,2) + 10*pow(tx,2)*(-3*pow(ty,2) + pow(tz,2)) - 9*pow(tz,4)) + 
            20*ry*tx*ty*(pow(tx,2)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
            3*rx*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)))) + 
         nx*(3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
            2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
               2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))))),
   -4*nx*nz*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    pow(ny,2)*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
          5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)))) - 
    pow(nz,2)*(6*ty*(-2*rz*tz*(pow(ty,4) + 5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
          rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))) + 
       ry*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))
      + pow(nx,2)*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) - 
    2*ny*(nz*(2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
             ry*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
          rz*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
             4*pow(tz,6))) + nx*(6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
          6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
          rx*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   pow(ny,2)*(3*rz*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
       2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          2*ry*ty*(15*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
    2*ny*(2*nx*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
          tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
       nz*(3*ry*(-(pow(tz,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4))) + 
          2*ty*(rz*tz*(-3*pow(ty,4) - 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
             3*rx*tx*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4))))) - 
    pow(nz,2)*(2*tz*(3*rx*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          ry*ty*(3*pow(ty,4) + 30*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4))) + 
       rz*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 15*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))
      + 2*nx*nz*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nx,2)*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))),
   pow(ny,2)*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
       rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(-2*nx*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
          rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
       nz*(6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
          rx*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)))) - 
    2*ny*(nz*(6*ry*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
          2*rx*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4))) + 
       nx*(2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
          ry*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   2*(pow(ny,2)*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
         tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
      ny*(nz*(-6*rz*tx*tz*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 6*ry*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
            rx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
         nx*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
            rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)))) - 
      nz*(nz*(3*rz*tx*ty*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) + 
            tz*(3*ry*tx*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + rx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)))) + 
         nx*(2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
            ry*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6))))),
   pow(ny,2)*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) - 
    2*ny*nz*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nz,2)*(2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   2*ny*nz*(-2*rz*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       ry*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    pow(ny,2)*(2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    pow(nz,2)*(-2*ry*tz*(3*pow(ty,5) - 10*pow(ty,3)*pow(tz,2) + 3*ty*pow(tz,4)) + 
       rz*(-pow(ty,6) + 15*pow(ty,4)*pow(tz,2) - 15*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   8*nz*tx*tz*(-pow(tx,6) + 7*pow(tx,4)*pow(tz,2) - 7*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    nx*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   8*ty*(nz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
       nx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    ny*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   8*nx*tx*tz*(pow(tx,6) - 7*pow(tx,4)*pow(tz,2) + 7*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    nz*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   4*(2*ny*tx*ty*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)) - 
      2*nz*tx*tz*(pow(tx,6) + 7*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) + 
         pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 21*pow(tz,4)) - 4*pow(tz,6)) + 
      nx*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   8*(nz*tx*ty*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)) + 
      tz*(nx*ty*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
         ny*(pow(tx,7) - 7*pow(tx,5)*pow(tz,2) + 7*pow(tx,3)*pow(tz,4) - tx*pow(tz,6)))),
   4*(2*ty*(7*nx*tx*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
         nz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 7*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) - 7*pow(ty,2)*pow(tz,4) + 
            4*pow(tz,6))) + ny*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   4*(2*tz*(nx*tx*(7*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 14*pow(tx,2)*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) - 
            3*pow(tz,6)) + ny*ty*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
      nz*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   2*(28*ny*tx*ty*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
      4*nz*tx*tz*(35*pow(ty,4)*pow(tz,2) + 7*pow(tx,4)*(-3*pow(ty,2) + pow(tz,2)) - 63*pow(ty,2)*pow(tz,4) - 
         7*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
      nx*(35*pow(ty,4)*pow(tz,4) + 35*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 42*pow(ty,2)*pow(tz,6) - 
         42*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 3*pow(tz,8))),
   8*(tz*(ny*tx*(7*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 14*pow(tx,2)*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
         nx*ty*(35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 3*pow(tz,6))) + 
      7*nz*tx*ty*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))),
   2*(4*ty*(7*nx*tx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
         nz*tz*(-35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) - 
            7*pow(tx,2)*(3*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 6*pow(tz,6))) + 
      ny*(35*pow(ty,4)*pow(tz,4) + 35*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 42*pow(ty,2)*pow(tz,6) - 
         42*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 3*pow(tz,8))),
   8*tz*(nx*tx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)) + 
       ny*ty*(35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 3*pow(tz,6))) + 
    nz*(70*pow(ty,4)*pow(tz,4) + 70*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 84*pow(ty,2)*pow(tz,6) - 
       84*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 6*pow(tz,8)),
   4*(-2*nz*tx*tz*(7*pow(ty,6) - 70*pow(ty,4)*pow(tz,2) + 63*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         4*pow(tz,6)) + 14*ny*tx*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 
         pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
      nx*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
         21*pow(ty,2)*pow(tz,6) + pow(tz,8))),8*(tz*(ny*tx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 
            7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)) + 
         nx*ty*(-7*pow(ty,4)*pow(tz,2) + 14*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 3*pow(tz,6))) + 
      7*nz*tx*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))),
   4*(2*ty*(nx*tx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) + 
         nz*tz*(-pow(ty,6) + 14*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) - 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
            4*pow(tz,6))) + ny*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 
         7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 21*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   4*(2*tz*(ny*ty*(-7*pow(ty,4)*pow(tz,2) + 14*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
            3*pow(tz,6)) + nx*tx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      nz*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
         21*pow(ty,2)*pow(tz,6) + pow(tz,8))),8*ny*tx*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) + 
    8*nz*tx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    nx*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   8*(tz*(nx*ty*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
         ny*tx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      nz*tx*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))),
   8*nz*ty*tz*(-pow(ty,6) + 7*pow(ty,4)*pow(tz,2) - 7*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    ny*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   8*ny*ty*tz*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    nz*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   nx*(rz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
       rx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    nz*(rx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) - 
       rz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))),
   ny*(rz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
       rx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    nx*(7*ty*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
          rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       ry*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    nz*(-14*rx*tx*ty*tz*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) - 
       7*rz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
       ry*(-7*tz*pow(tx,6) + 35*pow(tx,4)*pow(tz,3) - 21*pow(tx,2)*pow(tz,5) + pow(tz,7))),
   nz*(rz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
       rx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    nx*(-(rx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))) + 
       rz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))),
   -(nz*(rz*tx*(pow(tx,6) + 21*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) - 105*pow(tx,2)*(2*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
            7*(15*pow(ty,2) - 4*pow(tz,2))*pow(tz,4)) + tz*(14*ry*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
            rx*(7*pow(tx,6) + 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) + 
               pow(tx,2)*(-210*pow(ty,2)*pow(tz,2) + 63*pow(tz,4)) - 4*pow(tz,6))))) + 
    nx*(7*rx*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       7*ry*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
          21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6))) + ny*(7*ty*
        (-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
          rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       ry*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))),
   nx*(tz*(14*rx*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       7*rz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(7*ty*(-2*rz*tz*(3*pow(tx,5) - 10*pow(tx,3)*pow(tz,2) + 3*tx*pow(tz,4)) + 
          rx*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       ry*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    ny*(-(rx*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6))) + 
       rz*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))),
   -(nz*(tz*(14*rx*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
            ry*(7*pow(tx,6) + 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) + 
               pow(tx,2)*(-210*pow(ty,2)*pow(tz,2) + 63*pow(tz,4)) - 4*pow(tz,6))) + 
         7*rz*ty*(pow(tx,6) + 5*pow(tx,4)*(pow(ty,2) - 6*pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-30*pow(ty,2)*pow(tz,2) + 45*pow(tz,4)) - 
            4*pow(tz,6)))) + 7*nx*(ty*(-2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       ry*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    ny*(7*rx*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       7*ry*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
          21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6))),nx*(tz*(14*ry*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       7*rz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    ny*(tz*(14*rx*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
       7*rz*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(7*rx*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       7*ry*ty*(pow(tx,6) - 15*pow(tx,4)*pow(tz,2) + 15*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 3*pow(tz,4)) - 
          21*pow(ty,2)*pow(tz,4) + 4*pow(tz,6))),7*ny*(ty*(-2*rz*tx*tz*
           (3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       ry*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    nx*(7*ry*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       rz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
          21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
       7*rx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))) + 
    nz*(-7*rz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) - 15*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 
          5*pow(tx,2)*(pow(ty,4) - 12*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 6*pow(tz,6)) + 
       tz*(-14*ry*tx*ty*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
             21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)))),
   nx*(tz*(28*rx*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       7*rz*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    7*nz*(ty*(-2*rz*tx*tz*(3*pow(tx,4) + 10*pow(tx,2)*(pow(ty,2) - 2*pow(tz,2)) - 10*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          rx*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       ry*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    ny*(tz*(14*ry*tx*ty*(3*pow(tx,4) - 10*pow(tx,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       7*rz*tx*(3*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 15*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))),
   ny*(7*ry*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       rz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
          21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
       7*rx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))) + 
    7*nx*(ty*(2*rz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          rx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
       ry*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))) + 
    nz*(-7*rz*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) - 3*pow(tz,2)*(pow(ty,4) - 5*pow(ty,2)*pow(tz,2) + 2*pow(tz,4)) + 
          3*pow(tx,2)*(pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 15*pow(tz,4))) + 
       tz*(-14*rx*tx*ty*(3*pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
          ry*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
             21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)))),
   ny*(tz*(28*rx*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) - 42*pow(tx,2)*(5*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
       7*rz*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))) + 
    nz*(7*ry*ty*(5*pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 30*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
       rz*tz*(-35*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 35*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 
          21*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
       7*rx*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))) + 
    nx*(tz*(28*ry*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
       7*rz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   -(nz*(tz*(14*ry*tx*ty*(3*pow(ty,4) + 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 
            rx*(7*pow(ty,6) - 70*pow(ty,4)*pow(tz,2) + 63*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
               4*pow(tz,6))) + 7*rz*tx*(pow(ty,6) - 30*pow(ty,4)*pow(tz,2) + 45*pow(ty,2)*pow(tz,4) + 
            5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 4*pow(tz,6)))) + 
    nx*(7*ry*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
       7*rx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          4*pow(tz,6))) + 7*ny*(ty*(2*rz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          rx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
       ry*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   nx*(tz*(14*rx*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
       7*rz*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
    7*nz*(ty*(2*rz*tx*tz*(-3*pow(ty,4) - 10*pow(tx,2)*(pow(ty,2) - pow(tz,2)) + 20*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) + 
          rx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
       ry*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))) + 
    ny*(tz*(28*ry*tx*ty*(5*pow(tx,2)*(pow(ty,2) - pow(tz,2)) - 5*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
       7*rz*tx*(5*pow(tx,2)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)))),
   -(nz*(tz*(14*rx*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
            ry*(7*pow(ty,6) - 70*pow(ty,4)*pow(tz,2) + 63*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
               4*pow(tz,6))) + rz*ty*(pow(ty,6) - 42*pow(ty,4)*pow(tz,2) + 105*pow(ty,2)*pow(tz,4) + 
            21*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 28*pow(tz,6)))) + 
    nx*(ty*(-14*rz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
       7*ry*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(7*ry*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
       7*rx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          4*pow(tz,6))),ny*(tz*(14*rx*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          ry*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6))) + 
       7*rz*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))) + 
    nx*(tz*(14*ry*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
       7*rz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(7*ry*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 3*pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
       7*rx*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       rz*tz*(-7*pow(ty,6) + 70*pow(ty,4)*pow(tz,2) - 63*pow(ty,2)*pow(tz,4) - 21*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) + 
          4*pow(tz,6))),ny*(ty*(-14*rz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
       7*ry*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nz*(-14*ry*tx*ty*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
       7*rz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
       rx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
    nx*(rz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
       ry*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))),
   nz*(ty*(-14*rz*tx*tz*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))) + 
       7*ry*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    ny*(tz*(14*ry*tx*ty*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
          rx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
       7*rz*tx*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
    nx*(-(ry*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
       rz*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))),
   ny*(rz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
       ry*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))) + 
    nz*(ry*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) - 
       rz*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))),
   nz*(rz*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
       ry*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))) + 
    ny*(-(ry*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))) + 
       rz*(pow(ty,7) - 21*pow(ty,5)*pow(tz,2) + 35*pow(ty,3)*pow(tz,4) - 7*ty*pow(tz,6))),
   pow(tx,9) - 36*pow(tx,7)*pow(tz,2) + 126*pow(tx,5)*pow(tz,4) - 84*pow(tx,3)*pow(tz,6) + 9*tx*pow(tz,8),
   9*ty*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   9*tz*pow(tx,8) - 84*pow(tx,6)*pow(tz,3) + 126*pow(tx,4)*pow(tz,5) - 36*pow(tx,2)*pow(tz,7) + pow(tz,9),
   36*tx*(pow(tx,6)*(pow(ty,2) - pow(tz,2)) + 7*pow(tx,4)*pow(tz,2)*(-3*pow(ty,2) + pow(tz,2)) + 7*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
      7*pow(ty,2)*pow(tz,6) + pow(tz,8)),72*tx*ty*tz*(pow(tx,6) - 7*pow(tx,4)*pow(tz,2) + 7*pow(tx,2)*pow(tz,4) - pow(tz,6)),
   12*ty*(7*pow(tx,6)*(pow(ty,2) - 3*pow(tz,2)) + 105*pow(tx,4)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 
      21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + 3*pow(tz,8)),
   4*(21*pow(tx,6)*(3*tz*pow(ty,2) - pow(tz,3)) - 63*pow(tx,4)*(5*pow(ty,2)*pow(tz,3) - pow(tz,5)) + 
      27*pow(tx,2)*(7*pow(ty,2)*pow(tz,5) - pow(tz,7)) - 9*pow(ty,2)*pow(tz,7) + pow(tz,9)),
   18*tx*(35*pow(ty,4)*pow(tz,4) + 7*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 42*pow(ty,2)*pow(tz,6) - 
      14*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 3*pow(tz,8)),
   24*tx*ty*tz*(21*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 9*pow(tz,6)),
   18*ty*(7*pow(ty,4)*pow(tz,4) + 7*pow(tx,4)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 14*pow(ty,2)*pow(tz,6) - 
      14*pow(tx,2)*(3*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + 3*pow(tz,6)) + 3*pow(tz,8)),
   6*(21*pow(ty,4)*pow(tz,5) + 21*pow(tx,4)*(5*tz*pow(ty,4) - 10*pow(ty,2)*pow(tz,3) + pow(tz,5)) - 18*pow(ty,2)*pow(tz,7) - 
      6*pow(tx,2)*(35*pow(ty,4)*pow(tz,3) - 42*pow(ty,2)*pow(tz,5) + 3*pow(tz,7)) + pow(tz,9)),
   12*tx*(7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
      3*pow(tz,2)*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6))),
   24*tx*ty*tz*(-21*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 9*pow(tz,6)),
   36*ty*(-(pow(ty,6)*pow(tz,2)) + 7*pow(ty,4)*pow(tz,4) + pow(tx,2)*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) - 
      7*pow(ty,2)*pow(tz,6) + pow(tz,8)),4*(-21*pow(ty,6)*pow(tz,3) + 63*pow(ty,4)*pow(tz,5) + 
      9*pow(tx,2)*(7*tz*pow(ty,6) - 35*pow(ty,4)*pow(tz,3) + 21*pow(ty,2)*pow(tz,5) - pow(tz,7)) - 27*pow(ty,2)*pow(tz,7) + pow(tz,9)),
   9*tx*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   72*tx*ty*tz*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)),
   pow(ty,9) - 36*pow(ty,7)*pow(tz,2) + 126*pow(ty,5)*pow(tz,4) - 84*pow(ty,3)*pow(tz,6) + 9*ty*pow(tz,8),
   9*tz*pow(ty,8) - 84*pow(ty,6)*pow(tz,3) + 126*pow(ty,4)*pow(tz,5) - 36*pow(ty,2)*pow(tz,7) + pow(tz,9),
   8*rz*tx*tz*(-pow(tx,6) + 7*pow(tx,4)*pow(tz,2) - 7*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
    rx*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   8*ty*(rz*tz*(-7*pow(tx,6) + 35*pow(tx,4)*pow(tz,2) - 21*pow(tx,2)*pow(tz,4) + pow(tz,6)) + 
       rx*(pow(tx,7) - 21*pow(tx,5)*pow(tz,2) + 35*pow(tx,3)*pow(tz,4) - 7*tx*pow(tz,6))) + 
    ry*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   8*rx*tx*tz*(pow(tx,6) - 7*pow(tx,4)*pow(tz,2) + 7*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
    rz*(pow(tx,8) - 28*pow(tx,6)*pow(tz,2) + 70*pow(tx,4)*pow(tz,4) - 28*pow(tx,2)*pow(tz,6) + pow(tz,8)),
   4*(2*ry*tx*ty*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)) - 
      2*rz*tx*tz*(pow(tx,6) + 7*pow(tx,4)*(3*pow(ty,2) - 2*pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) + 
         pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 21*pow(tz,4)) - 4*pow(tz,6)) + 
      rx*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   8*(rz*tx*ty*(pow(tx,6) - 21*pow(tx,4)*pow(tz,2) + 35*pow(tx,2)*pow(tz,4) - 7*pow(tz,6)) + 
      tz*(rx*ty*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6)) + 
         ry*(pow(tx,7) - 7*pow(tx,5)*pow(tz,2) + 7*pow(tx,3)*pow(tz,4) - tx*pow(tz,6)))),
   4*(2*ty*(7*rx*tx*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
         rz*tz*(-7*pow(tx,6) - 35*pow(tx,4)*(pow(ty,2) - 2*pow(tz,2)) + 7*pow(tx,2)*(10*pow(ty,2)*pow(tz,2) - 9*pow(tz,4)) - 7*pow(ty,2)*pow(tz,4) + 
            4*pow(tz,6))) + ry*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   4*(2*tz*(rx*tx*(7*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 14*pow(tx,2)*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) - 
            3*pow(tz,6)) + ry*ty*(7*pow(tx,6) - 35*pow(tx,4)*pow(tz,2) + 21*pow(tx,2)*pow(tz,4) - pow(tz,6))) + 
      rz*(7*pow(tx,6)*(pow(ty,2) - pow(tz,2)) - 35*pow(tx,4)*(3*pow(ty,2)*pow(tz,2) - pow(tz,4)) + 
         21*pow(tx,2)*(5*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 7*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   2*(28*ry*tx*ty*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
      4*rz*tx*tz*(35*pow(ty,4)*pow(tz,2) + 7*pow(tx,4)*(-3*pow(ty,2) + pow(tz,2)) - 63*pow(ty,2)*pow(tz,4) - 
         7*pow(tx,2)*(5*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 6*pow(tz,6)) + 
      rx*(35*pow(ty,4)*pow(tz,4) + 35*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 42*pow(ty,2)*pow(tz,6) - 
         42*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 3*pow(tz,8))),
   8*(tz*(ry*tx*(7*pow(tx,4)*(3*pow(ty,2) - pow(tz,2)) + 14*pow(tx,2)*pow(tz,2)*(-5*pow(ty,2) + pow(tz,2)) + 21*pow(ty,2)*pow(tz,4) - 3*pow(tz,6)) + 
         rx*ty*(35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 3*pow(tz,6))) + 
      7*rz*tx*ty*(pow(tx,4)*(pow(ty,2) - 3*pow(tz,2)) + 10*pow(tx,2)*pow(tz,2)*(-pow(ty,2) + pow(tz,2)) + 5*pow(ty,2)*pow(tz,4) - 3*pow(tz,6))),
   2*(4*ty*(7*rx*tx*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
         rz*tz*(-35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) - 
            7*pow(tx,2)*(3*pow(ty,4) - 20*pow(ty,2)*pow(tz,2) + 9*pow(tz,4)) + 6*pow(tz,6))) + 
      ry*(35*pow(ty,4)*pow(tz,4) + 35*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 42*pow(ty,2)*pow(tz,6) - 
         42*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 3*pow(tz,8))),
   8*tz*(rx*tx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)) + 
       ry*ty*(35*pow(tx,4)*(pow(ty,2) - pow(tz,2)) + 7*pow(ty,2)*pow(tz,4) + pow(tx,2)*(-70*pow(ty,2)*pow(tz,2) + 42*pow(tz,4)) - 3*pow(tz,6))) + 
    rz*(70*pow(ty,4)*pow(tz,4) + 70*pow(tx,4)*(pow(ty,4) - 6*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 84*pow(ty,2)*pow(tz,6) - 
       84*pow(tx,2)*(5*pow(ty,4)*pow(tz,2) - 10*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 6*pow(tz,8)),
   4*(-2*rz*tx*tz*(7*pow(ty,6) - 70*pow(ty,4)*pow(tz,2) + 63*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 
         4*pow(tz,6)) + 14*ry*tx*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + 
         pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6)) + 
      rx*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
         21*pow(ty,2)*pow(tz,6) + pow(tz,8))),8*(tz*(ry*tx*(-35*pow(ty,4)*pow(tz,2) + 42*pow(ty,2)*pow(tz,4) + 
            7*pow(tx,2)*(5*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + pow(tz,4)) - 3*pow(tz,6)) + 
         rx*ty*(-7*pow(ty,4)*pow(tz,2) + 14*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 3*pow(tz,6))) + 
      7*rz*tx*ty*(-3*pow(ty,4)*pow(tz,2) + 10*pow(ty,2)*pow(tz,4) + pow(tx,2)*(pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 5*pow(tz,4)) - 3*pow(tz,6))),
   4*(2*ty*(rx*tx*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) + 
         rz*tz*(-pow(ty,6) + 14*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) - 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) + 
            4*pow(tz,6))) + ry*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 
         7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 21*pow(ty,2)*pow(tz,6) + pow(tz,8))),
   4*(2*tz*(ry*ty*(-7*pow(ty,4)*pow(tz,2) + 14*pow(ty,2)*pow(tz,4) + 7*pow(tx,2)*(3*pow(ty,4) - 10*pow(ty,2)*pow(tz,2) + 3*pow(tz,4)) - 
            3*pow(tz,6)) + rx*tx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      rz*(-7*pow(ty,6)*pow(tz,2) + 35*pow(ty,4)*pow(tz,4) + 7*pow(tx,2)*(pow(ty,6) - 15*pow(ty,4)*pow(tz,2) + 15*pow(ty,2)*pow(tz,4) - pow(tz,6)) - 
         21*pow(ty,2)*pow(tz,6) + pow(tz,8))),8*ry*tx*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6)) + 
    8*rz*tx*tz*(-7*pow(ty,6) + 35*pow(ty,4)*pow(tz,2) - 21*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    rx*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   8*(tz*(rx*ty*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
         ry*tx*(7*pow(ty,6) - 35*pow(ty,4)*pow(tz,2) + 21*pow(ty,2)*pow(tz,4) - pow(tz,6))) + 
      rz*tx*ty*(pow(ty,6) - 21*pow(ty,4)*pow(tz,2) + 35*pow(ty,2)*pow(tz,4) - 7*pow(tz,6))),
   8*rz*ty*tz*(-pow(ty,6) + 7*pow(ty,4)*pow(tz,2) - 7*pow(ty,2)*pow(tz,4) + pow(tz,6)) + 
    ry*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8)),
   8*ry*ty*tz*(pow(ty,6) - 7*pow(ty,4)*pow(tz,2) + 7*pow(ty,2)*pow(tz,4) - pow(tz,6)) + 
    rz*(pow(ty,8) - 28*pow(ty,6)*pow(tz,2) + 70*pow(ty,4)*pow(tz,4) - 28*pow(ty,2)*pow(tz,6) + pow(tz,8));


    
  }