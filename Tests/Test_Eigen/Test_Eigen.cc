// in the following file we test the eigen solver of eigen which is showing errors
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

main()
{
	SparseMatrix<double,Eigen::RowMajor> A(3,3);
	A.coeffRef(0,0) = 1;
	A.coeffRef(1,1) = 1;
        A.coeffRef(2,2) = 2;

        EigenSolver<MatrixXd> ES(A);
        MatrixXd vecs = ES.pseudoEigenvectors();
        VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

       std::cout << vecs << std::endl;
}
