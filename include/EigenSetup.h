//#define EIGEN_USE_MKL_ALL
//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE MKL_INT
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// setup relevant vector/matrix datatypes from Eigen

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
typedef SparseMatrix<double,Eigen::RowMajor> Sparse_matrix;
typedef Eigen::Triplet<double> triplet;
typedef Eigen::MatrixXd Full_matrix;

		struct system_matrix
		{
			std::vector<triplet> Row_Col_Value;
			Sparse_matrix matrix;
		};
