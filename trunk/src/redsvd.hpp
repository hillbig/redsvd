#ifndef REDSVD_HPP__
#define REDSVD_HPP__

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace REDSVD {

void sampleTwoGaussian(float& f1, float& f2);
void sampleGaussianMat(Eigen::MatrixXf& x);
void processGramSchmidt(Eigen::MatrixXf& mat);

template <class Mat>
void runSVD(Mat& A, const int origR,
	    Eigen::MatrixXf& Uret, Eigen::VectorXf& Sret, Eigen::MatrixXf& Vret){
  if (A.cols() == 0 || A.rows() == 0) return;
  int r = (origR < A.cols()) ? origR : A.cols();
  r = (r < A.rows()) ? r : A.rows();
  
  // Gaussian Random Matrix for A^T
  Eigen::MatrixXf O(A.rows(), r);
  sampleGaussianMat(O);

  // Compute Sample Matrix of A^T
  Eigen::MatrixXf Y = A.transpose() * O;

  // Orthonormalize Y
  processGramSchmidt(Y);

  // Range(B) = Range(A^T)
  Eigen::MatrixXf B = A * Y;

  // Gaussian Random Matrix
  Eigen::MatrixXf P(B.cols(), r);
  sampleGaussianMat(P);
  
  // Compute Sample Matrix of B
  Eigen::MatrixXf Z = B * P;
  
  // Orthonormalize Z
  processGramSchmidt(Z);

  // Range(C) = Range(B)
  Eigen::MatrixXf C = Z.transpose() * B; 

  Eigen::SVD<Eigen::MatrixXf> svdOfC(C);

  // C = USV^T
  // A = 
  // A = Z * U * S * V^T * Y^T()
  Uret = Z * svdOfC.matrixU();
  Sret = svdOfC.singularValues();
  Vret = Y * svdOfC.matrixV();
}

}

#endif // REDSVD_HPP__
