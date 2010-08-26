#ifndef REDSVD_HPP__
#define REDSVD_HPP__

#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace REDSVD {

void sampleGaussianMat(Eigen::MatrixXf& x);
void processGramSchmidt(Eigen::MatrixXf& mat);

class RedSVD {
public:
  RedSVD(){}

  template <class Mat>
  RedSVD(Mat& A, const int rank){
    run(A, rank);
  }

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
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
    // A = Z * U * S * V^T * Y^T()
    matU_ = Z * svdOfC.matrixU();
    matS_ = svdOfC.singularValues();
    matV_ = Y * svdOfC.matrixV();
  }
  
  const Eigen::MatrixXf& matrixU() const {
    return matU_;
  }

  const Eigen::VectorXf& singularValues() const {
    return matS_;
  }

  const Eigen::MatrixXf& matrixV() const {
    return matV_;
  }

private:
  Eigen::MatrixXf matU_;
  Eigen::VectorXf matS_;
  Eigen::MatrixXf matV_;
};

class RedSymEigen {
public:
  RedSymEigen(){}

  template <class Mat>
  RedSymEigen(Mat& A, const int rank){
    run(A, rank);
  }  

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
    r = (r < A.rows()) ? r : A.rows();
    
    // Gaussian Random Matrix
    Eigen::MatrixXf O(A.rows(), r);
    sampleGaussianMat(O);
    
    // Compute Sample Matrix of A
    Eigen::MatrixXf Y = A.transpose() * O;
    
    // Orthonormalize Y
    processGramSchmidt(Y);

    Eigen::MatrixXf B = Y.transpose() * A * Y;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenOfB(B);
    
    eigenValues_ = eigenOfB.eigenvalues();
    eigenVectors_ = Y * eigenOfB.eigenvectors();
  }
  
  const Eigen::MatrixXf& eigenVectors() const {
    return eigenVectors_;
  }

  const Eigen::VectorXf& eigenValues() const {
    return eigenValues_;
  }

private:
  Eigen::VectorXf eigenValues_;
  Eigen::MatrixXf eigenVectors_;
};

class RedPCA {
public:
  RedPCA(){}

  template <class Mat>
  RedPCA(const Mat& A, const int rank) {
    run(A, rank);
  }

  template <class Mat> 
  void run(const Mat& A, const int rank) {
    RedSVD redsvd;
    redsvd.run(A, rank);
    const Eigen::VectorXf& S = redsvd.singularValues();
    principalComponents_ = redsvd.matrixV();
    scores_              = redsvd.matrixU() * S.asDiagonal();
  }

  const Eigen::MatrixXf& principalComponents() const {
    return principalComponents_;
  }

  const Eigen::MatrixXf& scores() const {
    return scores_;
  }

 private:
  Eigen::MatrixXf principalComponents_;
  Eigen::MatrixXf scores_;
};

}

#endif // REDSVD_HPP__
