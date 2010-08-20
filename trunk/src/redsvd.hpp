#ifndef REDSVD_HPP__
#define REDSVD_HPP__

#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

namespace REDSVD {

class RedSVD {
public:
  RedSVD(){}

  template <class Mat>
  void runSVD(Mat& A, const int rank){
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
  void sampleTwoGaussian(float& f1, float& f2);
  void sampleGaussianMat(Eigen::MatrixXf& x);
  void processGramSchmidt(Eigen::MatrixXf& mat);
  
  Eigen::MatrixXf matU_;
  Eigen::VectorXf matS_;
  Eigen::MatrixXf matV_;
};



class RedPCA {
public:
  template <class Mat> 
  void runPCA(const Mat& A, const int rank) {
    RedSVD redsvd;
    redsvd.runSVD(A, rank);
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
