#include <iostream>
#include "redsvd.hpp"

using namespace std;
using namespace Eigen;

namespace REDSVD {

namespace {

const float SVD_EPS = 0.0001f;
void sampleTwoGaussian(float& f1, float& f2){
  float v1 = (float)(rand()+1) / ((float)RAND_MAX+1);
  float v2 = (float)(rand()+1) / ((float)RAND_MAX+1);
  float len = sqrt(-2.f * log(v1));
  f1 = len * cos(2.f * M_PI * v2);
  f2 = len * sin(2.f * M_PI * v2);
}

}

void sampleGaussianMat(MatrixXf& mat){
  for (int i = 0; i < mat.rows(); ++i){
    int j = 0;
    for ( ; j+1 < mat.cols(); j += 2){
      float f1, f2;
      sampleTwoGaussian(f1, f2);
      mat(i,j  ) = f1;
      mat(i,j+1) = f2;
    }
    for (; j < mat.cols(); j ++){
      float f1, f2;
      sampleTwoGaussian(f1, f2);
      mat(i, j)  = f1;
    }
  }
} 

void processGramSchmidt(MatrixXf& mat){
  for (int i = 0; i < mat.cols(); ++i){
    for (int j = 0; j < i; ++j){
      float r = mat.col(i).dot(mat.col(j));
      mat.col(i) -= r * mat.col(j);
    }
    float norm = mat.col(i).norm();
    if (norm < SVD_EPS){
      for (int k = i; k < mat.cols(); ++k){
	mat.col(k).setZero();
      } 
    }
    mat.col(i) *= (1.f / norm);
  }
}

}
