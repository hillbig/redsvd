#include <gtest/gtest.h>
#include "redsvd.hpp"

using namespace std;
using namespace Eigen;

const float EPS = 0.001f;

TEST(redsvd, trivial){
  MatrixXf A;
  REDSVD::RedSVD redsvd;
  redsvd.runSVD(A, 0);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();
  ASSERT_EQ(0, U.rows());
  ASSERT_EQ(0, U.cols());
  ASSERT_EQ(0, V.rows());
  ASSERT_EQ(0, V.cols());
  ASSERT_EQ(0, S.rows());
  ASSERT_EQ(1, S.cols());
}


TEST(redsvd, artificial){
  int rowN = 3;
  int colN = 5;
  MatrixXf A(rowN, colN);
  A << 1,  2,  3,  4,  5,
       6,  7,  8,  9, 10, 
      11, 12, 13, 14, 15;
  
  int r = 2;
  REDSVD::RedSVD redsvd;
  redsvd.runSVD(A, r);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();

  ASSERT_EQ(rowN, U.rows());
  ASSERT_EQ(r,    U.cols());
  ASSERT_EQ(colN, V.rows());
  ASSERT_EQ(r,    V.cols());
  ASSERT_EQ(r,    S.rows());
  ASSERT_EQ(1,    S.cols());

  for (int i = 0; i < U.cols(); ++i){
    ASSERT_NEAR(1.f, U.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, U.col(i).dot(U.col(j)), EPS);
    }
  }

  for (int i = 0; i < V.cols(); ++i){
    ASSERT_NEAR(1.f, V.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, V.col(i).dot(V.col(j)), EPS);
    }
  }

  cout << U * S.asDiagonal() * V.transpose() << endl;
  ASSERT_NEAR(0.f, (A - U * S.asDiagonal() * V.transpose()).norm(), EPS);
}

TEST(redsvd, random){
  int rowN = 20;
  int colN = 70;
  MatrixXf A = MatrixXf::Random(rowN, colN);

  int r = 20;
  REDSVD::RedSVD redsvd;
  redsvd.runSVD(A, r);
  MatrixXf U = redsvd.matrixU();
  VectorXf S = redsvd.singularValues();
  MatrixXf V = redsvd.matrixV();

  ASSERT_EQ(rowN, U.rows());
  ASSERT_EQ(r,    U.cols());
  ASSERT_EQ(colN, V.rows());
  ASSERT_EQ(r,    V.cols());
  ASSERT_EQ(r,    S.rows());
  ASSERT_EQ(1,    S.cols());

  for (int i = 0; i < U.cols(); ++i){
    ASSERT_NEAR(1.f, U.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, U.col(i).dot(U.col(j)), EPS);
    }
  }

  for (int i = 0; i < V.cols(); ++i){
    ASSERT_NEAR(1.f, V.col(i).norm(), EPS);
    for (int j = 0; j < i; ++j){
      ASSERT_NEAR(0.f, V.col(i).dot(V.col(j)), EPS);
    }
  }

  ASSERT_NEAR(0.f, (A - U * S.asDiagonal() * V.transpose()).norm(), EPS);
}


