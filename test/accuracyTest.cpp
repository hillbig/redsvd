/* 
 *  Copyright (c) 2010 Daisuke Okanohara
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#include <iostream>
#include "../src/redsvd.hpp"

using namespace std;
using namespace Eigen;

void test(int row, int col, int actualRank, int estimateRank){
  cout << row << "\t" << col << "\t" << actualRank << "\t" << estimateRank << endl;
  MatrixXf U = MatrixXf::Random(row, actualRank);
  MatrixXf V=  MatrixXf::Random(col, actualRank);
  REDSVD::processGramSchmidt(U);
  REDSVD::processGramSchmidt(V);
  VectorXf S(actualRank);
  for (int i = 0; i < actualRank; ++i){
    S(i) = pow(0.9f, i);
  }
  MatrixXf A = U * S.asDiagonal() * V.transpose();
  REDSVD::RedSVD svdOfA(A, estimateRank);
  

  for (int i = 0; i < estimateRank; ++i){
    cout << i << "\t" << log(S(i)) << "\t" <<  log(svdOfA.singularValues()(i)) << "\t" 
	 << fabs(U.col(i).dot(svdOfA.matrixU().col(i))) << "\t"
	 << fabs(V.col(i).dot(svdOfA.matrixV().col(i))) << endl;
  }


}

int main(int argc, char* argv[]){
  test(100,  100,  10, 10);
  test(100,  100,  20, 10);
  test(100,  100, 100, 10);

  test(1000,  1000,   10, 10);
  test(1000,  1000,   20, 10);
  test(1000,  1000,  100, 10);
  test(1000,  1000, 1000, 10);
  
  return 0;
}
