#ifndef REDSVDFILE_HPP__
#define REDSVDFILE_HPP__

#include <string>
#include <iostream>
#include "redsvd.hpp"

namespace REDSVD{

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrixXf;

double getSec();

void readMatrix(const std::string& fn, SMatrixXf& A);
void readMatrix(const std::string& fn, Eigen::MatrixXf& A);

void writeMatrix(const std::string& fn, const RedSVD& A);
void writeMatrix(const std::string& fn, const RedPCA& A);
void writeMatrix(const std::string& fn, const RedSymEigen& A);

template<class Mat, class RetMat>
void fileProcess(const std::string& inputFileName,
		 const std::string& outputFileName,
		 int rank){
  double startSec = getSec();
  std::cout << "read matrix from " << inputFileName << " ... " << std::flush;
  Mat A;
  readMatrix(inputFileName.c_str(), A);
  double endSec = getSec();
  std::cout << endSec - startSec << " sec." <<std:: endl;
  std::cout << "rows:\t" << A.rows() << std::endl
	    << "cols:\t" << A.cols() << std::endl
	    << "rank:\t" << rank  << std::endl;

  std::cout << "compute SVD... " << std::flush;
  startSec = getSec();
  RetMat retMat(A, rank);
  std::cout << endSec - startSec << " sec." << std::endl;
  
  startSec = getSec();
  writeMatrix(outputFileName, retMat);
  std::cout << getSec() - startSec << " sec." << std::endl
	    << "finished." << std::endl;
}

}

#endif // REDSVDFILE_HPP__
