#ifndef REDSVDFILE_HPP__
#define REDSVDFILE_HPP__

#include <string>
#include <iostream>
#include "redsvd.hpp"

namespace REDSVD{

double getSec();

void readMatrix(const char* fn, const std::string& formatType,  Eigen::SparseMatrix<float, Eigen::RowMajor>& A);
void readMatrix(const char* fn, const std::string& formatType, Eigen::MatrixXf& A);

template<class Mat>
void SVDMatrix(const std::string& inputFileName,
	       const std::string& formatType,
	       int rank,
	       RedSVD& redsvd){
  double startSec = getSec();
  std::cout << "read matrix from " << inputFileName << " ... " << std::flush;
  Mat A;
  readMatrix(inputFileName.c_str(), formatType, A);
  double endSec = getSec();
  std::cout << endSec - startSec << " sec." <<std:: endl;
  std::cout << "rows:\t" << A.rows() << std::endl
	    << "cols:\t" << A.cols() << std::endl
	    << "rank:\t" << rank  << std::endl;

  std::cout << "compute SVD... " << std::flush;
  startSec = getSec();
  redsvd.runSVD(A, rank);
  std::cout << endSec - startSec << " sec." << std::endl;
}

template<class Mat>
void PCAMatrix(const std::string& inputFileName,
	       const std::string& formatType,
	       int rank,
	       RedPCA& redpca){
  double startSec = getSec();
  std::cout << "read matrix from " << inputFileName << " ... " << std::flush;
  Mat A;
  readMatrix(inputFileName.c_str(), formatType, A);
  double endSec = getSec();
  std::cout << endSec - startSec << " sec." <<std:: endl;
  std::cout << "rows:\t" << A.rows() << std::endl
	    << "cols:\t" << A.cols() << std::endl
	    << "rank:\t" << rank  << std::endl;

  std::cout << "compute PCA... " << std::flush;
  startSec = getSec();
  redpca.runPCA(A, rank);
  std::cout << endSec - startSec << " sec." << std::endl;
}


void SVDfromFile(const std::string& inputFileName, 
		 const std::string& outputFileName,
		 const std::string& formatType, 
		int rank);

int PCAfromFile(const std::string& inputFileName, 
		const std::string& outputFileName,
		const std::string& formatType, 
		int rank);

}

#endif // REDSVDFILE_HPP__
