#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/time.h>
#include "redsvdFile.hpp"

using namespace std;
using namespace Eigen;

namespace {
void writeMatrix(const string& fn, const MatrixXf& M){
  cout << "write " << fn << endl;
  FILE* outfp = fopen(fn.c_str(), "wb");
  if (outfp == NULL){
    throw string("cannot open ") + fn;
  }

  for (int i = 0; i < M.rows(); ++i){
    for (int j = 0; j < M.cols(); ++j){
      fprintf(outfp, "%+f ",  M(i, j));
    }
    fprintf(outfp, "\n");
  }

  fclose(outfp);
}

void writeVector(const string& fn, const VectorXf& V){
  cout << "write " << fn << endl;
  FILE* outfp = fopen(fn.c_str(), "wb");
  if (outfp == NULL){
    throw string("cannot open ") + fn;
  }

  for (int i = 0; i < V.rows(); ++i){
    fprintf(outfp, "%+f\n", V(i));
  }

  fclose(outfp);
}


typedef vector<pair<int, float> > fv_t;

void convertFV2Mat(const vector<fv_t>& fvs, SparseMatrix<float, RowMajor>& A){
  int maxID = 0;
  size_t nonZeroNum = 0;
  for (size_t i = 0; i < fvs.size(); ++i){
    const fv_t& fv(fvs[i]);
    for (size_t j = 0; j < fv.size(); ++j){
      maxID = max(fv[j].first+1, maxID);
    }
    nonZeroNum += fv.size();
  }
  A.resize(fvs.size(), maxID);
  A.reserve(nonZeroNum);
  for (size_t i = 0; i < fvs.size(); ++i){
    A.startVec(i);
    const fv_t& fv(fvs[i]);
    for (size_t j = 0; j < fv.size(); ++j){
      A.insertBack(i, fv[j].first) = fv[j].second;
    }
  }
  A.finalize();
}

void readLine(const string& line,  
	      const size_t lineN, 
	      const string& formatType, 
	      fv_t& fv){
  istringstream is(line);
  if (formatType == "libsvm"){
    int label = -1;
    if (!(is >> label)){
      throw string("cannot read ") + line;
    }
  }

  int id;
  char sep;
  float val;
  while (is >> id >> sep >> val){
    fv.push_back(make_pair(id, val));
  }
  sort(fv.begin(), fv.end());
  fv.erase(unique(fv.begin(), fv.end()), fv.end());
}

}

namespace REDSVD{
double getSec(){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

void readMatrix(const char* fn, const string& formatType, 
		SparseMatrix<float, RowMajor>& A){
  vector<fv_t> fvs;
  ifstream ifs(fn);
  if (!ifs){
    throw string("failed to open") + fn;
  }

  size_t lineN = 0;
  for (string line; getline(ifs, line); ){
    fv_t fv;
    readLine(line, lineN++, formatType, fv);
    if (fv.size() == 0) continue;
    fvs.push_back(fv);
  }
  convertFV2Mat(fvs, A);
}

void readMatrix(const char* fn, const string& formatType, MatrixXf& A){
  ifstream ifs(fn);
  if (!ifs){
    throw string("failed to open " ) + fn;
  }

  vector< vector<float> > vs;
  for (string line; getline(ifs, line); ){
    istringstream is(line);
    vector<float> v; 
    float val;
    while (is >> val){
      v.push_back(val);
    }
    vs.push_back(v);
  }

  size_t rowN = vs.size();
  if (rowN == 0) return;
  size_t colN = vs[0].size();
  A.resize(rowN, colN);
  
  for (size_t i = 0; i < rowN; ++i){
    if (colN != vs[i].size()){
      cerr << "warning: " << i+1 << "-th row has " << vs[i].size() << " entries. " 
	   << colN << " entries are expected" << endl;
    }
    size_t colNmin = min(colN, vs[i].size());
    for (size_t j = 0; j < colNmin; ++j){
      A(i, j) = vs[i][j];
    }
  }
}


void SVDfromFile(const string& inputFileName, 
		 const string& outputFileName, 
		 const string& formatType, 
		 int rank){
  RedSVD redsvd;
  if (formatType == "dense"){
    SVDMatrix<MatrixXf>(inputFileName, formatType, rank, redsvd);
  } else if (formatType == "sparse" || formatType == "libsvm"){
    SVDMatrix<SparseMatrix<float, RowMajor> >(inputFileName, formatType, rank, redsvd);
  } else {
    throw string("unknown formatType:") + formatType; 
  }
  
  double startSec = getSec();
  cout << "write matrix to " << outputFileName << flush << endl;
  writeMatrix(outputFileName + ".U", redsvd.matrixU());
  writeVector(outputFileName + ".S", redsvd.singularValues());
  writeMatrix(outputFileName + ".V", redsvd.matrixV());
  cout << getSec() - startSec << " sec." << endl
       << "finished." << endl;
}

int PCAfromFile(const string& inputFileName, 
		const string& outputFileName, 
		const string& formatType, 
		int rank){
  RedPCA redpca;
  if (formatType == "dense"){
    PCAMatrix<MatrixXf>(inputFileName, formatType, rank, redpca);
  } else if (formatType == "sparse" || formatType == "libsvm"){
    PCAMatrix<SparseMatrix<float, RowMajor> >(inputFileName, formatType, rank, redpca);
  } else {
    throw string("unknown formatType:") + formatType; 
  }
  
  double startSec = getSec();
  writeMatrix(outputFileName + ".PC",    redpca.principalComponents());
  writeMatrix(outputFileName + ".score", redpca.scores());
  cout << getSec() - startSec << " sec." << endl
       << "finished." << endl;
  return 0;
}


}
