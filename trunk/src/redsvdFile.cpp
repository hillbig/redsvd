#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/time.h>
#include "redsvdFile.hpp"

using namespace std;
using namespace Eigen;

namespace REDSVD{

namespace {

void writeMatrix_(const string& fn, const MatrixXf& M){
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

void writeVector_(const string& fn, const VectorXf& V){
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

void readLine(const string& line,  
	      const size_t lineN, 
	      fv_t& fv){
  istringstream is(line);

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


double getSec(){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

void readMatrix(const std::string& fn, SMatrixXf& A){
  vector<fv_t> fvs;
  ifstream ifs(fn.c_str());
  if (!ifs){
    throw string("failed to open") + fn;
  }

  size_t lineN = 0;
  for (string line; getline(ifs, line); ){
    fv_t fv;
    readLine(line, lineN++, fv);
    if (fv.size() == 0) continue;
    fvs.push_back(fv);
  }
  convertFV2Mat(fvs, A);
}

void readMatrix(const std::string& fn, MatrixXf& A){
  ifstream ifs(fn.c_str());
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

void writeMatrix(const string& fn, const REDSVD::RedSVD& A){
  writeMatrix_(fn + ".U", A.matrixU());
  writeVector_(fn + ".S", A.singularValues());
  writeMatrix_(fn + ".V", A.matrixV());
}

void writeMatrix(const string& fn, const REDSVD::RedPCA& A){
  writeMatrix_(fn + ".pc",    A.principalComponents());
  writeVector_(fn + ".score", A.scores());
}

void writeMatrix(const string& fn, const REDSVD::RedSymEigen& A){
  writeMatrix_(fn + ".evec", A.eigenVectors());
  writeVector_(fn + ".eval", A.eigenValues());
}

}
