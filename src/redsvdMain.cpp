#include <string>
#include <fstream>
#include <sstream>
#include <sys/time.h>

#include "cmdline.h"
#include "redsvd.hpp"

using namespace std;
using namespace Eigen;

namespace {

double getSec(){
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int writeMatrix(const string& fn, const MatrixXf& M){
  FILE* outfp = fopen(fn.c_str(), "wb");
  if (outfp == NULL){
    cerr << "cannot open " << fn << endl;
    return -1;
  }

  for (int i = 0; i < M.rows(); ++i){
    for (int j = 0; j < M.cols(); ++j){
      fprintf(outfp, "%+f ",  M(i, j));
    }
    fprintf(outfp, "\n");
  }

  fclose(outfp);
  return 0;
}

int writeVector(const string& fn, const VectorXf& V){
  FILE* outfp = fopen(fn.c_str(), "wb");
  if (outfp == NULL){
    cerr << "cannot open " << fn << endl;
    return -1;
  }

  for (int i = 0; i < V.rows(); ++i){
    fprintf(outfp, "%+f\n", V(i));
  }

  fclose(outfp);
  return 0;
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

void readLine(const string& line,  const size_t lineN, int formatType, fv_t& fv){
  istringstream is(line);
  if (formatType == 1){
    int label = -1;
    if (!(is >> label)){
      cerr << "cannot read " << line << endl;
      throw; 
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

int readFVs(const char* fn, int formatType, vector<fv_t>& fvs){
  fvs.clear();
  ifstream ifs(fn);
  if (!ifs){
    cerr << "failed to open " << fn << endl;
    return -1;
  }

  size_t lineN = 0;
  for (string line; getline(ifs, line); ){
    fv_t fv;
    readLine(line, lineN++, formatType, fv);
    if (fv.size() == 0) continue;
    fvs.push_back(fv);
  }
  return 0;
}

int SVD_SparseMatrix(const string& input, 
			 int formatType, 
			 int rank,
			 MatrixXf& Uret,
			 VectorXf& Sret, 
			 MatrixXf& Vret){
  SparseMatrix<float, RowMajor> A;
  vector<fv_t> fvs;
  double startSec = getSec();
  cout << "read sparse matrix from " << input << " ... " << flush;
  if (readFVs(input.c_str(), formatType, fvs) == -1){
    return -1;
  }
  double endSec = getSec();
  cout << endSec - startSec << " sec." << endl;
  convertFV2Mat(fvs, A);
  vector<fv_t>().swap(fvs);
  cout << "   rows:\t" << A.rows()     << endl
       << "   cols:\t" << A.cols()     << endl
       << "nonzero:\t" << A.nonZeros() << endl
       << "   rank:\t" << rank         << endl;

  cout << "compute SVD... " << flush;
  startSec = getSec();
  REDSVD::runSVD(A, rank, Uret, Sret, Vret);
  endSec = getSec();
  cout << endSec - startSec << " sec." << endl;
  return 0;
}

int readMatrix(const char* fn, MatrixXf& A){
  ifstream ifs(fn);
  if (!ifs){
    cerr << "failed to open " << fn << endl;
    return -1;
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
  if (rowN == 0) return 0;
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
  return 0;
}

int SVD_DenseMatrix(const string& input, 
			 int formatType, 
			 int rank,
			 MatrixXf& Uret,
			 VectorXf& Sret, 
			 MatrixXf& Vret){
  MatrixXf A;
  double startSec = getSec();
  cout << "read dense matrix from " << input << " ... " << flush;
  if (readMatrix(input.c_str(), A) == -1){
    return -1;
  }
  double endSec = getSec();
  cout << endSec - startSec << " sec." << endl;
  cout << "SVD for a dense matrix" << endl
       << "rows:\t" << A.rows() << endl
       << "cols:\t" << A.cols() << endl
       << "rank:\t" << rank  << endl;

  cout << "compute SVD... " << flush;
  startSec = getSec();
  REDSVD::runSVD(A, rank, Uret, Sret, Vret);
  cout << endSec - startSec << " sec." << endl;
  return 0;
}

void setFooter(cmdline::parser& p){
  p.footer(
   "\n\n"
   "redsvd supports the three format type (one line for each row)\n\n"
   "[format=0] ((colum_id:value)+\\n)+\n"
   "0:1.5 2:-0.3 6:1.0 ... \n\n"
   "[format=1] (label (colum_id:value)+\\n)+ (label will be ignored)\n"
   "+1 0:1.5 2:-0.3 6:1.0 ... \n\n"
   "[format=2] (<value>+\\n)+\n"
   "1.5 0.0 -0.3 0.0 0.0 0.0 1.0 0.0\n\n"
   "Use DenseMatrix for format=2, and use SparseMatrix for others\n\n"
   "Example:\n"
   ">redsvd -i imat -o omat -r 10 -f 2\n"
   "compuate SVD for a matrix in imat and output omat.U omat.V, and omat.S\n"
   "with the 10 largest eigen values/vectors\n" 
  );
}

}



int main(int argc, char* argv[]){
  cmdline::parser p;
  p.add<string>("input",  'i', "input file", true);
  p.add<string>("output", 'o', "output file (.U .S .V will be generated)", true);
  p.add<int>   ("rank",   'r', "rank      ", false, 10);
  p.add<int>   ("format", 'f', "format type (=0, 1, 2) See example. ", false, 0);
  p.set_program_name("redsvd");
  setFooter(p);

  if (argc == 1){
    cerr << p.usage() << endl;
    return 0;
  }
 
  if (p.parse(argc, argv) == 0){
    cerr << "Error:" << p.error() << endl
	 << p.usage() << endl;
    return -1;
  }

  int formatType = p.get<int>("format");
  string input   = p.get<string>("input");
  int rank       = p.get<int>("rank");

  MatrixXf Uret;
  VectorXf Sret;
  MatrixXf Vret;
  if (formatType != 2){
    SVD_SparseMatrix(input, formatType, rank, Uret, Sret, Vret);
  } else if (formatType == 2){
    SVD_DenseMatrix(input, formatType, rank, Uret, Sret, Vret);
  }

  string ofn = p.get<string>("output");
  double startSec = getSec();
  cout << "write matrix to " << ofn << "(.U|.S|.V) ... " << flush;
  if (writeMatrix(ofn + ".U", Uret) == -1) return -1;
  if (writeVector(ofn + ".S", Sret) == -1) return -1;
  if (writeMatrix(ofn + ".V", Vret) == -1) return -1;
  cout << getSec() - startSec << " sec." << endl
       << "finished." << endl;
  
  return 0;
}
