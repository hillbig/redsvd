#include <string>
#include <fstream>
#include <sstream>

#include "cmdline.h"
#include "redsvdFile.hpp"

using namespace std;

namespace {

void setFooter(cmdline::parser& p){
  p.footer(
   "\n\n"
   "redsvd supports the following format types (one line for each row)\n\n"
   "[format=dense] (<value>+\\n)+\n"
   "[format=sparse] ((colum_id:value)+\\n)+\n"
   "Example:\n"
   ">redsvd -i imat -o omat -r 10 -f dense\n"
   "compuate SVD for a dense matrix in imat and output omat.U omat.V, and omat.S\n"
   "with the 10 largest eigen values/vectors\n" 
   ">redsvd -i imat -o omat -r 3 -f sparse -m PCA\n"
   "compuate PCA for a sparse matrix in imat and output omat.PC omat.SCORE\n"
   "with the 3 largest principal components\n" 
  );
}
}

int main(int argc, char* argv[]){
  cmdline::parser p;
  p.add<string>("input",  'i', "input file", true);
  p.add<string>("output", 'o', "output file's prefix", true);
  p.add<int>   ("rank",   'r', "rank      ", false, 10);
  p.add<string>("format", 'f', "format type (dense|sparse) See example. ", false, "dense");
  p.add<string>("method", 'm', "method (SVD|PCA|SymEigen)", false, "SVD");
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

  string input  = p.get<string>("input");
  string output = p.get<string>("output");
  string format = p.get<string>("format");
  int    rank   = p.get<int>   ("rank");
  string method = p.get<string>("method");
  bool isInputSparse = (format == "dense") ? false : true;

  try {
    if (method == "SVD"){
      if (isInputSparse){
	REDSVD::fileProcess<REDSVD::SMatrixXf, REDSVD::RedSVD>(input, output, rank);
      } else {
	REDSVD::fileProcess<Eigen::MatrixXf, REDSVD::RedSVD>(input, output, rank);
      }
    } else if (method == "PCA"){
      if (isInputSparse){
	REDSVD::fileProcess<REDSVD::SMatrixXf, REDSVD::RedPCA>(input, output, rank);
      } else {
	REDSVD::fileProcess<Eigen::MatrixXf, REDSVD::RedPCA>(input, output, rank);
      }
    } else if (method == "SymEigen"){ 
      if (isInputSparse){
	REDSVD::fileProcess<REDSVD::SMatrixXf, REDSVD::RedSymEigen>(input, output, rank);
      } else {
	REDSVD::fileProcess<Eigen::MatrixXf, REDSVD::RedSymEigen>(input, output, rank);
      }
    } else {
      cerr << "unknown method:" << method << endl;
      return -1;
    }
  } catch (const string& error){
    cerr << "Error: " << error << endl;
    return -1;
  }
  return 0;
}
