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
   "redsvd supports the three format type (one line for each row)\n\n"
   "[format=sparse] ((colum_id:value)+\\n)+\n"
   "[format=dense] (<value>+\\n)+\n"
   "1.5 0.0 -0.3 0.0 0.0 0.0 1.0 0.0\n\n"
   "0:1.5 2:-0.3 6:1.0 ... \n\n"
   "[format=libsvm] (label (colum_id:value)+\\n)+ (label will be ignored)\n"
   "+1 0:1.5 2:-0.3 6:1.0 ... \n\n"
   "Use DenseMatrix for dense, and use SparseMatrix for sparse and libsvm\n\n"
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
  p.add<string>("format", 'f', "format type (dense|sparse|libsvm) See example. ", false, "dense");
  p.add<string>("method", 'm', "method (SVD|PCA)", false, "SVD");
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

  string method = p.get<string>("method");
  try {
    if (method == "SVD"){
      REDSVD::SVDfromFile(p.get<string>("input"), 
			  p.get<string>("output"), 
			  p.get<string>("format"), 
			  p.get<int>("rank"));
    } else if (method == "PCA"){

      REDSVD::PCAfromFile(p.get<string>("input"), 
			  p.get<string>("output"), 
			  p.get<string>("format"), 
			  p.get<int>("rank"));
    } else {
      cerr << "unknown method:" << method << endl;
      return -1;
    }
  } catch (string& error){
    cerr << "Error: " << error << endl;
    return -1;
  }
  return 0;
}
