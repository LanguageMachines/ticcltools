#include <iostream>
#include <fstream>
#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace TiCC;

void usage( const string& name ){
  cerr << name << " --data=datafile [-n size] [FILES]" << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hn:", "data:" );
  try {
    opts.init(argc,argv);
  }
  catch( OptionError& e ){
    cerr << e.what() << endl;
    usage( opts.prog_name() );
    exit( EXIT_FAILURE );
  }
  if ( opts.extract( 'h' ) ){
    usage( opts.prog_name() );
    exit( EXIT_SUCCESS );
  }
  string dataFile;
  if ( !opts.extract( "data", dataFile ) ){
    cerr << "missing '--data' option" << endl;
    exit( EXIT_FAILURE );
  }
  int NN = 20;
  string value;
  if ( opts.extract( 'n', value ) ){
    NN = stringTo<int>(value);
  }
  auto fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing input file(s)" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.empty() ) {
    cerr << "unsupported options: " << opts.toString() << endl;
    exit( EXIT_FAILURE );
  }
  wordvec_tester WV;
  if ( !WV.fill( dataFile ) ){
    cerr << "fill failed from " << dataFile << endl;
    exit(EXIT_FAILURE);
  }
  else
    cerr << "filled with " << WV.size() << " vectors" << endl;
  for ( auto const& name : fileNames ){
    ifstream is( name );
    if ( !is ){
      cerr << "failed to read: " << name << endl;
      continue;
    }
    string outname = name + ".out";
    ofstream os( outname );
    if ( !os ){
      cerr << "failed to open: " << outname << endl;
      continue;
    }
    string line;
    while ( getline( is, line ) ){
      os << "NEIGHBORS of '" << line << "':" << endl;
      vector<word_dist> res;
      if ( WV.lookup( line, NN, res ) ){
	for ( auto const& i : res ){
	  os << "\t" << i.w << "\t" << i.d << endl;
	}
      }
      else {
	os << "\tNone" << endl;
      }
    }
    cerr << "results in: " << outname << endl;
  }
  return EXIT_SUCCESS;
}
