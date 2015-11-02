#include <iostream>
#include <fstream>
#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace TiCC;

void usage( const string& name ){
  cerr << name << " --data=datafile [FILES]" << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "h", "data:" );
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
    int err_cnt = 5;
    while ( getline( is, line ) ){
      vector<string> parts;
      if ( TiCC::split_at_first_of( line, parts, "\t#" ) != 2 ){
	cerr << "invalid line (expected 2 tab or # seperated words/sentences)."
	     << endl;
	if ( --err_cnt > 0 )
	  continue;
	else {
	  cerr << "skiping file " << name << " (too many errors)" << endl;
	  break;
	}
      }
      double cosine = 0.0;
      try {
	cosine = WV.distance( parts[0], parts[1] );
	os << line << "\t" << cosine << endl;
      }
      catch( ... ){
	os << line << "\tUNKNOWN word(s)" << endl;
      }
    }
    cerr << "results in: " << outname << endl;
  }
  return EXIT_SUCCESS;
}
