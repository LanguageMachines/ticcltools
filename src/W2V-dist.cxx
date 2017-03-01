#include <map>
#include <iostream>
#include <fstream>
#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace TiCC;

void usage( const string& name ){
  cerr << name << " --vectors=vectorfile [FILES]" << endl;
}

bool fill( const string& freqsFile, map<string,size_t>& freqs ){
  ifstream is( freqsFile );
  if ( !is.good() ){
    return false;
  }
  string line;
  while ( getline( is, line ) ){
    if ( line.empty() )
      continue;
    vector<string> parts;
    size_t num = split(line, parts);
    if ( num != 2 )
      continue;
    size_t freq;
    if ( !stringTo( parts[1], freq ) ){
      cerr << "error in line: '" << line << "' second part isn't an integer?"
	   << endl;
      continue;
    }
    freqs[parts[0]] = freq;
  }
  return true;
}

size_t lookup( const string& s, const map<string,size_t> frqs ){
  const auto it = frqs.find( s );
  if ( it != frqs.end() )
    return it->second;
  else
    return 0;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "h", "vectors:,freqs:" );
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
  string vectorsFile;
  if ( !opts.extract( "vector", vectorsFile ) ){
    cerr << "missing '--vectors' option" << endl;
    exit( EXIT_FAILURE );
  }
  string freqsFile;
  opts.extract( "freqs", freqsFile );
  auto fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing input file(s)" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.empty() ) {
    cerr << "unsupported options: " << opts.toString() << endl;
    exit( EXIT_FAILURE );
  }
  map<string,size_t> freqs;
  if ( !freqsFile.empty() ){
    if ( fill( freqsFile, freqs ) ){
      cerr << "filled a frequency hash with " << freqs.size()
	   << " entries from " << freqsFile << endl;
    }
    else {
      cerr << "problem reading freqs file: " << freqsFile << endl;
      exit( EXIT_FAILURE );
    }
  }
  wordvec_tester WV;
  if ( !WV.fill( vectorsFile ) ){
    cerr << "fill failed from " << vectorsFile << endl;
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
	if ( !freqs.empty() ){
	  size_t f1 = lookup( parts[0], freqs );
	  size_t f2 = lookup( parts[1], freqs );
	  os << parts[0] << "\t" << f1 << "\t"
	     << parts[1] << "\t" << f2 << "\t" << cosine << endl;
	}
	else {
	  os << line << "\t" << cosine << endl;
	}
      }
      catch( ... ){
	os << line << "\tUNKNOWN word(s)" << endl;
      }
    }
    cerr << "results in: " << outname << endl;
  }
  return EXIT_SUCCESS;
}
