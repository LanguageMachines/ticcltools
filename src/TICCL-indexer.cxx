#include <unistd.h>
#include <set>
#include <limits>
#include <algorithm>
#include <vector>
#include <climits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcl/unicode.h"

#include "config.h"

using namespace std;
typedef signed long int bitType;

void usage( const string& name ){
  cerr << name << endl;
  cerr << "options: " << endl;
  cerr << "\t--hash=<anahash>\tname of the anagram hashfile. (produced by TICCL-anahash)" << endl;
  cerr << "\t--charconf=<charconf>\tname of the character confusion file. (produced by TICCL-lexstat)" << endl;
  cerr << "\t-o <outputfile>\tname for the outputfile. " << endl;
  cerr << "\t--low=<low>\t skip entries from the anagram file shorter than "
       << endl;
  cerr << "\t\t'low' characters. (default = 5)" << endl;
  cerr << "\t--high=<high>\t skip entries from the anagram file longer than "
       << endl;
  cerr << "\t\t'high' characters. (default=35)" << endl;
  cerr << "\t-V show version " << endl;
  cerr << "\t-v verbosity " << endl;
  cerr << "\t-h this message " << endl;
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "charconf:,hash:,low:,high:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('h' ) ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( argc < 3	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  bool verbose = opts.extract( 'v' );
  string anahashFile;
  string confFile;
  string outFile;
  int lowValue = 5;
  int highValue = 35;
  opts.extract( "hash", anahashFile );
  opts.extract( "charconf", confFile );
  opts.extract( 'o', outFile );
  string value;
  if ( opts.extract("low", value ) ){
    if ( !TiCC::stringTo(value,lowValue) ) {
      cerr << "illegal value for --low (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract("high", value ) ){
    if ( !TiCC::stringTo(value,highValue) ) {
      cerr << "illegal value for --high (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  ifstream ana( anahashFile.c_str() );
  if ( !ana ){
    cerr << "problem opening anagram hashfile: " << anahashFile << endl;
    exit(1);
  }
  ifstream conf( confFile.c_str() );
  if ( !conf ){
    cerr << "problem opening charconfusion file: " << confFile << endl;
    exit(1);
  }
  if ( outFile.empty() ){
    outFile = anahashFile;
    string::size_type pos = outFile.rfind(".");
    if ( pos != string::npos ){
      outFile = outFile.substr(0,pos);
    }
    outFile += ".index";
  }
  else if ( !TiCC::match_back( outFile, ".index" ) ){
    outFile += ".index";
  }

  ofstream of( outFile.c_str() );
  if ( !of ){
    cerr << "problem opening outputfile: " << outFile << endl;
    exit(1);
  }
  cout << "reading anagram hash values" << endl;
  size_t skipped = 0;
  set<bitType> anaSet;
  string line;
  while ( getline( ana, line ) ){
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "~" ) > 1 ){
      bitType bit = TiCC::stringTo<bitType>( parts[0] );
      vector<string> parts2;
      if ( TiCC::split_at( parts[1], parts2, "#" ) > 0 ){
	UnicodeString firstItem = UTF8ToUnicode( parts2[0] );
	if ( firstItem.length() >= lowValue &&
	     firstItem.length() <= highValue ){
	  anaSet.insert( bit );
	}
	else {
	  if ( verbose ){
	    cerr << "skip " << parts2[0] << endl;
	  }
	  ++skipped;
	}
      }
    }
  }
  cout << "read " << anaSet.size() << " anagram values" << endl;
  cout << "skipped " << skipped << " out-of-band anagram values" << endl;

  cout << "reading character confusion values" << endl;
  set<bitType> confSet;
  size_t count = 0;
  while ( getline( conf, line ) ){
    vector<string> parts;
    if ( ++count % 1000 == 0 ){
      cout << ".";
      cout.flush();
      if ( count % 50000 == 0 ){
	cout << endl << count << endl;;
      }
    }
    if ( TiCC::split_at( line, parts, "#" ) > 0 ){
      bitType bit = TiCC::stringTo<bitType>( parts[0] );
      confSet.insert( bit );
    }
    else {
      cerr << "problems with line " << line << endl;
      cerr << "bail out " << endl;
      exit(1);
    }
  }
  cout << "read " << confSet.size() << " confusion values" << endl;

  bitType vorige = 0;
  cout << "processing all confusion values" << endl;
  bitType totalShift = 0;
  count = 0;
  for ( const auto& bit : confSet ){
    bitType diff = bit - vorige;
    totalShift += diff;
    if ( ++count % 100 == 0 ){
      cout << ".";
      cout.flush();
      if ( count % 5000 == 0 ){
	cout << endl << count << endl;;
      }
    }
    auto it1 = anaSet.begin();
    auto it2 = anaSet.begin();
    bool doKomma = false;
    bool first = true;
    while ( it1 != anaSet.end() && it2 != anaSet.end() ){
      bitType v1 = *it1;
      bitType v2 = *it2 - totalShift;
      if ( v1 == v2 ){
	if ( doKomma ){
	  of << ",";
	}
	else {
	  doKomma = true;
	}
	if ( first ){
	  // avoid empty entries
	  of << bit << "#";
	  first = false;
	}
	of << v1;
	++it1;
	++it2;
      }
      else if ( v1 < v2 ){
	++it1;
      }
      else
	++it2;
    }
    if ( !first )
      of << endl;
    vorige = bit;
  }

}
