/*
  Copyright (c) 2006 - 2017
  CLST  - Radboud University
  ILK   - Tilburg University

  This file is part of ticcltools

  ticcltools is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  ticcltools is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

  For questions and suggestions, see:
      https://github.com/LanguageMachines/ticcltools/issues
  or send mail to:
      lamasoftware (at ) science.ru.nl

*/
#include <unistd.h>
#include <set>
#include <map>
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
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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
  cerr << "\t--foci=<focifile>\tname of the file produced by the --artifrq parameter of TICCL-anahash." << endl;
  cerr << "\t\tThis file is used to limit the searchspace" << endl;
  cerr << "\t-t <threads>\t\trun on 'threads' threads." << endl;
  cerr << "\t-V or --version show version " << endl;
  cerr << "\t-v verbosity " << endl;
  cerr << "\t-h or --help this message " << endl;
}

struct experiment {
  bitType vorige;
  bitType totalShift;
  set<bitType>::const_iterator start;
  set<bitType>::const_iterator finish;
};

size_t init( vector<experiment>& exps,
	     const set<bitType>& hashes,
	     size_t threads ){
  exps.clear();
  size_t partsize = hashes.size() / threads;
  if ( partsize < 1 ){
    experiment e;
    e.vorige = 0;
    e.totalShift = 0;
    e.start = hashes.begin();
    e.finish = hashes.end();
    exps.push_back( e );
    return 1;
  }
  set<bitType>::const_iterator s = hashes.begin();
  for ( size_t i=0; i < threads; ++i ){
    experiment e;
    e.start = s;
    for ( size_t j=0; j < partsize && s != hashes.end(); ++j ){
      ++s;
    }
    e.finish = s;
    exps.push_back( e );
  }
  if ( s != hashes.end() ){
    exps[exps.size()-1].finish = hashes.end();
  }
  return threads;
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "charconf:,hash:,low:,high:,help,version,foci:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('h') || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V') || opts.extract("version") ){
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
  string fociFile;
  string outFile;
  int lowValue = 5;
  int highValue = 35;
  int threads = 1;
  opts.extract( "hash", anahashFile );
  opts.extract( "charconf", confFile );
  opts.extract( "foci", fociFile );
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
  if ( opts.extract('t', value ) ){
#ifdef HAVE_OPENMP
    if ( !TiCC::stringTo(value,threads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
#else
    cerr << "You don't have OpenMP support. Setting -t is useless!" << endl;
    exit( EXIT_FAILURE );
#endif
  }

  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  ifstream ana( anahashFile );
  if ( !ana ){
    cerr << "problem opening anagram hashfile: " << anahashFile << endl;
    exit(1);
  }
  ifstream conf( confFile );
  if ( !conf ){
    cerr << "problem opening charconfusion file: " << confFile << endl;
    exit(1);
  }

  set<bitType> focSet;
  if ( !fociFile.empty() ){
    ifstream foc( fociFile );
    if ( !foc ){
      cerr << "problem opening foci file: " << fociFile << endl;
      exit(1);
    }
    while ( foc ){
      bitType bit;
      foc >> bit;
      foc.ignore( INT_MAX, '\n' );
      focSet.insert( bit );
    }
    cout << "read " << focSet.size() << " foci values" << endl;
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

  ofstream of( outFile );
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

  map<bitType,set<bitType> > result;
  vector<experiment> experiments;
  size_t expsize = init( experiments, confSet, threads );
#ifdef HAVE_OPENMP
  omp_set_num_threads( expsize );
#endif

  cout << "processing all confusion values" << endl;
  count = 0;
#pragma omp parallel for shared( experiments, anaSet, focSet )
  for ( size_t i=0; i < expsize; ++i ){
    set<bitType>::const_iterator bit = experiments[i].start;
    while ( bit != experiments[i].finish ){
      bitType confusion = *bit;
      bitType diff = confusion - experiments[i].vorige;
      experiments[i].totalShift += diff;
#pragma omp critical
      {
	if ( ++count % 100 == 0 ){
	  cout << ".";
	  cout.flush();
	  if ( count % 5000 == 0 ){
	    cout << endl << count << endl;;
	  }
	}
      }
      auto it1 = anaSet.begin();
      auto it2 = anaSet.begin();
      while ( it1 != anaSet.end() && it2 != anaSet.end() ){
	bitType v1 = *it1;
	bitType v2 = *it2;
	bool foc = true;
	if ( !focSet.empty() ){
	  // do we have to some focus?
	  foc = !( focSet.find( v1 ) == focSet.end()
		   && focSet.find( v2 ) == focSet.end() );
	}
	v2 -= experiments[i].totalShift;
	if ( v1 == v2 ){
	  if ( foc ){
#pragma omp critical
	    {
	      cerr << "YES " << confusion << endl;
	      result[confusion].insert(v1);
	    }
	  }
	  ++it1;
	  ++it2;
	}
	else if ( v1 < v2 ){
	  ++it1;
	}
	else {
	  ++it2;
	}
      }
      experiments[i].vorige = confusion;
      ++bit;
    }
  }

  for ( auto const& rit : result ){
    of << rit.first << "#";
    set<bitType>::const_iterator it = rit.second.begin();
    while ( it != rit.second.end() ){
      of << *it;
      ++it;
      if ( it != rit.second.end() ){
	of << ",";
      }
    }
    of << endl;
  }

}
