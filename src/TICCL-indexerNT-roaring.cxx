/*
  Copyright (c) 2006 - 2021
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
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <climits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "config.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"
#include "roaring/roaring64map.hh"
#include "config.h"

using namespace std;
using namespace icu;

//#define TRANSPOSE_TEST 1

void usage( const string& name ){
  cerr << name << endl;
  cerr << "options: " << endl;
  cerr << "\t--hash=<anahash>\tname of the anagram hashfile. (produced by TICCL-anahash)" << endl;
  cerr << "\t--charconf=<charconf>\tname of the character confusion file. (produced by TICCL-lexstat)" << endl;
  cerr << "\t--foci=<focifile>\tname of the file produced by the --artifrq parameter of TICCL-anahash" << endl;
  cerr << "\t-o <outputfile>\tname for the outputfile. " << endl;
  cerr << "\t--low=<low>\t skip entries from the anagram file shorter than "
       << endl;
  cerr << "\t\t'low' characters. (default = 5)" << endl;
  cerr << "\t--high=<high>\t skip entries from the anagram file longer than "
       << endl;
  cerr << "\t\t'high' characters. (default=35)" << endl;
  cerr << "\t-t <threads>\t\trun on 'threads' threads." << endl;
  cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  cerr << "\t\t\t reasonable value. (OMP_NUM_TREADS - 2)" << endl;
  cerr << "\t-V show version " << endl;
  cerr << "\t-h this message " << endl;
}

struct experiment {
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
    e.start = hashes.begin();
    e.finish = hashes.end();
    exps.push_back( e );
    return 1;
  }
  set<bitType>::const_iterator s = hashes.begin();
  for ( size_t i=0; i < threads; ++i ){
    experiment e;
    e.start = s;
    for ( size_t j=0; j < partsize; ++j )
      ++s;
    e.finish = s;
    exps.push_back( e );
  }
  if ( s != hashes.end() ){
    exps[exps.size()-1].finish = hashes.end();
  }
  return threads;
}

void handle_exp( const experiment& exp,
		 size_t& count,
		 const set<bitType>& hashSet,
		 const set<bitType>& confSet,
		 map<bitType,Roaring64Map>& r_result ){
  bitType max = *confSet.rbegin();
  auto it1 = exp.start;
  while ( it1 != exp.finish ){
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
    set<bitType>::const_iterator it3 = hashSet.find( *it1 );
    if ( it3 != hashSet.end() ){
      set<bitType>::const_reverse_iterator it2( it3 );
      while ( it2 != hashSet.rend() ){
	bitType diff = *it1 - *it2;
	if ( diff > max )
	  break;
	set<bitType>::const_iterator sit = confSet.find( diff );
	if ( sit != confSet.end() ){
#pragma omp critical
	  {
#ifdef TRANSPOSE_TEST
	    r_result[*it2].add( diff );
#else
	    r_result[diff].add( *it2 );
#endif
	  }
	}
	++it2;
      }
      // it3 is already set at hashSet.find( *it1 );
      ++it3;
      while ( it3 != hashSet.end() ){
	bitType diff = *it3 - *it1;
	if ( diff > max )
	  break;
	set<bitType>::const_iterator sit = confSet.find( diff );
	if ( sit != confSet.end() ){
#pragma omp critical
	  {
#ifdef TRANSPOSE_TEST
	    r_result[*it1].add( diff );
#else
	    r_result[diff].add( *it1 );
#endif
	  }
	}
	++it3;
      }
    }
    ++it1;
  }
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "charconf:,hash:,low:,high:,foci:,help,version,threads:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('h' ) || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V' ) || opts.extract("version") ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( argc < 4	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  bool verbose = opts.extract( 'v' );
  string anahashFile;
  string confFile;
  string outFile;
  string fociFile;
  int lowValue = 5;
  int highValue = 35;
  int num_threads = 1;
  if ( !opts.extract( "hash", anahashFile ) ){
    cerr << "missing --hash option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "charconf", confFile ) ){
    cerr << "missing --charconf option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "foci", fociFile ) ){
    cerr << "missing --foci option" << endl;
    exit( EXIT_FAILURE );
  }
  opts.extract( 'o', outFile );
  string value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
#ifdef HAVE_OPENMP
  if ( TiCC::lowercase(value) == "max" ){
    num_threads = omp_get_max_threads() - 2;
  }
  else {
    if ( !TiCC::stringTo(value,num_threads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
#else
  if ( value != "1" ){
    cerr << "unable to set number of threads!.\nNo OpenMP support available!"
	 <<endl;
    exit(EXIT_FAILURE);
  }
#endif
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

  ifstream cwav( anahashFile );
  if ( !cwav ){
    cerr << "problem opening corpus word anagram hash file: "
	 << anahashFile << endl;
    exit(1);
  }
  if ( outFile.empty() ){
    outFile = anahashFile;
    string::size_type pos = outFile.rfind(".");
    if ( pos != string::npos ){
      outFile = outFile.substr(0,pos);
    }
    outFile += ".indexNT.R";
  }
  else if ( !TiCC::match_back( outFile, ".indexNT.R" ) ){
    outFile += ".indexNT.R";
  }
#ifdef TRANSPOSE_TEST
  outFile += ".T";
#endif

  ofstream of( outFile );
  if ( !of ){
    cerr << "problem opening output file: " << outFile << endl;
    exit(1);
  }

  ifstream conf( confFile );
  if ( !conf ){
    cerr << "problem opening character confusion anagram file: "
	 << confFile << endl;
    exit(1);
  }

  ifstream foc( fociFile );
  if ( !foc ){
    cerr << "problem opening foci file: " << fociFile << endl;
    exit(1);
  }

  cout << "reading corpus word anagram hash values" << endl;
  size_t skipped = 0;
  set<bitType> hashSet;
  string line;
  while ( getline( cwav, line ) ){
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "~" ) > 1 ){
      bitType bit = TiCC::stringTo<bitType>( parts[0] );
      vector<string> parts2;
      if ( TiCC::split_at( parts[1], parts2, "#" ) > 0 ){
	UnicodeString firstItem = TiCC::UnicodeFromUTF8( parts2[0] );
	if ( firstItem.length() >= lowValue &&
	     firstItem.length() <= highValue ){
	  hashSet.insert( bit );
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
  cout << "read " << hashSet.size() << " corpus word anagram values" << endl;
  cout << "skipped " << skipped << " out-of-band corpus word values" << endl;

  set<bitType> focSet;
  while ( foc ){
    bitType bit;
    foc >> bit;
    foc.ignore( INT_MAX, '\n' );
    focSet.insert( bit );
  }
  cout << "read " << focSet.size() << " foci values" << endl;

  set<bitType> confSet;
  while ( getline( conf, line ) ){
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "#" ) > 0 ){
      bitType bit = TiCC::stringTo<bitType>( parts[0] );
      confSet.insert(bit);
    }
    else {
      cerr << "problems with line " << line << endl;
      cerr << "bail out " << endl;
      exit(1);
    }
  }
  cout << "read " << confSet.size()
       << " character confusion anagram values" << endl;

  vector<experiment> experiments;
  size_t expsize = init( experiments, focSet, num_threads );

  cout << "created " << expsize << " separate experiments" << endl;

#ifdef HAVE_OPENMP
  omp_set_num_threads( expsize );
#endif

  size_t count = 0;
  map<bitType,Roaring64Map> r_result;
#pragma omp parallel for shared(experiments, count , r_result )
  for ( size_t i=0; i < expsize; ++i ){
    handle_exp( experiments[i], count, hashSet, confSet, r_result );
  }

  for ( auto const& rit : r_result ){
    of << rit.first;
    of << "#";
    uint64_t expectedsize = rit.second.getSizeInBytes();
    char *serializedbytes = new char [expectedsize];
    uint64_t uit = rit.second.write(serializedbytes);
    of << uit << " " << string(serializedbytes,uit) << endl;
    delete [] serializedbytes;
  }
}
