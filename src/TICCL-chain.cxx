/*
  Copyright (c) 2006 - 2018
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
#include <functional>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>
#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/Unicode.h"
#include "ticcl/unicode.h"
#include "ticcl/word2vec.h"

using namespace std;
typedef signed long int bitType;

const int RANK_COUNT=13;

bool verbose = false;

void usage( const string& name ){
  cerr << "usage: " << name << endl;
  exit( EXIT_FAILURE );
}

unsigned int ld( const string& in1, const string& in2 ){
  UnicodeString s1 = TiCC::UnicodeFromUTF8(in1);
  UnicodeString s2 = TiCC::UnicodeFromUTF8(in2);
  const size_t len1 = s1.length(), len2 = s2.length();
  vector<unsigned int> col(len2+1), prevCol(len2+1);
  for ( unsigned int i = 0; i < prevCol.size(); ++i ){
    prevCol[i] = i;
  }
  for ( unsigned int i = 0; i < len1; ++i ) {
    col[0] = i+1;
    for ( unsigned int j = 0; j < len2; ++j )
      col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]),
		      prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
    col.swap(prevCol);
  }
  unsigned int result = prevCol[len2];
  return result;
}

void calc_chain( ostream& os,
		 string root,
		 size_t root_frq,
		 string s2,
		 map<string, set<string>>& table,
		 map<string, size_t>& var_freq,
		 set<string>& done ){
  if ( verbose ){
    using TiCC::operator<<;
    cerr << "doorzoek met:" << s2 << " " << table[s2] << endl;
  }
  for ( const auto& it : table[s2] ){
    if ( table[it].empty() ){
      if ( done.find( it ) != done.end() ){
	continue;
      }
      os << it << "#" << var_freq[it] << "#" << root << "#"
	 << root_frq << "#" << ld( it, s2 ) << "#C" << endl;
      done.insert( it );
    }
    else {
      calc_chain( os, root, root_frq, it, table, var_freq, done );
    }
  }
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( argc < 2	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  if ( opts.extract('h' ) ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  bool verbose = opts.extract( 'v' );
  int numThreads=1;
  string outFile;
  opts.extract( 'o', outFile );
  string value;
  if ( opts.extract( 't', value ) ){
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  vector<string> fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing an inputfile" << endl;
    exit(EXIT_FAILURE);
  }
  if ( fileNames.size() > 1 ){
    cerr << "only one inputfile may be provided." << endl;
    exit(EXIT_FAILURE);
  }
  string inFile = fileNames[0];
  if ( !TiCC::match_back( inFile, ".ranked" ) ){
    cerr << "inputfile must have extension .ranked" << endl;
    exit(EXIT_FAILURE);
  }
  if ( !outFile.empty() ){
    if ( !TiCC::match_back( outFile, ".chained" ) )
      outFile += ".chained";
  }
  else {
    outFile = inFile + ".chained";
  }
  if ( outFile == inFile ){
    cerr << "same filename for input and output!" << endl;
    exit(EXIT_FAILURE);
  }

  ifstream input( inFile );
  if ( !input ){
    cerr << "problem opening input file: " << inFile << endl;
    exit(1);
  }

  std::multimap< size_t, std::string, std::greater<size_t> > desc_freq;
  std::map< std::string, size_t > var_freq;
  set<string> done;
  map<string, set<string> > table;
  ofstream os( outFile );
  string line;
  while( getline( input, line ) ){
    vector<string> parts = TiCC::split_at( line, "#" );
    if ( parts.size() != 6 ){
      cerr << "invalid line: '" << line << "'" << endl;
    }
    else {
      string variant1 = parts[0];
      size_t freq1 = TiCC::stringTo<size_t>(parts[1]);
      var_freq[variant1] = freq1;
      string variant2 = parts[2];
      size_t freq2 = TiCC::stringTo<size_t>(parts[3]);
      // size_t ld    = TiCC::stringTo<int>(parts[4]);
      // double ld_rank = TiCC::stringTo<double>(parts[5]);
      if ( done.find( variant2 ) == done.end() ){
	desc_freq.insert( make_pair(freq2,variant2) );
	done.insert( variant2 );
      }
      table[variant2].insert( variant1 );
    }
  }
  ofstream db( outFile + ".debug" );
  if ( verbose ){
    using TiCC::operator<<;
    for ( const auto& val : desc_freq ){
      db << val.first << " " << val.second << " " << table[val.second] << endl;
    }
    cout << "debug results in " << outFile + ".debug" << endl;
  }
  for ( const auto& val : desc_freq ){
    done.clear();
    calc_chain( os, val.second, val.first, val.second, table, var_freq, done );
  }

  cout << "results in " << outFile << endl;
}
