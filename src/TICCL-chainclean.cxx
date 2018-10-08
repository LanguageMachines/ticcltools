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
using TiCC::operator<<;

void usage( const string& name ){
  cerr << "usage: " << name << "[options] chainfile " << endl;
  cerr << "\t\t The chainfiles is an outputfile from TICCL-chain." << endl;
  cerr << "-t--lexicon A validated lexicon." << endl;
  cerr << "-t--artifrq The artifreq. Default 100000000 ." << endl;
  cerr << "\t-o <outputfile> name of the outputfile." << endl;
  cerr << "\t-h or --help this message." << endl;
  cerr << "\t-v be verbose, repeat to be more verbose. " << endl;
  cerr << "\t-V or --version show version. " << endl;
  exit( EXIT_FAILURE );
}

const string SEPARATOR = "_";

class record {
public:
  string variant;
  vector<string> v_parts;
  string v_freq;
  string cc;
  string cc_freq;
  string ld;
};

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "lexicon:,artifreq:" );
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
  int verbosity = 0;
  while( opts.extract( 'v' ) ){
    ++verbosity;
  }
  unsigned int artifreq = 100000000;
  string value;
  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
  string lex_name;
  opts.extract( "lexicon", lex_name );
  if ( lex_name.empty() ){
    cerr << "missing --lexcion options"<< endl;
    exit(EXIT_FAILURE);
  }
  string out_name;
  opts.extract( 'o', out_name );
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
  string in_name = fileNames[0];
  if ( !out_name.empty() ){
    if ( out_name == in_name ){
      cerr << "same filename for input and output!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  ifstream input( in_name );
  if ( !input ){
    cerr << "problem opening input file: " << in_name << endl;
    exit(1);
  }
  set<string> valid_words;
  ifstream lexicon( lex_name );
  string line;
  while ( getline( lexicon, line ) ){
    if ( line.size() == 0 || line[0] == '#' )
      continue;
    vector<string> vec = TiCC::split( line );
    if ( vec.size() < 2 ){
      cerr << progname << ": invalid line '" << line << "' in " << lex_name << endl;
      exit( EXIT_FAILURE );
    }
    unsigned int freq = 0;
    if ( !TiCC::stringTo(vec[1], freq) ) {
      cerr << progname << ": invalid frequency in '" << line << "' in " << lex_name << endl;
      exit( EXIT_FAILURE );
    }
    if ( freq >= artifreq ){
      valid_words.insert( vec[0] );
    }
  }
  cout << "read " << valid_words.size() << " validated words from "
       << lex_name << endl;
  cout << "start reading chained results" << endl;
  list<record> records;
  while ( getline( input, line ) ){
    vector<string> vec = TiCC::split_at( line, "#" );
    if ( vec.size() != 6 ){
      cerr << progname << ": chained file should have 6 items per line: '" << line << "' in " << lex_name << endl;
      cerr << "\t found " << vec.size() << endl;
      exit( EXIT_FAILURE );
    }
    record rec;
    rec.variant = vec[0];
    rec.v_freq = vec[1];
    rec.cc = vec[2];
    rec.cc_freq = vec[3];
    rec.ld = vec[4];
    records.push_back( rec );
  }
  cout << "start processing " << records.size() << " chained results" << endl;
  map<string,int> parts_freq;
  for ( auto& rec : records ){
    rec.v_parts = TiCC::split_at( rec.variant, SEPARATOR );
    if ( rec.v_parts.size() == 1 ){
      continue;
    }
    for ( const auto& p : rec.v_parts ){
      if ( valid_words.find( p ) == valid_words.end() ){
	++parts_freq[p];
      }
    }
  }
  cout << "found " << parts_freq.size() << " unknown parts" << endl;
  list<record> output_records;
  for ( const auto& part : parts_freq ) {
    map<string,int> cc_freqs;
    auto it = records.begin();
    while ( it != records.end() ){
      if ( it->v_parts.size() != 1 ){
	bool match = false;
	for ( const auto& p : it->v_parts ){
	  if ( p == part.first ){
	    match = true;
	    break;
	  }
	}
	if ( match ){
	  ++cc_freqs[it->cc];
	}
      }
      ++it;
    }
    cerr << "found " << cc_freqs.size() << " CC's for: " << part.first << endl;
    cerr << cc_freqs << endl;
    for ( const auto& cc : cc_freqs ){
      auto it = records.begin();
      while ( it != records.end() ){
	if ( cc.first == it->cc ){
	  bool match = false;
	  for ( const auto& p : it->v_parts ){
	    if ( p == part.first ){
	      match = true;
	      break;
	    }
	  }
	  if ( match ){
	    output_records.push_back( *it );
	  }
	}
	++it;
      }
    }
  }
  ofstream os( out_name );
  for ( const auto& it : output_records ){
    os << it.variant << "#" << it.v_freq << "#" << it.cc << "#"
       << it.cc_freq << "#" << it.ld << "#C" << endl;
  }
  cout << "results in " << out_name << endl;
}
