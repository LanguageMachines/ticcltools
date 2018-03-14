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

unsigned int ldCompare( const UnicodeString& s1, const UnicodeString& s2 ){
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

unsigned int ld( const string& in1, const string& in2, bool caseless ){
  UnicodeString s1 = TiCC::UnicodeFromUTF8(in1);
  UnicodeString s2 = TiCC::UnicodeFromUTF8(in2);
  if ( caseless ){
    s1.toLower();
    s2.toLower();
  }
  return ldCompare( s1, s2 );
}

class chain_class {
public:
  chain_class(): chain_class( 0,false ){};
  chain_class( int v, bool c ): verbosity(v), caseless(c){};
  bool fill( const string& );
  void debug_info( const string& );
  void output( const string& );
  void calc_chain( ostream&, const string&, size_t,
		   const set<string>&, set<string>& ) const;
private:
  multimap< size_t, string, std::greater<size_t> > desc_freq;
  map<string, set<string>> table;
  map< string, size_t > var_freq;
  set<string> best_trans;
  set<string> done;
  int verbosity;
  bool caseless;

};

bool chain_class::fill( const string& line ){
  vector<string> parts = TiCC::split_at( line, "#" );
  if ( parts.size() != 6 ){
    return false;
  }
  else {
    string variant1 = parts[0]; // a possibly correctable word
    size_t freq1 = TiCC::stringTo<size_t>(parts[1]);
    var_freq[variant1] = freq1;
    string variant2 = parts[2]; // a Correction Candidate
    size_t freq2 = TiCC::stringTo<size_t>(parts[3]);
    var_freq[variant2] = freq2;
    if ( done.find( variant2 ) == done.end() ){
      // save CC's frequency only once
      desc_freq.insert( make_pair(freq2,variant2) );
      done.insert( variant2 );
    }
    // is this word already stored?
    if ( best_trans.find( variant1 ) == best_trans.end() ){
      // NO, so add is to the list of the CC
      table[variant2].insert( variant1 );
      best_trans.insert( variant1 );
    }
    return true;
  }
}

void chain_class::debug_info( const string& name ){
  string out_file = name + ".debug";
  ofstream db( out_file );
  using TiCC::operator<<;
  for ( const auto& val : desc_freq ){
    if ( table.find( val.second ) != table.end() ){
      db << val.first << " " << val.second
	 << " " << table[val.second] << endl;
    }
  }
  cout << "debug info stored in " << out_file << endl;
}

void chain_class::output( const string& out_file ){
  ofstream os( out_file );
  set<string> done;
#pragma omp parallel for shared(done)
  for ( size_t i=0; i < desc_freq.size(); ++i ){
    auto desc_it = desc_freq.begin();
    advance( desc_it, i );
    auto const it = table.find( desc_it->second );
    if ( it != table.end() ) {
      calc_chain( os, desc_it->second, desc_it->first, it->second, done );
    }
  }
}

void chain_class::calc_chain( ostream& os,
			      const string& root,
			      size_t root_frq,
			      const set<string>& candidates,
			      set<string>& done ) const {
  if ( verbosity > 2 ){
#pragma omp critical
    {
      using TiCC::operator<<;
      cerr << "doorzoek " << candidates << endl;
    }
  }
  for ( const auto& candidate: candidates ){
    if ( verbosity > 2 ){
#pragma omp critical
      {
	cerr << "kandidaat :" << candidate << endl;
      }
    }
    map<string, set<string>>::const_iterator table_it;
#pragma omp critical
    {
      table_it = table.find(candidate);
    }
    // do they have CC themselves?
    if ( table_it != table.end() ){
      if ( verbosity > 2 ){
#pragma omp critical
	{
	  cerr << " recurse for: " << candidate << endl;
	}
      }
      calc_chain( os, root, root_frq, table_it->second, done );
    }
    else {
      if ( verbosity > 2 ){
#pragma omp critical
	{
	  cerr << " no more entries for: " << candidate << endl;
	}
      }
    }
#pragma omp critical
    {
      if ( done.find( candidate ) == done.end() ){
	// Didn't we yet output a translation?
	done.insert( candidate );
	os << candidate << "#" << var_freq.at(candidate) << "#" << root << "#"
	   << root_frq << "#" << ld( root, candidate, caseless ) << "#C" << endl;
	if ( verbosity > 2 ){
	  cerr << "      output: " << candidate << " ==> " << root << endl;
	}
      }
    }
  }
}

void usage( const string& name ){
  cerr << "usage: " << name << endl;
  // cerr << "\t-t <threads>\n\t--threads <threads> Number of threads to run on." << endl;
  // cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  // cerr << "\t\t\t reasonable value. (which can be set with OMP_NUM_TREADS environment variable.)" << endl;
  cerr << "\t--caseless Calculate the Levensthein (or edit) distance ignoring case." << endl;
  cerr << "\t-o <outputfile> name of the outputfile." << endl;
  cerr << "\t-h or --help this message." << endl;
  cerr << "\t-v be verbose, repeat to be more verbose. " << endl;
  cerr << "\t-V or --version show version. " << endl;
  exit( EXIT_FAILURE );
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "caseless,threads:" );
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
  bool caseless = opts.extract( "caseless" );
  int numThreads=1;
  string out_file;
  opts.extract( 'o', out_file );
  string value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
  // if ( value != "1" ){
  //   cerr << "-t or --threads options not supported!" << endl;
  //   exit( EXIT_FAILURE );
  // }
#ifdef HAVE_OPENMP
  if ( TiCC::lowercase(value) == "max" ){
    numThreads = omp_get_max_threads();
    omp_set_num_threads( numThreads );
    cout << "runing on " << numThreads << " threads." << endl;
  }
  else {
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    omp_set_num_threads( numThreads );
    cout << "runing on " << numThreads << " threads." << endl;
  }
#else
  if ( value != "1" ){
    cerr << "unable to set number of threads!.\nNo OpenMP support available!"
	 <<endl;
    exit(EXIT_FAILURE);
  }
#endif
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
  string in_file = fileNames[0];
  if ( !TiCC::match_back( in_file, ".ranked" ) ){
    cerr << "inputfile must have extension .ranked" << endl;
    exit(EXIT_FAILURE);
  }
  if ( !out_file.empty() ){
    if ( !TiCC::match_back( out_file, ".chained" ) )
      out_file += ".chained";
  }
  else {
    out_file = in_file + ".chained";
  }
  if ( out_file == in_file ){
    cerr << "same filename for input and output!" << endl;
    exit(EXIT_FAILURE);
  }

  ifstream input( in_file );
  if ( !input ){
    cerr << "problem opening input file: " << in_file << endl;
    exit(1);
  }

  chain_class chains( verbosity, caseless );
  string line;
  while( getline( input, line ) ){
    if ( !chains.fill( line ) ){
      cerr << "invalid line: '" << line << "'" << endl;
    }
  }
#ifdef HAVE_OPENMP
  omp_set_num_threads( 1 );
#endif
  chains.debug_info( out_file );
#ifdef HAVE_OPENMP
  omp_set_num_threads( numThreads );
#endif
  chains.output( out_file );
  cout << "results in " << out_file << endl;
}
