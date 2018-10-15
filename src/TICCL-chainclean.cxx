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

set<string> follow_words;

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
  vector<string> cc_parts;
  string cc_freq;
  string ld;
};

ostream& operator<<( ostream& os, const record& rec ){
  os << rec.variant << "#" << rec.v_freq << "#" << rec.cc << "#" << rec.cc_freq << "#" << rec.ld << "#C";
  return os;
}

ostream& operator<<( ostream& os, const record *rec ){
  if ( rec ){
    os << *rec;
  }
  else {
    os << "NULL";
  }
  return os;
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "lexicon:,artifreq:,follow:" );
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

  while ( opts.extract( "follow", value ) ){
    vector<string> parts = TiCC::split_at( value, "," );
    for ( const auto& p : parts ){
      follow_words.insert( p );
    }
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
    rec.cc_parts = TiCC::split_at( rec.cc, SEPARATOR );
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
  //
  // Maybe (but why) sorting parts_freq on freqency is needed?
  //
  multimap<int,string,std::greater<int>> desc_parts_freq;
  // sort on highest frequency first.
  // DOES IT REALLY MATTER???
  for ( const auto& cc : parts_freq ){
    desc_parts_freq.insert( make_pair(cc.second,cc.first) );
  }
  if ( verbosity > 0 ){
    cerr << "The unknown parts:" << endl;
    for ( const auto& it : desc_parts_freq ){
      cerr << it.first << "\t" << it.second << endl;
    }
  }

  list<record*> copy_records;
  for ( auto& rec : records ){
    copy_records.push_back( &rec );
  }
  bool show = false;
  set<record*> done_records;
  map<string,string> done;
  for ( const auto& part : desc_parts_freq ) {
    show = (verbosity>0 )
    || follow_words.find( part.second ) != follow_words.end();
    if ( show ){
      cerr << "\n  Loop for part: " << part.second << endl;
    }
    map<string,int> cc_freqs;
    auto it = records.begin();
    while ( it != records.end() ){
      if ( it->v_parts.size() != 1 ){
	bool match = false;
	for ( const auto& p : it->v_parts ){
	  if ( p == part.second ){
	    match = true;
	    break;
	  }
	}
	if ( match ){
	  for ( const auto& cp : it->cc_parts ){
	    ++cc_freqs[cp];
	    if ( show ){
	      cerr << "for: " << part.second << " increment " << cp << endl;
	    }
	  }
	}
      }
      ++it;
    }
    multimap<int,string,std::greater<int>> desc_cc;
    // sort on highest frequency first.
    // DOES IT REALLY MATTER???
    for ( const auto& cc : cc_freqs ){
      desc_cc.insert( make_pair(cc.second,cc.first) );
    }
    if ( verbosity > 0 ){
      cerr << "found " << desc_cc.size() << " CC's for: " << part.second << endl;
      for ( const auto& it : desc_cc ){
	cerr << it.first << "\t" << it.second << endl;
      }
    }
    for ( const auto& dcc : desc_cc ){
      if ( show ){
	cerr << "BEKIJK: " << dcc.second << "[" << dcc.first << "]" << endl;
      }
      auto it = copy_records.begin();
      while ( it != copy_records.end() ){
	record* rec = *it;
	if ( !rec ){
	  ++it;
	  continue;
	}
	if ( done_records.find( rec ) != done_records.end() ){
	  ++it;
	  continue;
	}
	string key = rec->variant + rec->cc;
	if ( rec->v_parts.size() > 1 ){
	  bool local_show = verbosity > 0;
	  for ( const auto& p : rec->v_parts ){
	    local_show |= follow_words.find( p ) != follow_words.end();
	  }
	  // if ( local_show ){
	  //   cerr << "bekijk met " << dcc.second << ":" << rec << endl;
	  // }
	  bool match = false;
	  for( const auto& cp : rec->cc_parts ){
	    if ( dcc.second == cp ){
	      // CC match
	      for ( const auto& p : rec->v_parts ){
		if ( p == part.second ){
		  // variant match too
		  match = true;
		  break;
		}
	      }
	      if ( match ){
		if ( local_show ){
		  cerr << "both " << cp << " and " << part.second
		       << " matched in: " << rec << endl;
		}
		if ( done.find( cp ) != done.end() ){
		  string v = done[cp];
		  if ( rec->variant.find(v ) != string::npos ){
		    if ( local_show ){
		      cerr << "IGNORE: " << rec << endl;
		    }
		    *it = 0;
		  }
		  else {
		    if ( local_show ){
		      cerr << "INSERT: " << rec << endl;
		    }
		    done[rec->cc] = rec->variant;
		    done_records.insert(rec);
		  }
		}
		else {
		  if ( local_show ){
		    cerr << "INSERT: " << rec << endl;
		  }
		  done[rec->cc] = rec->variant;
		  done_records.insert(rec);
		}
		break;
	      }
	    }
	  }
	}
	else {
	  // if ( rec->variant == part.second ){
	  //   if ( show ){
	  //     cerr << "remove translation of unknown part: " << part.second
	  // 	   << " in " << rec << endl;
	  //   }
	  //   *it = 0;
	  // }
	  done[rec->cc] = rec->variant;
	  done_records.insert(rec);
	}
	++it;
      }
    }
  }
  ofstream os( out_name );
  cerr << "copy_records.size()= " << copy_records.size() << endl;
  int count = 0;
  for ( const auto it : copy_records ){
    if ( it != 0 ){
      ++count;
      os << it << endl;
    }
  }
  cerr << "copy_records.count()= " << count << endl;
  cout << "results in " << out_name << endl;
}
