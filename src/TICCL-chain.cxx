/*
  Copyright (c) 2006 - 2026
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
#include "config.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace icu;
using ticcl::bitType;
using TiCC::operator<<;

const UnicodeString high_101 = TiCC::toUnicodeString(ticcl::HonderdEenHash);

unsigned int ld( const UnicodeString& in1,
		 const UnicodeString& in2,
		 bool caseless ){
  UnicodeString s1 = in1;
  UnicodeString s2 = in2;
  if ( caseless ){
    s1.toLower();
    s2.toLower();
  }
  return ticcl::ldCompare( s1, s2 );
}

map<UChar,bitType> alphabet;

class chain_class {
public:
  chain_class(): chain_class( 0,false ){};
  chain_class( int v, bool c ): verbosity(v), caseless(c), cc_vals_present(true){};
  bool fill( const UnicodeString&, bool );
  void debug_info( ostream& );
  void output( const string& );
  UnicodeString top_head( const UnicodeString& );
  void final_merge();
 private:
  map<UnicodeString,UnicodeString> heads;
  map<UnicodeString, set<UnicodeString>> table;
  map<UnicodeString, size_t > var_freq;
  map<UnicodeString, UnicodeString> w_cc_conf;
  set<UnicodeString> processed;
  int verbosity;
  bool caseless;
  bool cc_vals_present;
};

UnicodeString chain_class::top_head( const UnicodeString& candidate ){
  UnicodeString result = heads[candidate];
  if ( !result.isEmpty() ){
    UnicodeString next = top_head( result );
    if ( !next.isEmpty() ){
      result = next;
    }
  }
  return result;
}

void chain_class::final_merge(){
  for ( auto& [word,word_set] : table ){
    if ( !word_set.empty() ){
      // for all entries that seem to be a 'head'
      UnicodeString head = top_head( word );
      assert( head != word );
      if ( !head.isEmpty() ){
	// so it has a higher head
	if ( verbosity > 3 ){
	  cerr << "merge: " << word << word_set << " into "
	       << head << table[head] << endl;
	}
	for ( const auto& s : word_set ){
	  table[head].insert( s );
	  heads[s] = head;
	}
	word_set.clear();
      }
    }
  }
}

UChar diff_char( const UnicodeString& in1, const UnicodeString& in2 ){
  UnicodeString s1 = in1;
  UnicodeString s2 = in2;
  s1.toLower();
  s2.toLower();
  UChar result = alphabet.begin()->first;
  for ( int i=0; i < s2.length(); ++i ){
    if ( s1.indexOf(s2[i] ) == -1 ){
      result = s2[i];
      break;
    }
  }
  return result;
}

bool chain_class::fill( const UnicodeString& line, bool nounk ){
  vector<UnicodeString> parts = TiCC::split_exact_at( line, "#" );
  if ( parts.size() < 6
       || parts.size() > 7 ){
    return false;
  }
  else {
    if ( parts.size() == 6 ){
      cc_vals_present = false;
    }
    else if ( cc_vals_present == false ){
      cerr << "conflicting data in chained file, didn't expect cc_val entries" << endl;
      exit(EXIT_FAILURE);
    }
    UnicodeString a_word = parts[0]; // a possibly correctable word
    if ( processed.find(a_word) != processed.end() ){
      // we have already seen this word. probably ranked with a clip >1
      // just ignore!
      //      cerr << "ignore extra entry for: " << a_word << endl;
      return true;
    }
    else {
      processed.insert(a_word);
      // so a new word with Correction Candidate
      size_t freq1 = TiCC::stringTo<size_t>(parts[1]);
      UnicodeString candidate = parts[2];
      // a Correction Candidate
      size_t freq2 = TiCC::stringTo<size_t>(parts[3]);
      if ( cc_vals_present ){
	UnicodeString cc_val = parts[4];
	if ( nounk && cc_val == high_101 ){
	  //	  cerr << "diff?? " << a_word << " " << candidate << endl;
	  // one character difference
	  if ( candidate.length() > a_word.length() ){
	    UChar diff = diff_char( a_word, candidate );
	    //	    cerr << "diff=" << diff << endl;
	    if ( alphabet.find( diff ) == alphabet.end() ){
	      // this will not do. ignore!
	      cerr << "skip " << a_word << " --> " << candidate << endl;
	      return true;
	    }
	  }
	}
	UnicodeString key = a_word+candidate;
	w_cc_conf[key] = cc_val;
      }
      var_freq[a_word] = freq1;
      var_freq[candidate] = freq2;
      if ( verbosity > 3 ){
	cerr << endl << "word=" << a_word << " CC=" << candidate << endl;
      }
      UnicodeString head = heads[a_word];
      if ( head.isEmpty() ){
	// this word does not have a 'head' yet
	if ( verbosity > 3 ){
	  cerr << "word: " << a_word << " NOT in heads " << endl;
	}
	UnicodeString head2 = heads[candidate];
	if ( head2.isEmpty() ){
	  // the correction candidate also has no head
	  // we add it as a new head for a_word, with a table
	  heads[a_word] = candidate;
	  table[candidate].insert( a_word );
	  if ( verbosity > 3 ){
	    cerr << "candidate : " << candidate << " not in heads too." << endl;
	    cerr << "add " << candidate << " to heads[" << a_word << "]" << endl;
	    cerr << "add " << a_word << " to table of " << candidate
		 << " ==> " << table[candidate] << endl;
	  }
	}
	else {
	  // the candidate knows its head already
	  // add the word to the table of that head, and also register
	  // the head as an (intermediate) head of a_word
	  heads[a_word] = head2;
	  table[head2].insert( a_word );
	  if ( verbosity > 3 ){
	    cerr << "BUT: Candidate " << candidate << " has head: "
		 << head2 << endl;
	    cerr << "add " << a_word << " to table[" << head2 << "]"
		 << " ==> " << table[head2] << endl;
	    cerr << "AND add " << head2 << " as a head of " << a_word << endl;
	  }
	}
      }
      else {
	// the word has a head
	if ( verbosity > 3 ){
	  cerr << "word: " << a_word << " IN heads " << head << endl;
	}
	auto const tit = table.find( head );
	if ( tit != table.end() ){
	  // there MUST be some candidates registered for the head
	  if ( verbosity > 3 ){
	    cerr << "lookup " << a_word << " in " << tit->second << endl;
	  }
	  if ( tit->second.find( a_word ) == tit->second.end() ){
	    string msg = "Error: " + TiCC::UnicodeToUTF8(a_word)
	      + " has a heads entry, but no table entry!";
	    throw logic_error( msg );
	  }
	}
	else {
	  string msg = "Error: " +  TiCC::UnicodeToUTF8(a_word)
	    + " has no head entry!";
	  throw logic_error( msg );
	}
      }
    }
    if ( verbosity > 4 ){
      cerr << endl;
      debug_info( cerr );
    };
    return true;
  }
}

void chain_class::debug_info( ostream& db ){
  for ( const auto& [word,head] : heads ){
    db << "head[" << word << "]=" << head << endl;
  }
  for ( const auto& [word,word_set] : table ){
    db << var_freq[word] << " " << word
       << " " << word_set << endl;
  }
}

void chain_class::output( const string& out_file ){
  ofstream os( out_file );
  multimap<size_t, string,std::greater<size_t>> out_map;
  for ( const auto& [word,word_set] : table ){
    for ( const auto& s : word_set ){
      stringstream oss;
      oss << s << "#" << var_freq[s] << "#" << word
	  << "#" << var_freq[word];
      if ( cc_vals_present ){
	UnicodeString val = w_cc_conf[s+word];
	if ( val.isEmpty() ){
	  //	  cerr << "GEEN waarde voor " << s+word << endl;
	  bitType h1 = ticcl::hash(s, alphabet );
	  bitType h2 = ticcl::hash(word, alphabet );
	  bitType h_val;
	  if ( h1 > h2 ){
	    h_val = h1 - h2;
	  }
	  else {
	    h_val = h2 - h1;
	  }
	  //	  cerr << "h_val=" << h_val << endl;
	  w_cc_conf[s+word] = TiCC::toUnicodeString(h_val);
	  //	  cerr << "nieuwe waarde voor " << s+word << "=" << w_cc_conf[s+word] << endl;
	}
	oss << "#" + w_cc_conf[s+word];
      }
      oss << "#" << ld( word, s, caseless ) << "#C";
      out_map.insert( make_pair( var_freq[word], oss.str() ) );
    }
  }
  for ( const auto& t_it : out_map ){
    os << t_it.second << endl;
  }
}

void usage( const string& name ){
  cerr << "usage: " << name << endl;
  cerr << "\t--caseless Calculate the Levensthein (or edit) distance ignoring case." << endl;
  cerr << "\t--alph <alphafile> name of the alphabet file." << endl;
  cerr << "\t--nounk Skip Correction Candidates that are equal except for one extra UNK character." << endl;
  cerr << "\t-o <outputfile> name of the output file." << endl;
  cerr << "\t-h or --help this message." << endl;
  cerr << "\t-v be verbose, repeat to be more verbose. " << endl;
  cerr << "\t-V or --version show version. " << endl;
  exit( EXIT_FAILURE );
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.add_short_options( "vVho:" );
    opts.add_long_options( "caseless,alph:,nounk" );
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
  bool nounk = opts.extract( "nounk" );
  string alphabet_name;
  opts.extract( "alph", alphabet_name );
  if ( alphabet_name.empty() ){
    cerr << "missing --alph option" << endl;
    exit(EXIT_FAILURE);
  }
  else {
    ifstream is( alphabet_name );
    cout << "start reading alphabet: " << alphabet_name << endl;
    ticcl::fillAlphabet( is, alphabet, 0 );
    cout << "finished reading alphabet. (" << alphabet.size() << " characters)"
	 << endl;
  }
  string out_file;
  opts.extract( 'o', out_file );

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
  UnicodeString line;
  while( TiCC::getline( input, line ) ){
    if ( !chains.fill( line, nounk ) ){
      cerr << "invalid line: '" << line << "'" << endl;
    }
  }
  chains.final_merge();
  if ( verbosity > 0 ){
    string db_file = out_file + ".debug";
    ofstream db( db_file );
    chains.debug_info( db );
    cout << endl << "debug info stored in " << out_file << endl;
  }
  chains.output( out_file );
  cout << "results in " << out_file << endl;
}
