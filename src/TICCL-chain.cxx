/*
  Copyright (c) 2006 - 2019
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
using namespace icu;

typedef signed long int bitType;
using TiCC::operator<<;

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
  chain_class( int v, bool c ): verbosity(v), caseless(c), cc_vals_present(true){};
  bool fill( const string& );
  void debug_info( const string& );
  void output( const string& );
private:
  map<string,string> heads;
  map<string, set<string>> table;
  map< string, size_t > var_freq;
  map<string,string> w_cc_conf;
  int verbosity;
  bool caseless;
  bool cc_vals_present;
};

bool chain_class::fill( const string& line ){
  vector<string> parts;
  TiCC::split_at( line, parts, "#", true );
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
    string a_word = parts[0]; // a possibly correctable word
    size_t freq1 = TiCC::stringTo<size_t>(parts[1]);
    var_freq[a_word] = freq1;
    string candidate = parts[2]; // a Correction Candidate
    size_t freq2 = TiCC::stringTo<size_t>(parts[3]);
    string cc_val;
    if ( cc_vals_present ){
      cc_val = parts[4];
      string key = a_word+candidate;
      w_cc_conf[key] = cc_val;
    }
    var_freq[candidate] = freq2;
    if ( verbosity > 3 ){
      cerr << "word=" << a_word << " CC=" << candidate << endl;
    }
    string head = heads[a_word];
    if ( head.empty() ){
      // this word does not have a 'head' yet
      if ( verbosity > 3 ){
	cerr << "word: " << a_word << " NOT in heads " << endl;
      }
      string head2 = heads[candidate];
      if ( head2.empty() ){
	// the correction candidate also has no head
	// we add it as a new head for a_word, with a table
	heads[a_word] = candidate;
	table[candidate].insert( a_word );
	if ( verbosity > 3 ){
	  cerr << "candidate : " << candidate << " not in heads too." << endl;
	  cerr << "add (" << a_word << "," << candidate << ") to heads " << endl;
	  cerr << "add " << a_word << " to table of " << candidate
	       << " ==> " << table[candidate] << endl;
	}
      }
      else {
	// the candidate knows its head already
	// add the word to the table of that head, and also register
	// the head as an (intermediate) head of a_word
	if ( verbosity > 3 ){
	  cerr << "BUT: Candidate " << candidate << " has head: "
	       << head2 << endl;
	  cerr << "add " << a_word << " to table[" << head2 << "]" << endl;
	  cerr << "AND add " << head2 << " as a head of " << a_word << endl;
	}
	heads[a_word] = head2;
	table[head2].insert( a_word );
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
	  string msg = "Error: " + a_word
	    + " has a heads entry, but no table entry!";
	  throw logic_error( msg );
	}
      }
      else {
	string msg = "Error: " + a_word
	  + " has no head entry!";
	throw logic_error( msg );
      }
    }
    return true;
  }
}

void chain_class::debug_info( const string& name ){
  string out_file = name + ".debug";
  ofstream db( out_file );
  for ( const auto& it : table ){
    db << var_freq[it.first] << " " << it.first
       << " " << it.second << endl;
  }
  cout << "debug info stored in " << out_file << endl;
}

bitType high_five( int val ){
  bitType result = val;
  result *= val;
  result *= val;
  result *= val;
  result *= val;
  return result;
}

map<UChar,bitType> alphabet;

bool fillAlpha( istream& is,
		map<UChar,bitType>& alphabet,
		int clip ){
  cout << "start reading alphabet." << endl;
  string line;
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> v;
    int n = TiCC::split_at( line, v, "\t" );
    if ( n != 3 ){
      cerr << "unsupported format for alphabet file" << endl;
      exit(EXIT_FAILURE);
    }
    int freq = TiCC::stringTo<int>( v[1] );
    if ( freq > clip || freq == 0 ){
      // freq = 0 is special, for separator
      UnicodeString v0 = TiCC::UnicodeFromUTF8( v[0] );
      bitType hash = TiCC::stringTo<bitType>( v[2] );
      alphabet[v0[0]] = hash;
    }
  }
  cout << "finished reading alphabet. (" << alphabet.size() << " characters)"
       << endl;
  return true;
}

bitType hash( const UnicodeString& s,
	      map<UChar,bitType>& alphabet ){
  static bitType HonderdEenHash = 0;
  static bitType HonderdHash = 0;
  if ( HonderdEenHash == 0 ){
    HonderdHash = high_five( 100 );
    HonderdEenHash = high_five( 101 );
  }
  UnicodeString us = s;
  us.toLower();
  bitType result = 0;
  bool multPunct = false;
  bool klets = false; //( s == "voörnaamfle" );
  for( int i=0; i < us.length(); ++i ){
    map<UChar,bitType>::const_iterator it = alphabet.find( us[i] );
    if ( it != alphabet.end() ){
      result += it->second;
      if ( klets ){
	cerr << "  CHAR, add " << UnicodeString( us[i] ) << " "
	     << it->second << " ==> " << result << endl;
      }
    }
    else {
      int8_t charT = u_charType( us[i] );
      if ( u_isspace( us[i] ) ){
	continue;
      }
      else if ( ticc_ispunct( charT ) ){
	if ( !multPunct ){
	  result += HonderdHash;
	  if ( klets ){
	    cerr << "PUNCT, add " << UnicodeString( us[i] ) << " "
		 << HonderdHash	 << " ==> " << result << endl;
	  }
	  multPunct = true;
	}
      }
      else {
	result += HonderdEenHash;
	if ( klets ){
	  cerr << "   UNK, add " << UnicodeString( us[i] ) << " "
	       << HonderdHash << " ==> " << result << endl;
	}
      }
    }
  }
  return result;
}

void chain_class::output( const string& out_file ){
  ofstream os( out_file );
  multimap<size_t, string,std::greater<size_t>> out_map;
  for ( const auto& t_it : table ){
    for ( const auto& s : t_it.second ){
      stringstream oss;
      oss << s << "#" << var_freq[s] << "#" << t_it.first
	  << "#" << var_freq[t_it.first];
      if ( cc_vals_present ){
	string val = w_cc_conf[s+t_it.first];
	if ( val.empty() ){
	  //	  cerr << "GEEN waarde voor " << s+t_it.first << endl;
	  bitType h_val = abs( ::hash(TiCC::UnicodeFromUTF8(s), alphabet )
			       - ::hash(TiCC::UnicodeFromUTF8(t_it.first), alphabet) );
	  //	  cerr << "h_val=" << h_val << endl;
	  w_cc_conf[s+t_it.first] = TiCC::toString(h_val);
	  //	  cerr << "nieuwe waarde voor " << s+t_it.first << "=" << w_cc_conf[s+t_it.first] << endl;
	}
	oss << "#" + w_cc_conf[s+t_it.first];
      }
      oss << "#" << ld( t_it.first, s, caseless ) << "#C";
      out_map.insert( make_pair( var_freq[t_it.first], oss.str() ) );
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
  cerr << "\t-o <outputfile> name of the output file." << endl;
  cerr << "\t-h or --help this message." << endl;
  cerr << "\t-v be verbose, repeat to be more verbose. " << endl;
  cerr << "\t-V or --version show version. " << endl;
  exit( EXIT_FAILURE );
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "caseless,alph:" );
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
  string alphabet_name;
  opts.extract( "alph", alphabet_name );
  if ( alphabet_name.empty() ){
    cerr << "missing --alph option" << endl;
    exit(EXIT_FAILURE);
  }
  else {
    ifstream is( alphabet_name );
    fillAlpha( is, alphabet, 0 );
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
  string line;
  while( getline( input, line ) ){
    if ( !chains.fill( line ) ){
      cerr << "invalid line: '" << line << "'" << endl;
    }
  }
  if ( verbosity > 0 ){
    chains.debug_info( out_file );
  }
  chains.output( out_file );
  cout << "results in " << out_file << endl;
}
