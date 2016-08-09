/*
  Copyright (c) 2006 - 2016
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
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcl/unicode.h"

#include "config.h"

using namespace	std;

const string SEPARATOR = "_";

typedef unsigned long int bitType;

bitType high_five( int val ){
  bitType result = val;
  result *= val;
  result *= val;
  result *= val;
  result *= val;
  return result;
}

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
    if ( freq > clip ){
      UnicodeString us = UTF8ToUnicode( v[0] );
      bitType hash = TiCC::stringTo<bitType>( v[2] );
      alphabet[us[0]] = hash;
    }
  }
  cout << "finished reading alphabet." << endl;
  return true;
}

bitType hash( const string& s,
	      map<UChar,bitType>& alphabet ){
  static bitType HonderdEenHash = 0;
  static bitType HonderdHash = 0;
  if ( HonderdEenHash == 0 ){
    HonderdHash = high_five( 100 );
    HonderdEenHash = high_five( 101 );
  }
  UnicodeString us = UTF8ToUnicode( s );
  us.toLower();
  bitType result = 0;
  bool multPunct = false;
  bool klets = false; // ( s == "Engeĳclíe" );
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

void create_output( ostream& os,
		    map<bitType, set<string> >& anagrams ){
  for ( const auto& it : anagrams ){
    bitType val = it.first;
    os << val << "~";
    for ( auto const&  s : it.second ){
      os << s;
      if ( &s != &(*it.second.crbegin()) )
	os << "#";
    }
    os << endl;
  }
  os << endl;
}

string filter_tilde_hashtag( const string& w ){
  // assume that we cannot break UTF8 by replacing # or ~ by _
  string result = w;
  for ( size_t i=0; i < result.length(); ++i ){
    if ( result[i] == '~' || result[i] == '#' ){
      result[i] = '_';
    }
  }
  return result;
}

void usage( const string& name ){
  cerr << "usage:" << name << " [options] <clean frequencyfile>" << endl;
  cerr << "\t" << name << " will read a wordfrequency list (in FoLiA-stats format) " << endl;
  cerr << "\t\t which is assumed to be 'clean', but may consiste of n-grams with different arities." << endl;
  cerr << "\t\t The output will be an anagram hash file." << endl;
  cerr << "\t\t When a background corpus is specified, we also produce" << endl;
  cerr << "\t\t a new (merged) frequency file. " << endl;
  cerr << "\t--alph='file'\t name of the alphabet file" << endl;
  cerr << "\t--background='file'\t name of the background corpus" << endl;
  cerr << "\t--clip=<clip> : cut off of the alphabet." << endl;
  cerr << "\t-h or --help\t this message " << endl;
  cerr << "\t--artifrq='value': if value > 0, create a separate list of anagram" << endl;
  cerr << "\t\t values that don't have the lexical frequency 'artifrq' " << endl;
  cerr << "\t\t for n-grams, only those n-grams are written where at least one" << endl;
  cerr << "\t\t of the composing parts does not have the lexical frequency artifrq. " << endl;
  cerr << "\t--ngrams\t When the frequency file contains n-grams. (not necessary of equal arity)" << endl;
  cerr << "\t\t we split them into 1-grams and do a frequency lookup per part for the artifreq value." << endl;
  cerr << "\t-V ot --version\t show version " << endl;
  cerr << "\t-v\t verbose (not used yet) " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVh" );
    opts.set_long_options( "alph:,background:,artifrq:,clip:,help,version,ngrams" );
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
  string alphafile;
  string backfile;
  int clip = 0;
  size_t artifreq = 0;
  if ( opts.extract('h' ) || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V') || opts.extract("version") ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  bool verbose = opts.extract( 'v' );
  opts.extract( "alph", alphafile );
  opts.extract( "background", backfile );
  string value;
  if ( opts.extract( "clip", value ) ){
    if ( !TiCC::stringTo(value,clip) ) {
      cerr << "illegal value for --clip (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  bool do_ngrams = opts.extract( "ngrams" );
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  if ( verbose ){
    cout << "artifrq= " << artifreq << endl;
  }
  vector<string> fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing input file" << endl;
    exit(EXIT_FAILURE);
  }
  else if ( fileNames.size() > 1 ){
    cerr << "only one input file is possible." << endl;
    cerr << "but found: " << endl;
    for ( const auto& s : fileNames ){
      cerr << s << endl;
    }
    exit(EXIT_FAILURE);
  }
  string file_name = fileNames[0];
  ifstream is( file_name.c_str() );
  if ( !is ){
    cerr << "unable to open corpus frequency file: " << file_name << endl;
    exit(EXIT_FAILURE);
  }
  bool doMerge = false;
  if ( !backfile.empty() ){
    ifstream bs( backfile.c_str() );
    if ( !bs ){
      cerr << "unable to open background frequency file: " << backfile << endl;
      exit(EXIT_FAILURE);
    }
    doMerge = true;
  }
  if ( alphafile.empty() ){
    cerr << "We need an alphabet file!" << endl;
    exit(EXIT_FAILURE);
  }
  ifstream as( alphafile.c_str() );
  if ( !as ){
    cerr << "unable to open alphabet file: " << alphafile << endl;
    exit(EXIT_FAILURE);
  }
  map<UChar,bitType> alphabet;
  if ( !fillAlpha( as, alphabet, clip ) ){
    cerr << "serious problems reading alphabet file: " << alphafile << endl;
    exit(EXIT_FAILURE);
  }

  string out_file_name = file_name + ".anahash";
  ofstream os( out_file_name.c_str() );
  if ( !os ){
    cerr << "unable to open output file: " << out_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream fos;
  string foci_file_name = file_name + ".corpusfoci";
  if ( artifreq > 0 ){
    fos.open( foci_file_name.c_str() );
    if ( !fos ){
      cerr << "unable to open foci file: " << foci_file_name << endl;
      exit(EXIT_FAILURE);
    }
  }

  map<string,bitType> merged;
  map<bitType, set<string> > anagrams;
  set<string> focus_words;
  cout << "start hashing from the corpus frequency file." << endl;
  string line;
  while ( getline( is, line ) ){
    vector<string> v;
    int n = TiCC::split_at( line, v, "\t" );
    if ( n != 2 ){
      cerr << "frequency file in wrong format!" << endl;
      cerr << "offending line: " << line << endl;
      exit(EXIT_FAILURE);
    }
    string word = filter_tilde_hashtag( v[0] );
    bitType h = ::hash( word, alphabet );
    anagrams[h].insert( word );
    if ( artifreq > 0 ){
      bitType freq = TiCC::stringTo<bitType>( v[1] );
      if ( freq >= artifreq ){
	focus_words.insert( word );
      }
      if ( doMerge ){
	merged[v[0]] = freq;
      }
    }
  }
  set<bitType> foci;
  if ( artifreq > 0 ){
    is.clear();
    is.seekg(ios_base::beg);
    string line;
    while ( getline( is, line ) ){
      vector<string> v;
      TiCC::split_at( line, v, "\t" );
      string word = filter_tilde_hashtag( v[0] );
      bitType h = ::hash( word, alphabet );
      if ( do_ngrams ){
	vector<string> parts;
	if ( TiCC::split_at( word, parts, SEPARATOR ) ){
	  bool accept = false;
	  for ( auto const& part: parts ){
	    if ( focus_words.find(part) == focus_words.end() ){
	      accept = true;
	    }
	  }
	  if ( accept ){
	    foci.insert( h );
	  }
	}
      }
      else {
	if ( focus_words.find(word) == focus_words.end() ){
	  foci.insert( h );
	}
      }
    }
  }
  if ( artifreq > 0 ){
    cout << "generating foci file: " << foci_file_name << " with " << foci.size() << " entries" << endl;
    for ( const auto& f : foci ){
      fos << f << endl;
    }
  }
  if ( doMerge ){
    cerr << "merge background corpus: " << backfile << endl;
    ifstream bs( backfile.c_str() );
    while ( getline( bs, line ) ){
      vector<string> v;
      int n = TiCC::split_at( line, v, "\t" );
      if ( n != 2 ){
	cerr << "background file in wrong format!" << endl;
	cerr << "offending line: " << line << endl;
	exit(EXIT_FAILURE);
      }
      string word = filter_tilde_hashtag(v[0] );
      bitType h = ::hash( word, alphabet );
      anagrams[h].insert( word );
      bitType freq = TiCC::stringTo<bitType>( v[1] );
      merged[v[0]] += freq;
    }
    string merge_file_name = file_name + ".merged";
    ofstream ms( merge_file_name.c_str() );
    for ( const auto& it : merged ){
      ms << it.first << "\t" << it.second << endl;
    }
    cerr << "stored merged corpus in " << merge_file_name << endl;

  }

  cout << "generating output file: " << out_file_name << endl;
  create_output( os, anagrams );

  cout << "done!" << endl;
}
