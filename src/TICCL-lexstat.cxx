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
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/Unicode.h"

#include "config.h"

typedef signed long int bitType;

bool verbose = false;

using namespace	std;

bitType high_five( int val ){
  bitType result = val;
  result *= val;
  result *= val;
  result *= val;
  result *= val;
  return result;
}

void create_output( string& name, const map<UChar,size_t>& chars,
		    string& orig, map<icu::UnicodeString,bitType>& hashes,
		    int clip, icu::UnicodeString& separator ){
  ofstream os( name );
  if ( !os ){
    cerr << "unable to open output file: " << name << endl;
    exit(EXIT_FAILURE);
  }
  multimap<size_t,UChar> reverse;
  map<UChar,size_t>::const_iterator it = chars.begin();
  bitType count = 0;
  bitType out_count = 0;
  while ( it != chars.end() ){
    if ( clip >= 0 && it->second < (size_t)clip ){
      out_count += it->second;
    }
    else {
      count += it->second;
      reverse.insert( make_pair(it->second,it->first) );
    }
    ++it;
  }
  os << "## Alphabetsize: " << reverse.size() << endl;
  os << "## Original file : " << orig << " with " << count
     << " accepted characters and " << out_count << " clipped characters."
     << endl;
  int start = 100;
  bitType hash = high_five( start );
  hashes.insert( make_pair( "*", hash ) );
  os << "# *\tdigits_and_punctuation\t" << hash << endl;
  start = 101;
  hash = high_five( start );
  hashes.insert( make_pair( "$", hash ) );
  os << "# $\tunknown_characters\t" << hash << endl;
  start = 102;
  int spec_cnt = 2;
  if ( !separator.isEmpty() ){
    ++spec_cnt;
    hash = high_five( start );
    hashes.insert( make_pair( separator, hash ) );
    os << separator << "\t0\t\t" << hash << endl;
    start = 103;
  }
  multimap<size_t,UChar>::const_reverse_iterator rit = reverse.rbegin();
  int out_cnt = 0;
  while ( rit != reverse.rend() ){
    hash = high_five( start );
    icu::UnicodeString us( rit->second );
    hashes.insert( make_pair( us, hash ) );
    os << us << "\t" << rit->first << "\t" << hash << endl;
    ++out_cnt;
    ++start;
    ++rit;
  }
  cerr << "created outputfile " << name
       << " with " << out_cnt << " hased lowercase characters, and " << spec_cnt
       << " special hashes for (*,$" << (separator.isEmpty()?"":","+separator)
       << ")" << endl;
}

void create_dia_file( const string& filename,
		      const map<UChar,size_t>& chars,
		      const map<icu::UnicodeString,bitType>& hashes ){
  ofstream os( filename );
  map<UChar,size_t>::const_iterator it = chars.begin();
  while ( it != chars.end() ){
    icu::UnicodeString us;
    us += it->first;
    icu::UnicodeString ss = TiCC::filter_diacritics( us );
    if ( ss != us ){
      map<icu::UnicodeString,bitType>::const_iterator hit = hashes.find( us );
      if ( hit == hashes.end() ){
	if ( verbose ){
	  cerr << "problem: " << us << " not in the hashes?" << endl;
	}
	++it;
	continue;
      }
      bitType h1 = hit->second;
      hit = hashes.find( ss );
      if ( hit == hashes.end() ){
	if ( verbose ){
	  cerr << "problem: " << ss << " not in the hashes?" << endl;
	}
	++it;
	continue;
      }
      bitType h2 = hit->second;
      bitType h;
      if ( h1 > h2 )
	h = h1 - h2;
      else
	h = h2 - h1;
      os << h << "#" << us << "~" << ss << endl;
    }
    ++it;
  }
  cerr << "created a diacritic confusion file: " << filename << endl;
}

void meld_botsing( const multimap<bitType,icu::UnicodeString>& mm, bitType h ){
  multimap<bitType,icu::UnicodeString>::const_iterator it;
  map<set<UChar>,icu::UnicodeString > ref;
  for ( it = mm.lower_bound(h); it != mm.upper_bound(h); ++it ){
    set<UChar> st;
    icu::UnicodeString s = it->second;
    for( int i=0; i < s.length(); ++i ){
      st.insert( s[i] );
    }
    ref.insert( make_pair(st,s) );
  }
  cerr << "collision at hash " << h << ": ";
  map<set<UChar>,icu::UnicodeString >::const_iterator sit=ref.begin();
  while ( sit != ref.end() ){
    cerr << sit->second << ";";
    ++sit;
  }
  cerr << endl;
}

void conditionally_insert( multimap<bitType,icu::UnicodeString>& confusions,
			   bitType key,
			   const icu::UnicodeString& value,
			   bool full ){
  if ( full ){
    confusions.insert( make_pair( key, value ) );
  }
  else if ( confusions.find( key ) == confusions.end() ){
    confusions.insert( make_pair( key, value ) );
  }
}

void generate_confusion( const string& name,
			 const map<icu::UnicodeString,bitType>& hashes,
			 int depth,
			 bool full ){
  ofstream os( name );
  if ( !os ){
    cerr << "unable to open output file: " << name << endl;
    exit(EXIT_FAILURE);
  }
  multimap<bitType,icu::UnicodeString> confusions;
  cerr << "start : " << hashes.size() << " iterations " << endl;
  map<icu::UnicodeString,bitType>::const_iterator it1 = hashes.begin();
  while ( it1 != hashes.end() ){
    cerr << "iteration " << it1->first << endl;
    // deletions/inserts of 1
    icu::UnicodeString s1 = it1->first + "~";
    conditionally_insert( confusions, it1->second, s1, full );
    map<icu::UnicodeString,bitType>::const_iterator it2 = hashes.begin();
    while ( it2 != hashes.end() ){
      if ( it2 != it1 ){
	// 1-1 substitutions
	bitType div = it2->second - it1->second;
	if ( div < 0 )
	  div = -div;
	icu::UnicodeString s2 = it1->first + "~" + it2->first;
	conditionally_insert( confusions, div, s2, full );
      }
      ++it2;
    }
    if ( depth > 1 ){
      it2 = hashes.begin();
      while ( it2 != hashes.end() ){
	// 2-0 substitutions
	icu::UnicodeString s20 = it1->first + it2->first + "~";
	bitType div = it2->second + it1->second;
	if ( div < 0 )
	  div = -div;
	conditionally_insert( confusions, div, s20, full );
	map<icu::UnicodeString,bitType>::const_iterator it3 = hashes.begin();
	while ( it3 != hashes.end() ){
	  if ( it3 != it2 && it3 != it1 ){
	    // 2-1 substitutions
	    icu::UnicodeString s = it1->first + it2->first + "~" + it3->first;
	    bitType div =  it1->second + it2->second - it3->second;
	    if ( div < 0 )
	      div = -div;
	    conditionally_insert( confusions, div, s, full );
	  }
	  if ( it2 != it1 && it3 != it1 ){
	    // 1-2 substitutions
	    icu::UnicodeString s = it1->first + "~" + it2->first + it3->first;
	    bitType div =  it1->second - it2->second - it3->second;
	    if ( div < 0 )
	      div = -div;
	    conditionally_insert( confusions, div, s, full );
	  }
	  map<icu::UnicodeString,bitType>::const_iterator it4 = hashes.begin();
	  while ( it4 != hashes.end() ){
	    if ( it3 != it1 && it3 != it2 && it4 != it1 && it4 != it2 ){
	      // 2-2 substitutions
	      icu::UnicodeString s = it1->first + it2->first + "~"
		+ it3->first + it4->first;
	      bitType div = it1->second + it2->second
		- it3->second - it4->second;
	      if ( div < 0 )
		div = -div;
	      conditionally_insert( confusions, div, s, full );
	    }
	    if ( depth > 2 ){
	      // 3-0 substitutions
	      icu::UnicodeString s30 = it1->first + it2->first + it3->first + "~";
	      bitType div = it1->second + it2->second + it3->second;
	      if ( div < 0 )
		div = -div;
	      conditionally_insert( confusions, div, s30, full );
	      if ( it4 != it1 && it4 != it2 && it4 != it3 ){
		// 3-1 substitutions
		icu::UnicodeString s = it1->first + it2->first + it3->first + "~"
		  + it4->first;
		bitType div = it4->second - it3->second
		  - it2->second - it1->second;
		if ( div < 0 )
		  div = -div;
		conditionally_insert( confusions, div, s, full );
	      }
	      if ( it2 != it1 && it3 != it1 && it4 != it1 ){
		// 1-3 substitutions
		icu::UnicodeString s = it1->first + "~"
		  + it2->first + it3->first + it4->first;
		bitType div = it1->second - it2->second - it3->second - it4->second;
		if ( div < 0 )
		  div = -div;
		conditionally_insert( confusions, div, s, full );
	      }
	      map<icu::UnicodeString,bitType>::const_iterator it5 = hashes.begin();
	      while ( it5 != hashes.end() ){
		if ( it4 != it1 && it4 != it2 && it4 != it3 &&
		     it5 != it1 && it5 != it2 && it5 != it3 ){
		  // 3-2 substitutions
		  icu::UnicodeString s = it1->first + it2->first + it3->first + "~"
		    + it4->first + it5->first;
		  bitType div = it5->second + it4->second - it3->second
		    - it2->second - it1->second;
		  if ( div < 0 )
		    div = -div;
		  conditionally_insert( confusions, div, s, full );
		}
		if ( it3 != it1 && it3 != it2 &&
		     it4 != it1 && it4 != it2 &&
		     it5 != it1 && it5 != it2 ){
		  // 2-3 substitutions
		  icu::UnicodeString s = it1->first + it2->first + "~" + it3->first
		    +it4->first + it5->first;
		  bitType div = it5->second + it4->second + it3->second
		    - it2->second - it1->second;
		  if ( div < 0 )
		    div = -div;
		  conditionally_insert( confusions, div, s, full );
		}
		map<icu::UnicodeString,bitType>::const_iterator it6 = hashes.begin();
		while ( it6 != hashes.end() ){
		  if ( it6 != it1 && it6 != it2 && it6 != it3 &&
		       it5 != it1 && it5 != it2 && it5 != it3 &&
		       it4 != it1 && it4 != it2 && it4 != it3 ){
		    // 3-3 substitutions
		    icu::UnicodeString s = it1->first + it2->first + it3->first + "~"
		      +it4->first + it5->first + it6->first;
		    bitType div = it6->second + it5->second + it4->second
		      - it3->second - it2->second - it1->second;
		    if ( div < 0 )
		      div = -div;
		    conditionally_insert( confusions, div, s, full );
		  }
		  ++it6;
		}
		++it5;
	      }
	    }
	    ++it4;
	  }
	  ++it3;
	}
	++it2;
      }
    }
    ++it1;
  }
  multimap<bitType,icu::UnicodeString>::const_iterator cit = confusions.begin();
  if ( full ){
    bitType start=0;
    set<icu::UnicodeString> unique;
    while ( cit != confusions.end() ){
      if ( cit->first != start ){
	// a new KWC starts
	if ( !unique.empty() ){
	  if ( unique.size() > 8 ){
	    meld_botsing( confusions, start );
	  }
	  os << start;
	  for ( const auto& un : unique ){
	    os << "#" << un;
	  }
	  os << endl;
	  unique.clear();
	}
	start = cit->first;
      }
      unique.insert( cit->second );
      ++cit;
    }
  }
  else {
    while ( cit != confusions.end() ){
      os << cit->first << "#" << cit->second << endl;
      ++cit;
    }
  }
  cout << "generated confusion file " << name << endl;
}

void usage( const string& name ){
  cerr << "Usage:\t" << name << " [options] dictionary" << endl;
  cerr << "\t" << name << " will create a lowercased character frequency" << endl
       << "\t\tlist from a dictionary file," << endl
       << "\t\tand a character confusion file, based on that list." << endl;
  cerr << "\t-h\t this message " << endl;
  cerr << "\t-o 'name'\t create outputfile(s) with prefix 'name'" << endl;
  cerr << "\t--diac produces an extra diacritics confusion file (extension .diac)" << endl;
  cerr << "\t--clip 'clip' truncates the character file at frequency 'clip'" << endl;
  cerr << "\t--LD depth 1, 2 or 3. (default 2) The characterlength of the confusions." << endl;
  cerr << "\t\tWhen LD=0 only a frequency list is generated." << endl;
  cerr << "\t--separator=<sep> Add the 'sep' symbol to the alphabet." << endl;
  cerr << "\t--all\tfull output. Show ALL variants in the confusions file." << endl;
  cerr << "\t\tNormally only the first is shown." << endl;
  cerr << "\t-V\tshow version " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "LD:,clip:,diac,all,separator:" );
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
  verbose = opts.extract('v');
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( argc < 2	){
    usage( progname );
    exit(EXIT_FAILURE);
  }

  string output_name;
  opts.extract( 'o', output_name );
  int depth = 2;
  string depthS = "2";
  int clip = -1;
  string clipS = "0";
  string diafile;
  bool full = opts.extract( "all" );
  bool stripdia = opts.extract( "diac" );
  string value;
  if ( opts.extract( "clip", clipS ) ){
    if ( !TiCC::stringTo(clipS,clip) ) {
      cerr << "illegal value for --clip (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( opts.extract( "LD", depthS ) ){
    if ( !TiCC::stringTo(depthS,depth) ) {
      cerr << "illegal value for --LD (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    if ( depth <0 || depth > 3 ){
      cerr << "Invalid depth --LD " << depth << endl;
      exit(EXIT_FAILURE);
    }
  }
  icu::UnicodeString separator;
  if ( opts.extract( "separator", separator ) ){
    if ( separator.length() != 1 ){
      cerr << "invalid separator, should be 1 unicode point!" << endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

  vector<string> fileNames = opts.getMassOpts();
  if ( fileNames.size() == 0 ) {
    cerr << "missing a dictionary inputfile!" << endl;
    exit( EXIT_FAILURE );
  }
  if ( fileNames.size() > 1 ) {
    cerr << "only one inputfile my be provided!" << endl;
    exit( EXIT_FAILURE );
  }
  string file_name = fileNames[0];
  ifstream is( file_name );
  if ( !is ){
    cerr << "unable to open dictionary file: " << file_name << endl;
    exit(EXIT_FAILURE);
  }
  string orig = TiCC::basename( file_name );
  if ( output_name.empty() ){
    output_name = file_name;
  }
  string lc_file_name = output_name +  ".clip" + clipS + ".lc.chars";
  string confusion_file_name = output_name + ".clip" + clipS + ".ld" + depthS + ".charconfus";
  if ( stripdia ){
    diafile = output_name + ".lc.diac";
  }

  map<UChar,size_t> lchars;
  string line;
  while ( getline( is, line ) ){
    icu::UnicodeString us = TiCC::UnicodeFromUTF8( line );
    us.toLower();
    for ( int i = 0; i < us.length(); ++i ){
      ++lchars[us[i]];
    }
  }
  cout << "done reading" << endl;
  map<icu::UnicodeString,bitType> hashes;
  create_output( lc_file_name, lchars, orig, hashes, clip, separator );
  if ( stripdia ){
    create_dia_file( diafile, lchars, hashes );
  }
  if ( depth > 0 ){
    generate_confusion( confusion_file_name, hashes, depth, full );
  }
  cout << "done!" << endl;
}
