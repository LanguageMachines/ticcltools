/*
  Copyright (c) 2006 - 2023
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
#include "ticcl/ticcl_common.h"

#include "config.h"

bool verbose = false;

using namespace	std;
using namespace icu;
using ticcl::bitType;

void create_output( string& name, const map<UChar,size_t>& chars,
		    string& orig, map<UnicodeString,bitType>& hashes,
		    int clip, UnicodeString& separator ){
  ofstream os( name );
  if ( !os ){
    cerr << "unable to open output file: " << name << endl;
    exit(EXIT_FAILURE);
  }
  multimap<size_t,UChar> reverse;
  bitType count = 0;
  bitType out_count = 0;
  for ( const auto& it : chars ){
    if ( clip >= 0 && it.second < (size_t)clip ){
      out_count += it.second;
    }
    else {
      count += it.second;
      reverse.insert( make_pair(it.second,it.first) );
    }
  }
  os << "## Alphabetsize: " << reverse.size() + (separator.isEmpty()?0:1)
     << endl;
  os << "## Original file : " << orig << " with " << count
     << " accepted characters and " << out_count << " clipped characters."
     << endl;
  int start = 100;
  bitType hash_val = ticcl::high_five( start );
  hashes.insert( make_pair( "*", hash_val ) );
  os << "# *\tdigits_and_punctuation\t" << hash_val << endl;
  start = 101;
  hash_val = ticcl::high_five( start );
  hashes.insert( make_pair( "$", hash_val ) );
  os << "# $\tunknown_characters\t" << hash_val << endl;
  start = 102;
  int spec_cnt = 2;
  if ( !separator.isEmpty() ){
    ++spec_cnt;
    hash_val = ticcl::high_five( start );
    hashes.insert( make_pair( separator, hash_val ) );
    os << separator << "\t0\t\t" << hash_val << endl;
    start = 103;
  }
  int out_cnt = 0;
  auto rit = reverse.rbegin();
  while ( rit != reverse.rend() ){
    hash_val = ticcl::high_five( start );
    UnicodeString us( rit->second );
    hashes.insert( make_pair( us, hash_val ) );
    os << us << "\t" << rit->first << "\t" << hash_val << endl;
    ++out_cnt;
    ++start;
    ++rit;
  }
  cout << "created outputfile " << name
       << " with " << out_cnt << " hased lowercase characters, and " << spec_cnt
       << " special hashes for (*,$" << (separator.isEmpty()?"":","+separator)
       << ")" << endl;
}

void create_dia_file( const string& filename,
		      const map<UChar,size_t>& chars,
		      const map<UnicodeString,bitType>& hashes ){
  ofstream os( filename );
  for ( const auto& it : chars ){
    UnicodeString us;
    us += it.first;
    UnicodeString ss = TiCC::filter_diacritics( us );
    if ( ss != us ){
      auto hit = hashes.find( us );
      if ( hit == hashes.end() ){
	if ( verbose ){
	  cerr << "problem: " << us << " not in the hashes?" << endl;
	}
	continue;
      }
      bitType h1 = hit->second;
      hit = hashes.find( ss );
      if ( hit == hashes.end() ){
	if ( verbose ){
	  cerr << "problem: " << ss << " not in the hashes?" << endl;
	}
	continue;
      }
      bitType h2 = hit->second;
      bitType h;
      if ( h1 > h2 ){
	h = h1 - h2;
      }
      else {
	h = h2 - h1;
      }
      os << h << "#" << us << "~" << ss << endl;
    }
  }
  cout << "created a diacritic confusion file: " << filename << endl;
}

void meld_botsing( const multimap<bitType,UnicodeString>& mm, bitType h ){
  map<set<UChar>,UnicodeString > ref;
  for ( auto it = mm.lower_bound(h); it != mm.upper_bound(h); ++it ){
    set<UChar> st;
    UnicodeString s = it->second;
    for( int i=0; i < s.length(); ++i ){
      st.insert( s[i] );
    }
    ref.insert( make_pair(st,s) );
  }
  cerr << "collision at hash " << h << ": ";
  for ( const auto& sit : ref ){
    cerr << sit.second << ";";
  }
  cerr << endl;
}

void conditionally_insert( multimap<bitType,UnicodeString>& confusions,
			   bitType key,
			   const UnicodeString& value,
			   bool full ){
  if ( full ){
    confusions.insert( make_pair( key, value ) );
  }
  else if ( confusions.find( key ) == confusions.end() ){
    confusions.insert( make_pair( key, value ) );
  }
}

void generate_confusion( const string& name,
			 const map<UnicodeString,bitType>& hashes,
			 int depth,
			 bool full ){
  ofstream os( name );
  if ( !os ){
    cerr << "unable to open output file: " << name << endl;
    exit(EXIT_FAILURE);
  }
  multimap<bitType,UnicodeString> confusions;
  cout << "start : " << hashes.size() << " iterations " << endl;
  auto it1 = hashes.cbegin();
  while ( it1 != hashes.cend() ){
    cout << it1->first << " ";
    // deletions/inserts of 1
    UnicodeString s1 = it1->first + "~";
    conditionally_insert( confusions, it1->second, s1, full );
    auto it2 = hashes.cbegin();
    while ( it2 != hashes.cend() ){
      if ( it2 != it1 ){
	// 1-1 substitutions
	bitType div11;
	if ( it2->second > it1->second ){
	  div11 = it2->second - it1->second;
	}
	else {
	  div11 = it1->second - it2->second;
	}
	UnicodeString s2 = it1->first + "~" + it2->first;
	conditionally_insert( confusions, div11, s2, full );
      }
      ++it2;
    }
    if ( depth > 1 ){
      it2 = hashes.begin();
      while ( it2 != hashes.end() ){
	// 2-0 substitutions
	UnicodeString s20 = it1->first + it2->first + "~";
	bitType div20 = it2->second + it1->second;
	conditionally_insert( confusions, div20, s20, full );
	auto it3 = hashes.cbegin();
	while ( it3 != hashes.cend() ){
	  if ( it3 != it2 && it3 != it1 ){
	    // 2-1 substitutions
	    UnicodeString s = it1->first + it2->first + "~" + it3->first;
	    bitType div21;
	    bitType d1 = it1->second + it2->second;
	    bitType d2 = it3->second;
	    if ( d1 > d2 ){
	      div21 = d1 - d2;
	    }
	    else {
	      div21 = d2 - d1;
	    }
	    conditionally_insert( confusions, div21, s, full );
	  }
	  if ( it2 != it1 && it3 != it1 ){
	    // 1-2 substitutions
	    UnicodeString s = it1->first + "~" + it2->first + it3->first;
	    bitType div12;
	    bitType d1 = it1->second;
	    bitType d2 = it2->second + it3->second;
	    if ( d1 > d2 ){
	      div12 = d1 - d2;
	    }
	    else {
	      div12 = d2 - d1;
	    }
	    conditionally_insert( confusions, div12, s, full );
	  }
	  auto it4 = hashes.cbegin();
	  while ( it4 != hashes.cend() ){
	    if ( it3 != it1 && it3 != it2 && it4 != it1 && it4 != it2 ){
	      // 2-2 substitutions
	      UnicodeString s = it1->first + it2->first + "~"
		+ it3->first + it4->first;
	      bitType div22;
	      bitType d1 = it1->second + it2->second;
	      bitType d2 = it3->second + it4->second;
	      if ( d1 > d2 ){
		div22 = d1 - d2;
	      }
	      else {
		div22 = d2 - d1;
	      }
	      conditionally_insert( confusions, div22, s, full );
	    }
	    if ( depth > 2 ){
	      // 3-0 substitutions
	      UnicodeString s30 = it1->first + it2->first + it3->first + "~";
	      bitType div30 = it1->second + it2->second + it3->second;
	      conditionally_insert( confusions, div30, s30, full );
	      if ( it4 != it1 && it4 != it2 && it4 != it3 ){
		// 3-1 substitutions
		UnicodeString s = it1->first + it2->first + it3->first + "~"
		  + it4->first;
		bitType div31;
		bitType d1 = it4->second;
		bitType d2 = it3->second + it2->second + it1->second;
		if ( d1 > d2 ){
		  div31 = d1 - d2;
		}
		else {
		  div31 = d2 - d1;
		}
		conditionally_insert( confusions, div31, s, full );
	      }
	      if ( it2 != it1 && it3 != it1 && it4 != it1 ){
		// 1-3 substitutions
		UnicodeString s = it1->first + "~"
		  + it2->first + it3->first + it4->first;
		bitType div13;
		bitType d1 = it1->second;
		bitType d2 = it2->second + it3->second + it4->second;
		if ( d1 > d2 ){
		  div13 = d1 - d2;
		}
		else {
		  div13 = d2 - d1;
		}
		conditionally_insert( confusions, div13, s, full );
	      }
	      auto it5 = hashes.cbegin();
	      while ( it5 != hashes.cend() ){
		if ( it4 != it1 && it4 != it2 && it4 != it3 &&
		     it5 != it1 && it5 != it2 && it5 != it3 ){
		  // 3-2 substitutions
		  UnicodeString s = it1->first + it2->first + it3->first + "~"
		    + it4->first + it5->first;
		  bitType div32;
		  bitType d1 = it5->second + it4->second;
		  bitType d2 = it3->second + it2->second + it1->second;
		  if ( d1 > d2 ){
		    div32 = d1 - d2;
		  }
		  else {
		    div32 = d2 - d1;
		  }
		  conditionally_insert( confusions, div32, s, full );
		}
		if ( it3 != it1 && it3 != it2 &&
		     it4 != it1 && it4 != it2 &&
		     it5 != it1 && it5 != it2 ){
		  // 2-3 substitutions
		  UnicodeString s = it1->first + it2->first + "~" + it3->first
		    +it4->first + it5->first;
		  bitType div23;
		  bitType d1 = it5->second + it4->second + it3->second;
		  bitType d2 = it2->second + it1->second;
		  if ( d1 > d2 ){
		    div23 = d1 - d2;
		  }
		  else {
		    div23 = d2 - d1;
		  }
		  conditionally_insert( confusions, div23, s, full );
		}
		auto it6 = hashes.cbegin();
		while ( it6 != hashes.cend() ){
		  if ( it6 != it1 && it6 != it2 && it6 != it3 &&
		       it5 != it1 && it5 != it2 && it5 != it3 &&
		       it4 != it1 && it4 != it2 && it4 != it3 ){
		    // 3-3 substitutions
		    UnicodeString s = it1->first + it2->first + it3->first + "~"
		      +it4->first + it5->first + it6->first;
		    bitType d1 = it6->second + it5->second + it4->second;
		    bitType d2 = it3->second + it2->second + it1->second;
		    bitType div33;
		    if ( d1 > d2 ){
		      div33 = d1 - d2;
		    }
		    else {
		      div33 = d2 - d1;
		    }
		    conditionally_insert( confusions, div33, s, full );
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
  auto cit = confusions.begin();
  if ( full ){
    bitType start=0;
    set<UnicodeString> unique;
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
  cout << endl << "generated confusion file " << name << endl;
}

void usage( const string& name ){
  cerr << "Usage:\t" << name << " [options] dictionary" << endl;
  cerr << "\t" << name << " will create a lowercased character frequency" << endl
       << "\t\tlist from a dictionary file," << endl
       << "\t\tand a character confusion file, based on that list." << endl;
  cerr << "\t-o 'name'\t create outputfile(s) with prefix 'name'" << endl;
  cerr << "\t--diac\t produces an extra diacritics confusion file (extension .diac)" << endl;
  cerr << "\t--clip\t 'clip' truncates the character file at frequency 'clip'" << endl;
  cerr << "\t--LD depth\t 1, 2 or 3. (default 2) The characterlength of the confusions." << endl;
  cerr << "\t\tWhen LD=0 only a frequency list is generated." << endl;
  cerr << "\t--separator=<sep> Add the 'sep' symbol to the alphabet." << endl;
  cerr << "\t--all\tfull output. Show ALL variants in the confusions file." << endl;
  cerr << "\t\tNormally only the first is shown." << endl;
  cerr << "\t-h or --help\t this message " << endl;
  cerr << "\t-v or --verbose\t give more details during run." << endl;
  cerr << "\t-V or --version\t show version " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "LD:,clip:,diac,all,separator:,help,verbose,version" );
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
  verbose = opts.extract('v') || opts.extract("verbose");
  if ( opts.extract('V' ) || opts.extract("version") ){
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
  UnicodeString separator;
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
  if ( !TiCC::createPath( lc_file_name ) ){
    cerr << "unable to acces the path for '" << lc_file_name << "'" << endl;
    return EXIT_FAILURE;
  }
  string confusion_file_name = output_name + ".clip" + clipS + ".ld" + depthS + ".charconfus";
  if ( stripdia ){
    diafile = output_name + ".lc.diac";
  }

  map<UChar,size_t> lchars;
  UnicodeString line;
  while ( TiCC::getline( is, line ) ){
    UnicodeString us = line;
    us.toLower();
    for ( int i = 0; i < us.length(); ++i ){
      ++lchars[us[i]];
    }
  }
  cout << "done reading" << endl;
  map<UnicodeString,bitType> hashes;
  create_output( lc_file_name, lchars, orig, hashes, clip, separator );
  if ( stripdia ){
    create_dia_file( diafile, lchars, hashes );
  }
  if ( depth > 0 ){
    generate_confusion( confusion_file_name, hashes, depth, full );
  }
  cout << "done!" << endl;
}
