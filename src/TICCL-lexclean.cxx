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
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/FileUtils.h"
#include "ticcl/unicode.h"

#include "config.h"

using namespace	std;
using namespace	TiCC;

void create_wf_list( const map<string, unsigned int>& wc,
		     const string& filename, unsigned int totalIn,
		     bool doperc ){
  unsigned int total = totalIn;
  ofstream os( filename.c_str() );
  if ( !os ){
    cerr << "failed to create outputfile '" << filename << "'" << endl;
    exit(EXIT_FAILURE);
  }
  map<unsigned int, set<string> > wf;
  for ( const auto& cit : wc ){
    wf[cit.second].insert( cit.first );
  }
  unsigned int sum=0;
  unsigned int types=0;
  map<unsigned int, set<string> >::const_reverse_iterator wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for ( const auto& s : wit->second ){
      sum += wit->first;
      os << s << "\t" << wit->first;
      if ( doperc ){
	os << "\t" << sum << "\t" << std::setprecision(8) << 100 * double(sum)/total;
      }
      os << endl;
      ++types;
    }
    ++wit;
  }
  cout << "created cleaned list '" << filename << "'" << endl;
  cout << "with " << total << " words." << endl;
}

void dump_quarantine( const string& filename,
		      const map<string, unsigned int>& qw ){
  ofstream os( filename.c_str() );
  if ( !os ){
    cerr << "failed to create outputfile '" << filename << "'" << endl;
    exit(EXIT_FAILURE);
  }
  for ( const auto& it : qw ){
    os << it.first << "\t" << it.second << endl;
  }
  cout << "created quarantine list '" << filename << "'" << endl;
  cout << "with " << qw.size() << " items. " << endl;
}

bool isClean( const string& s, const set<UChar>& alp ){
  UnicodeString us = UTF8ToUnicode( s );
  for ( int i=0; i < us.length(); ++i ){
    if ( alp.find( us[i] ) == alp.end() )
      return false;
  }
  return true;
}

bool fillAlpha( const string& file, set<UChar>& alphabet ){
  string line;
  ifstream is( file.c_str() );
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> v;
    int n = TiCC::split( line, v );
    if ( n != 3 ){
      cerr << "unsupported format for alphabet file" << endl;
      exit(EXIT_FAILURE);
    }
    UnicodeString us = UTF8ToUnicode( v[0] );
    us.toLower();
    alphabet.insert( us[0] );
    us.toUpper();
    alphabet.insert( us[0] );
    // for now, we don't use the other fields
  }
  return true;
}


void usage(){
  cerr << "Usage: [options] file/dir" << endl;
  cerr << "\t-a\t alphabet file" << endl;
  cerr << "\t-h\t this message" << endl;
  cerr << "\t-t\t assume a POS tagged input" << endl;
  cerr << "\t-V\t show version " << endl;
  cerr << "\t FoLiA-lexclean will clean FoLiA-lexstat files" << endl;
  cerr << "\t The output will be a 2 or 4 columned tab sparated file, extension: .clean " << endl;
  cerr << "\t\t (4 columns when -p is specified)" << endl;
  cerr << "\t-p\t output percentages too. " << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hVpa:t", "" );
  try {
    opts.init(argc,argv);
  }
  catch( OptionError& e ){
    cerr << e.what() << endl;
    usage();
    exit( EXIT_FAILURE );
  }

  bool dopercentage = false;
  bool postagged = false;
  bool mood;
  string value;
  string alpha;
  if ( opts.extract('V', value, mood ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('h', value, mood ) ){
    usage();
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('a', value, mood ) ){
    alpha = value;
  }
  if ( opts.extract('p', value, mood ) ){
    dopercentage = true;
  }
  if ( opts.extract('t', value, mood ) ){
    postagged = true;
  }
  if ( !opts.empty() ){
    usage();
    exit(EXIT_FAILURE);
  }

  set<UChar> alphabet;
  if ( !alpha.empty() ){
    fillAlpha( alpha, alphabet );
    cerr << "read alphabet file with " << alphabet.size()
	 << " characters" << endl;
  }

  vector<string> fileNames = opts.getMassOpts();

  size_t toDo = fileNames.size();
  if ( toDo == 0 ){
    cerr << "no matching files found" << endl;
    exit(EXIT_SUCCESS);
  }
  map<string,unsigned int> wc;
  map<string,unsigned int> qw;
  for ( size_t fn=0; fn < fileNames.size(); ++fn ){
    string docName = fileNames[fn];
    ifstream is( docName.c_str() );
    string line;
    unsigned int word_total = 0;
    while ( getline( is, line ) ){
      vector<string> vec;
      int num = split_at( line, vec, "\t" );
      if ( num != 2 && num != 4 ){
	cerr << "unexpected line: '" << line << "' in " << docName << endl;
	continue;
      }
      string val = vec[0];
      vector<string> v2;
      if ( postagged ){
	if ( split( val, v2 ) > 1 ){
	  val = v2[0];
	}
	else {
	  cerr << "pos tagged files need a space separated value in the first column" << endl;
	  exit(EXIT_FAILURE);
	}
      }
      unsigned int freq = stringTo<unsigned int>( vec[1] );
      if ( isClean( val, alphabet ) ){
	wc[vec[0]] = freq;
	word_total += freq;
      }
      else {
	qw[vec[0]] = freq;
      }
    }
    string outname = docName + ".clean";
    create_wf_list( wc, outname, word_total, dopercentage );
    outname = docName + ".dirty";
    dump_quarantine( outname, qw );
  }
}
