/*
  Copyright (c) 2006 - 2018
  CLST  - Radboud University
  ILK   - Tilburg University

  This file is part of ticcltools

  foliatools is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  foliatools is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

  For questions and suggestions, see:
      https://github.com/LanguageMachines/frog/issues
  or send mail to:
      lamasoftware (at ) science.ru.nl

*/

#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/XMLtools.h"
#include "ticcutils/Unicode.h"

#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace	std;
using namespace	TiCC;

bool verbose = false;

void create_wf_list( const map<string, unsigned int>& wc,
		     const string& filename, unsigned int total_in, bool doperc ){
  ofstream os( filename );
  if ( !os ){
    cerr << "failed to create outputfile '" << filename << "'" << endl;
    exit(EXIT_FAILURE);
  }
  map<unsigned int, set<string> > wf;
  for( const auto& cit : wc ){
    wf[cit.second].insert( cit.first );
  }
  unsigned int sum=0;
  unsigned int types=0;
  map<unsigned int, set<string> >::const_reverse_iterator wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for( const auto& sit : wit->second ){
      sum += wit->first;
      os << sit << "\t" << wit->first;
      if ( doperc ){
	os << "\t" << sum << "\t" << 100 * double(sum)/total_in;
      }
      os << endl;
      ++types;
    }
    ++wit;
  }
#pragma omp critical
  {
    cout << "created WordFreq list '" << filename << "'" << endl
	 << "with " << total_in << " tokens and " << types
	 << " types. TTR= " << (double)types/total_in
	 << ", the angle is " << atan((double)types/total_in)*180/M_PI
	 << " degrees" << endl;
  }
}

size_t read_words( const string& doc_name, map<string,unsigned int>& wc ){
  size_t word_total = 0;
  ifstream is( doc_name );
  string line;
  while ( getline( is, line ) ){
    vector<string> v;
    if ( TiCC::split( line, v ) < 2 ){
      cerr << "invalid input: " << line << endl;
      continue;
    }
    string wrd = v[0];
    size_t frq = stringTo<int>(v[1]);
#pragma omp critical
    {
      wc[wrd] += frq;
    }
    word_total += frq;
  }
  return word_total;
}


void usage( const string& name ){
  cerr << "Usage: " << name << " [options] file/dir" << endl;
  cerr << "\t TICCL-mergelex will create a merged lexicond from range of" << endl;
  cerr << "\t lexicon files in TICCL-lexstat or FoLiA-stats format." << endl;
  cerr << "The output will be a 2 or 4 columned tab separated file, extension: *tsv " << endl;
  cerr << "\t (4 columns when -p is specified)" << endl;
  cerr << "\t-p\t output percentages too. " << endl;
  cerr << "\t-t\t number_of_threads" << endl;
  cerr << "\t-h\t this message" << endl;
  cerr << "\t-v\t very verbose output." << endl;
  cerr << "\t-V\t show version " << endl;
  cerr << "\t-e\t expr: specify the expression all input files should match with." << endl;
  cerr << "\t-o\t name of the output file(s) prefix." << endl;
  cerr << "\t-R\t search the dirs recursively (when appropriate)." << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hVve:t:o:Rp", "" );
  try {
    opts.init(argc,argv);
  }
  catch( OptionError& e ){
    cerr << e.what() << endl;
    usage(argv[0]);
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( argc < 2 ){
    usage( progname );
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_OPENMP
  int num_of_threads = 1;
#endif
  string expression;
  string out_prefix;
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('h' ) ){
    usage(progname);
    exit(EXIT_SUCCESS);
  }
  verbose = opts.extract( 'v' );
  bool recursiveDirs = opts.extract( 'R' );
  if ( !opts.extract( 'o', out_prefix ) ){
    cerr << "an output filename prefix is required. (-o option) " << endl;
    exit(EXIT_FAILURE);
  }
  string value;
  if ( opts.extract('t', value ) ){
#ifdef HAVE_OPENMP
    if ( !stringTo(value, num_of_threads ) ){
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
#else
    cerr << "OpenMP support is missing. -t option is not supported" << endl;
    exit( EXIT_FAILURE );
#endif
  }
  opts.extract('e', expression );
  bool dopercentage = opts.extract('p');
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_OPENMP
  if ( num_of_threads != 1 )
    omp_set_num_threads( num_of_threads );
#endif

  vector<string> mass_opts = opts.getMassOpts();
  if ( mass_opts.empty() ){
    cerr << "no file or dir specified!" << endl;
    exit(EXIT_FAILURE);
  }
  vector<string> file_names;
  string dir_name;
  if ( mass_opts.size() > 1 ){
    // assume a list of files
    file_names = mass_opts;
  }
  else {
    dir_name = mass_opts[0];
    try {
      file_names = searchFilesMatch( dir_name, expression, recursiveDirs );
    }
    catch ( const exception& e ){
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
  }
  size_t to_do = file_names.size();
  if ( to_do == 0 ){
    cerr << "no matching files found" << endl;
    exit(EXIT_SUCCESS);
  }

  string::size_type pos = out_prefix.find( "." );
  if ( pos != string::npos && pos == out_prefix.length()-1 ){
    // outputname ends with a .
    out_prefix = out_prefix.substr(0,pos);
  }
  pos = out_prefix.find( "/" );
  if ( pos != string::npos && pos == out_prefix.length()-1 ){
    // outputname ends with a /
    out_prefix += "ticclstats";
  }

  if ( to_do > 1 ){
    cout << "start processing of " << to_do << " files " << endl;
  }
  map<string,unsigned int> wc;
  unsigned int word_total =0;
#pragma omp parallel for shared(file_names,word_total,wc)
  for ( size_t fn=0; fn < file_names.size(); ++fn ){
    string doc_name = file_names[fn];
    unsigned int word_count = read_words( doc_name, wc );
    word_total += word_count;
#pragma omp critical
    {
      cout << "Processed :" << doc_name << " with " << word_count << " words,"
	   << " still " << --to_do << " files to go." << endl;
    }
  }
  if ( !dir_name.empty() ){
    cout << "done processsing directory '" << dir_name << "' in total "
	 << word_total << " words were found." << endl;
  }
  cout << "start outputting the results" << endl;
  string file_name = out_prefix + ".wordfreqlist.tsv";
  create_wf_list( wc, file_name, word_total, dopercentage );
  exit( EXIT_SUCCESS );
}