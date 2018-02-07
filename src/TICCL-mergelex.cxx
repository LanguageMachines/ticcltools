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
		     const string& filename, unsigned int totalIn ){
  unsigned int total = totalIn;
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
      os << sit << "\t" << wit->first << endl;
      ++types;
    }
    ++wit;
  }
#pragma omp critical
  {
    cout << "created WordFreq list '" << filename << "'" << endl
	 << "with " << total << " tokens and " << types
	 << " types. TTR= " << (double)types/total
	 << ", the angle is " << atan((double)types/total)*180/M_PI
	 << " degrees" << endl;
  }
}

size_t read_words( const string& docName, map<string,unsigned int>& wc ){
  size_t wordTotal = 0;
  ifstream is( docName );
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
    wordTotal += frq;
  }
  return wordTotal;
}


void usage( const string& name ){
  cerr << "Usage: " << name << " [options] file/dir" << endl;
  cerr << "\t TICCL-mergelex will create a merged lexicond from range of" << endl;
  cerr << "\t lexicon files in TICCL-lexstat or FoLiA-stats format." << endl;
  cerr << "\t The output will be a 2 columned tab separated file, extension: *tsv " << endl;
  cerr << "\t-t\t number_of_threads" << endl;
  cerr << "\t-h\t this message" << endl;
  cerr << "\t-v\t very verbose output." << endl;
  cerr << "\t-V\t show version " << endl;
  cerr << "\t-e\t expr: specify the expression all input files should match with." << endl;
  cerr << "\t-o\t name of the output file(s) prefix." << endl;
  cerr << "\t-R\t search the dirs recursively (when appropriate)." << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hVve:t:o:R", "" );
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
  int numThreads = 1;
#endif
  string expression;
  string outputPrefix;
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
  if ( !opts.extract( 'o', outputPrefix ) ){
    cerr << "an output filename prefix is required. (-o option) " << endl;
    exit(EXIT_FAILURE);
  }
  string value;
  if ( opts.extract('t', value ) ){
#ifdef HAVE_OPENMP
    if ( !stringTo(value, numThreads ) ){
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
#else
    cerr << "OpenMP support is missing. -t option is not supported" << endl;
    exit( EXIT_FAILURE );
#endif
  }
  opts.extract('e', expression );
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_OPENMP
  if ( numThreads != 1 )
    omp_set_num_threads( numThreads );
#endif

  vector<string> massOpts = opts.getMassOpts();
  if ( massOpts.empty() ){
    cerr << "no file or dir specified!" << endl;
    exit(EXIT_FAILURE);
  }
  string name = massOpts[0];
  vector<string> fileNames = searchFilesMatch( name, expression, recursiveDirs );
  size_t toDo = fileNames.size();
  if ( toDo == 0 ){
    cerr << "no matching files found" << endl;
    exit(EXIT_SUCCESS);
  }

  string::size_type pos = outputPrefix.find( "." );
  if ( pos != string::npos && pos == outputPrefix.length()-1 ){
    // outputname ends with a .
    outputPrefix = outputPrefix.substr(0,pos);
  }
  pos = outputPrefix.find( "/" );
  if ( pos != string::npos && pos == outputPrefix.length()-1 ){
    // outputname ends with a /
    outputPrefix += "ticclstats";
  }

  if ( toDo > 1 ){
    cout << "start processing of " << toDo << " files " << endl;
  }
  map<string,unsigned int> wc;
  unsigned int wordTotal =0;
#pragma omp parallel for shared(fileNames,wordTotal,wc)
  for ( size_t fn=0; fn < fileNames.size(); ++fn ){
    string docName = fileNames[fn];
    unsigned int word_count = read_words( docName, wc );
    wordTotal += word_count;
#pragma omp critical
    {
      cout << "Processed :" << docName << " with " << word_count << " words,"
	   << " still " << --toDo << " files to go." << endl;
    }
  }
  if ( toDo > 1 ){
    cout << "done processsing directory '" << name << "' in total "
	 << wordTotal << " words were found." << endl;
  }
  cout << "start outputting the results" << endl;
  string filename = outputPrefix + ".wordfreqlist.tsv";
  create_wf_list( wc, filename, wordTotal );
  exit( EXIT_SUCCESS );
}
