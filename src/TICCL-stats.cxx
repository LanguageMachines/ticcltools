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

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/XMLtools.h"
#include "ticcl/unicode.h"

#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace	std;
using namespace	TiCC;

bool verbose = false;

void create_wf_list( const map<string, unsigned int>& wc,
		     const string& filename, unsigned int totalIn,
		     unsigned int clip,
		     bool doperc ){
  unsigned int total = totalIn;
  ofstream os( filename );
  if ( !os ){
    cerr << "failed to create outputfile '" << filename << "'" << endl;
    exit(EXIT_FAILURE);
  }
  map<unsigned int, set<string> > wf;
  for( const auto& cit : wc ){
    if ( cit.second <= clip ){
      total -= cit.second;
    }
    else {
      wf[cit.second].insert( cit.first );
    }
  }
  unsigned int sum=0;
  unsigned int types=0;
  map<unsigned int, set<string> >::const_reverse_iterator wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for( const auto& sit : wit->second ){
      sum += wit->first;
      os << sit << "\t" << wit->first;
      if ( doperc ){
	os << "\t" << sum << "\t" << 100 * double(sum)/total;
      }
      os << endl;
      ++types;
    }
    ++wit;
  }
#pragma omp critical
  {
    cout << "created WordFreq list '" << filename << "'";
    if ( clip > 0 ){
      cout << endl << "with " << total << " words and " << types << " types. (" << totalIn - total
	   << " of the original " << totalIn << " words were clipped.)" << endl;
    }
    else {
      cout << " for " << total << " word tokens." << endl;
    }
  }
}

static int error_sink(void *mydata, xmlError *error ){
  int *cnt = (int*)mydata;
  if ( *cnt == 0 ){
    cerr << "\nXML-error: " << error->message << endl;
  }
  (*cnt)++;
  return 1;
}

bool is_emph( const string& data ){
  return (data.size() < 2) && isalnum(data[0]);
}

size_t tel( const xmlNode *node, bool lowercase,
	    size_t ngram, const string& sep,
	    map<string, unsigned int>& wc,
	    set<string>& emps ){
  vector<string> buffer(ngram);
  size_t cnt = 0;
  size_t buf_cnt = 0;
  bool in_emph = false;
  string emph_start;
  string emph_word;
  xmlNode *pnt = node->children;
  while ( pnt ){
    //    cerr << "bekijk label: " << (char*)pnt->name << endl;
    cnt += tel( pnt, lowercase, ngram, sep, wc, emps );
    if ( pnt->type == XML_TEXT_NODE ){
      string line  = (char*)( pnt->content );
      //      cerr << "text: " << line << endl;
      vector<string> v;
      TiCC::split( line, v );
      for ( const auto& word : v ){
	string wrd = word;
	if ( lowercase ){
	  UnicodeString us = UTF8ToUnicode( word );
	  us.toLower();
	  wrd = UnicodeToUTF8( us );
	}
	if ( is_emph( wrd ) ){
	  if ( in_emph ){
	    emph_word += "_" + wrd;
	  }
	  else {
	    emph_start = wrd;
	    in_emph = true;
	  }
	}
	else {
	  if ( in_emph && !emph_word.empty() ){
	    emps.insert( emph_start + emph_word );
	  }
	  in_emph = false;
	  emph_start.clear();
	  emph_word.clear();
	}
	buffer[buf_cnt++] = wrd;
	if ( buf_cnt == ngram ){
	  buf_cnt = ngram-1;
	  string gram;
	  for( size_t i=0; i < ngram -1; ++i ){
	    gram += buffer[i] + sep;
	    buffer[i] = buffer[i+1];
	  }
	  gram+= buffer[ngram-1];
#pragma omp critical
	  {
	    ++wc[gram];
	  }
	  ++cnt;
	}
      }
    }
    pnt = pnt->next;
  }
  return cnt;
}

size_t word_xml_inventory( const string& docName,
			   bool lowercase,
			   size_t ngram,
			   const string& sep,
			   map<string,unsigned int>& wc,
			   set<string>& emps ){
  xmlDoc *d = 0;
  int cnt = 0;
  xmlSetStructuredErrorFunc( &cnt, (xmlStructuredErrorFunc)error_sink );
  d = xmlReadFile( docName.c_str(), 0, XML_PARSE_NOBLANKS|XML_PARSE_HUGE );
  if ( !d || cnt > 0 ){
#pragma omp critical
    {
      cerr << "failed to load document '" << docName << "'" << endl;
    }
    return 0;
  }
  xmlNode *root = xmlDocGetRootElement( d );
  size_t wordTotal = tel( root, lowercase, ngram, sep, wc, emps );
  xmlFree( d );
  return wordTotal;
}

size_t word_inventory( const string& docName,
		       bool lowercase,
		       size_t ngram,
		       const string& sep,
		       map<string,unsigned int>& wc,
		       set<string>& emps ){
  vector<string> buffer(ngram);
  size_t buf_cnt = 0;
  size_t wordTotal = 0;
  ifstream is( docName );
  string line;
  bool in_emph = false;
  string emph_start;
  string emph_word;
  while ( getline( is, line ) ){
    vector<string> v;
    TiCC::split( line, v );
    for ( const auto& word : v ){
      string wrd = word;
      if ( lowercase ){
	UnicodeString us = UTF8ToUnicode( word );
	us.toLower();
	wrd = UnicodeToUTF8( us );
      }
      if ( is_emph( wrd ) ){
	if ( in_emph ){
	  emph_word += "_" + wrd;
	}
	else {
	  emph_start = wrd;
	  in_emph = true;
	}
      }
      else {
	if ( in_emph && !emph_word.empty() ){
	  emps.insert( emph_start + emph_word );
	}
	in_emph = false;
	emph_start.clear();
	emph_word.clear();
      }
      buffer[buf_cnt++] = wrd;
      if ( buf_cnt == ngram ){
	buf_cnt = ngram-1;
	string gram;
	for( size_t i=0; i < ngram -1; ++i ){
	  gram += buffer[i] + sep;
	  buffer[i] = buffer[i+1];
	}
	gram += buffer[ngram-1];
#pragma omp critical
	{
	  ++wc[gram];
	}
	++wordTotal;
      }
    }
  }
  return wordTotal;
}


void usage( const string& name ){
  cerr << "Usage: " << name << " [options] file/dir" << endl;
  cerr << "\t TICCL-stats will produce ngram statistics for a file, " << endl;
  cerr << "\t or a whole directory of files " << endl;
  cerr << "\t The output will be a 2 columned tab separated file, extension: *tsv " << endl;
  cerr << "\t--clip\t clipping factor. " << endl;
  cerr << "\t\t\t(entries with frequency <= this factor will be ignored). " << endl;
  cerr << "\t-p\t output percentages too. " << endl;
  cerr << "\t--lower\t Lowercase all words" << endl;
  cerr << "\t--ngram\t create an ngram list (default 1-gram)" << endl;
  cerr << "\t--underscore\t separate the ngram entries with a '_' (default space)" << endl;
  cerr << "\t--hemp=<file>. Create a historical emphasis file. " << endl;
  cerr << "\t-t\t number_of_threads" << endl;
  cerr << "\t-h\t this message" << endl;
  cerr << "\t-v\t very verbose output." << endl;
  cerr << "\t-V\t show version " << endl;
  cerr << "\t-e\t expr: specify the expression all input files should match with." << endl;
  cerr << "\t-o\t name of the output file(s) prefix." << endl;
  cerr << "\t-X\t the inputfiles are assumed to be XML. (all TEXT nodes are used)" << endl;
  cerr << "\t-R\t search the dirs recursively (when appropriate)." << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hVvpe:t:o:RX", "clip:,lower,ngram:,underscore,hemp:" );
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
  int clip = 0;
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
  bool doXML = opts.extract( 'X' );
  bool dopercentage = opts.extract('p');
  bool lowercase = opts.extract("lower");
  bool recursiveDirs = opts.extract( 'R' );
  bool do_under = opts.extract( "underscore" );
  string hempName;
  opts.extract("hemp", hempName );
  if ( !opts.extract( 'o', outputPrefix ) ){
    cerr << "an output filename prefix is required. (-o option) " << endl;
    exit(EXIT_FAILURE);
  }
  string value;
  if ( opts.extract("clip", value ) ){
    if ( !stringTo(value, clip ) ){
      cerr << "illegal value for --clip (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
  size_t ngram = 1;
  if ( opts.extract("ngram", value ) ){
    if ( !stringTo(value, ngram )
	 || (ngram < 1) ){
      cerr << "illegal value for --ngram (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
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
  string sep = " ";
  if ( do_under ){
    sep = "_";
  }

  set<string> hemp;
#pragma omp parallel for shared(fileNames,wordTotal,wc,hemp)
  for ( size_t fn=0; fn < fileNames.size(); ++fn ){
    string docName = fileNames[fn];
    unsigned int word_count =  0;
    if ( doXML ){
      word_count = word_xml_inventory( docName, lowercase, ngram, sep, wc, hemp );
    }
    else {
      word_count = word_inventory( docName, lowercase, ngram, sep, wc, hemp );
    }
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
  cout << "start calculating the results" << endl;
  if ( !hempName.empty() ){
    ofstream out( hempName );
    if ( out ){
      for( auto const& it : hemp ){
	out << it << endl;
      }
      cout << "historical emphasis stored in: " << hempName << endl;
    }
    else {
      cerr << "unable to create historical emphasis file: " << hempName << endl;
    }
  }
  string filename = outputPrefix + ".wordfreqlist";
  string ng = toString(ngram);
  filename += "." + ng + ".tsv";
  create_wf_list( wc, filename, wordTotal, clip, dopercentage );
  exit( EXIT_SUCCESS );
}
