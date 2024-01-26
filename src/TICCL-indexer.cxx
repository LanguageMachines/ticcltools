/*
  Copyright (c) 2006 - 2024
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
#include <cassert>
#include <unistd.h>
#include <set>
#include <map>
#include <limits>
#include <algorithm>
#include <vector>
#include <climits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/Unicode.h"
#include "ticcutils/FileUtils.h"
#include "ticcl/ticcl_common.h"

#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace std;
using namespace icu;

using ticcl::bitType;

set<bitType> follow_nums;

void usage( const string& name ){
  cerr << name << endl;
  cerr << "options: " << endl;
  cerr << "\t--hash=<anahash>\tname of the anagram hashfile. (produced by TICCL-anahash)" << endl;
  cerr << "\t--charconf=<charconf>\tname of the character confusion file. (produced by TICCL-lexstat)" << endl;
  cerr << "\t--foci=<focifile>\tname of the file produced by the --artifrq parameter of TICCL-anahash." << endl;
  cerr << "\t\t\t This file is used to limit the searchspace" << endl;
  cerr << "\t-o <outputfile>\t\tname for the outputfile. " << endl;
  cerr << "\t--confstats=<statsfile>\tcreate a list of confusion statistics"
       << endl;
  cerr << "\t--low=<low>\t skip entries from the anagram file shorter than "
       << endl;
  cerr << "\t\t\t'low' characters. (default = 5)" << endl;
  cerr << "\t--high=<high>\t skip entries from the anagram file longer than "
       << endl;
  cerr << "\t\t\t'high' characters. (default=35)" << endl;
  cerr << "\t-t <threads> or --threads <threads>\n\t\t\t Number of threads to run on." << endl;
  cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  cerr << "\t\t\t reasonable value. ($OMP_NUM_TREADS - 2)" << endl;
  cerr << "\t-v\t\t run verbose " << endl;
  cerr << "\t-V or --version\t show version " << endl;
  cerr << "\t-h or --help\t this message " << endl;
}

struct experiment {
  set<bitType>::const_iterator start;
  set<bitType>::const_iterator finish;
};


void handle_confs( const experiment& exp,
		   size_t& count,
		   const set<bitType>& anaSet,
		   const set<bitType>& focSet,
		   ostream &of,
		   ostream *csf ){
  bitType vorige = 0;
  bitType totalShift = 0;
  auto sit = exp.start;
  while ( sit != exp.finish ){
    set<bitType> result;
#pragma omp critical(count)
    {
      if ( ++count % 100 == 0 ){
	cout << ".";
	cout.flush();
	if ( count % 5000 == 0 ){
	  cout << endl << count << endl;
	}
      }
    }
    bitType confusie = *sit;
    if ( follow_nums.find(confusie) != follow_nums.end() ){
      cerr << "found confusion value: " << confusie << endl;
    }
    bitType diff = confusie - vorige;
    if ( follow_nums.find(diff) != follow_nums.end() ){
      cerr << "found a difference value: " << diff << endl;
    }
    totalShift += diff;
    auto it1 = anaSet.begin();
    auto it2 = it1;
    while ( it1 != anaSet.end() && it2 != anaSet.end() ){
      bitType v1 = *it1;
      bitType v2 = *it2;
      bitType v2_save = v2;
      if ( v2 >= totalShift ) {
	v2 -= totalShift;
      }
      else {
	v2 = 0;
      }
      if ( v1 == v2 ){
	// if ( follow_nums.find(v1) != follow_nums.end() ){
	//   cerr << "found a possible focus value: " << v1 << endl;
	// }
	bool foc = true;
	if ( !focSet.empty() ){
	  // do we have to focus?
	  foc = !( focSet.find( v1 ) == focSet.end()
		   && focSet.find( v2_save ) == focSet.end() );
	  // not if both values out of focus
	}
	if ( foc ){
	  auto emp = result.emplace(v1);
	  if ( emp.second
	       && follow_nums.find(v1) != follow_nums.end() ){
	    cerr << "stored a focus value: " << v1 << endl;
	  }
	}
	++it1;
	++it2;
      }
      else if ( v1 < v2 ){
	++it1;
      }
      else {
	++it2;
      }
    }
    vorige = confusie;
    ++sit;
    if ( !result.empty() ){
      stringstream ss;
      ss << confusie << "#";
      bool hit = false;
      for ( const auto& it : result ){
	if ( it != *result.begin() ){
	  ss << ",";
	}
	if ( follow_nums.find(it) != follow_nums.end() ){
	  cerr << "Store " << it << " for confusion: " << confusie
	       << endl;
	  hit = true;
	}
	ss << it;
      }
      if ( hit
	   || follow_nums.find(confusie) != follow_nums.end()){
	cerr << "Stored followed value(s) in: " << ss.str() << endl;
      }
#pragma omp critical(update)
      {
	of << ss.str() << endl;
	if ( csf ){
	  *csf << confusie << "#" << result.size() << endl;
	}
      }
    }
  }
}

size_t init( vector<experiment>& exps,
	     const set<bitType>& hashes,
	     size_t threads ){
  exps.clear();
  size_t partsize = hashes.size() / threads;
  if ( partsize < 1 ){
    experiment e;
    e.start = hashes.begin();
    e.finish = hashes.end();
    exps.push_back( e );
    return 1;
  }
  auto hash_it = hashes.begin();
  for ( size_t i=0; i < threads; ++i ){
    experiment e;
    e.start = hash_it;
    for ( size_t j=0; j < partsize && hash_it != hashes.end(); ++j ){
      ++hash_it;
    }
    e.finish = hash_it;
    exps.push_back( e );
  }
  if ( hash_it != hashes.end() ){
    exps[exps.size()-1].finish = hashes.end();
  }
  return threads;
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "charconf:,hash:,low:,high:,help,version,"
			   "foci:,threads:,confstats:,follow:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('h') || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V') || opts.extract("version") ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( argc < 3	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  bool verbose = opts.extract( 'v' );
  string anahashFile;
  string confFile;
  string fociFile;
  string outFile;
  string confstats_file;
  int lowValue = 5;
  int highValue = 35;
  opts.extract( "hash", anahashFile );
  opts.extract( "charconf", confFile );
  opts.extract( "confstats", confstats_file );
  opts.extract( "foci", fociFile );
  opts.extract( 'o', outFile );
  string value;
  while ( opts.extract( "follow", value ) ){
    bitType follow_num;
    if ( !TiCC::stringTo(value,follow_num) ) {
      cerr << "illegal value for --follow (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    follow_nums.insert( follow_num );
  }
  if ( opts.extract("low", value ) ){
    if ( !TiCC::stringTo(value,lowValue) ) {
      cerr << "illegal value for --low (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract("high", value ) ){
    if ( !TiCC::stringTo(value,highValue) ) {
      cerr << "illegal value for --high (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  int numThreads=1;
  value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
#ifdef HAVE_OPENMP
  if ( TiCC::lowercase(value) == "max" ){
    numThreads = omp_get_max_threads() - 2;
  }
  else {
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
#else
  if ( value != "1" ){
    cerr << "unable to set number of threads!.\nNo OpenMP support available!"
	 <<endl;
    exit(EXIT_FAILURE);
  }
#endif
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  if ( !TiCC::isFile(anahashFile) ){
    cerr << "problem opening corpus anagram hashfile: " << anahashFile << endl;
    exit(1);
  }
  if ( !TiCC::isFile(confFile) ){
    cerr << "problem opening charconfusion file: " << confFile << endl;
    exit(1);
  }

  set<bitType> focSet;
  if ( !fociFile.empty() ){
    ifstream foc( fociFile );
    if ( !foc ){
      cerr << "problem opening foci file: " << fociFile << endl;
      exit(1);
    }
    focSet = ticcl::read_bit_set( foc );
    cout << "read " << focSet.size() << " foci values" << endl;
  }

  if ( outFile.empty() ){
    outFile = anahashFile;
    string::size_type pos = outFile.rfind(".");
    if ( pos != string::npos ){
      outFile.resize(0,pos);
    }
    outFile += ".index";
  }
  else if ( !TiCC::match_back( outFile, ".index" ) ){
    outFile += ".index";
  }

  ofstream of( outFile );
  if ( !of ){
    cerr << "problem opening outputfile: " << outFile << endl;
    exit(1);
  }
  ofstream *csf = 0;
  if ( !confstats_file.empty() ){
    csf = new ofstream( confstats_file );
    if ( !csf ){
      cerr << "problem opening outputfile: " << confstats_file << endl;
      exit(1);
    }
  }
  cout << "reading corpus word anagram hash values" << endl;
  ifstream ana( anahashFile );
  size_t skipped = 0;
  set<bitType> anaSet = ticcl::read_anahash( ana,
					     lowValue,
					     highValue,
					     skipped,
					     verbose );
  cout << "read " << anaSet.size() << " corpus anagram values" << endl;
  cout << "skipped " << skipped << " out-of-band corpus anagram values" << endl;

  cout << "reading character confusion anagram values" << endl;
  ifstream conf( confFile );
  set<bitType> confSet = ticcl::read_confusions( conf );
  cout << endl << "read " << confSet.size()
       << " character confusion anagram values" << endl;

  vector<experiment> experiments;
  size_t expsize = init( experiments, confSet, numThreads );
#ifdef HAVE_OPENMP
  omp_set_num_threads( expsize );
  cout << "running on " << expsize << " threads." << endl;
#endif


  cout << "processing all character confusion values" << endl;
  size_t count = 0;
#pragma omp parallel for shared( experiments, of, csf )
  for ( size_t i=0; i < expsize; ++i ){
    handle_confs( experiments[i], count, anaSet, focSet, of, csf );
  }
  cout << "\nwrote indexes into: " << outFile << endl;
  if ( csf ){
    cout << "wrote confusion statistics into: " << confstats_file << endl;
    csf->close();
    delete csf;
  }
}
