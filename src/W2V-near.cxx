/*
  Copyright (c) 2019 - 2023
  CLST  - Radboud University

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
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  For questions and suggestions, see:
      https://github.com/LanguageMachines/ticcltools/issues
  or send mail to:
      lamasoftware (at ) science.ru.nl

*/
#include <iostream>
#include <fstream>
#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/Unicode.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace TiCC;

void usage( const string& name ){
  cerr << name << " --vectors=vectorfile [-n size] [FILES]" << endl;
}

int main( int argc, char *argv[] ){
  CL_Options opts( "hn:", "vectors:" );
  try {
    opts.init(argc,argv);
  }
  catch( OptionError& e ){
    cerr << e.what() << endl;
    usage( opts.prog_name() );
    exit( EXIT_FAILURE );
  }
  if ( opts.extract( 'h' ) ){
    usage( opts.prog_name() );
    exit( EXIT_SUCCESS );
  }
  string vectorsFile;
  if ( !opts.extract( "vectors", vectorsFile ) ){
    cerr << "missing '--vectors' option" << endl;
    exit( EXIT_FAILURE );
  }
  int NN = 20;
  string value;
  if ( opts.extract( 'n', value ) ){
    NN = stringTo<int>(value);
  }
  auto fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing input file(s)" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.empty() ) {
    cerr << "unsupported options: " << opts.toString() << endl;
    exit( EXIT_FAILURE );
  }
  wordvec_tester WV;
  if ( !WV.fill( vectorsFile ) ){
    cerr << "fill failed from " << vectorsFile << endl;
    exit(EXIT_FAILURE);
  }
  else
    cerr << "filled with " << WV.size() << " vectors" << endl;
  for ( auto const& name : fileNames ){
    ifstream is( name );
    if ( !is ){
      cerr << "failed to read: " << name << endl;
      continue;
    }
    string outname = name + ".out";
    ofstream os( outname );
    if ( !os ){
      cerr << "failed to open: " << outname << endl;
      continue;
    }
    UnicodeString line;
    while ( TiCC::getline( is, line ) ){
      os << "NEIGHBORS of '" << line << "':" << endl;
      vector<word_dist> res;
      if ( WV.lookup( TiCC::UnicodeToUTF8(line), NN, res ) ){
	for ( auto const& i : res ){
	  os << "\t" << i.w << "\t" << i.d << endl;
	}
      }
      else {
	os << "\tNone" << endl;
      }
    }
    cerr << "results in: " << outname << endl;
  }
  return EXIT_SUCCESS;
}
