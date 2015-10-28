#include <set>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <climits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
typedef signed long int bitType;

void usage(){
  cerr << "gebruik: prog< anagramfile> <confusionfile> [<outputfile>]" << endl;
}


size_t split( const string& src, vector<string>& results,
	      const string& sep ){
  // split a string into substrings, using seps as seperator
  // silently skip empty entries (e.g. when two or more seperators co-incide)
  results.clear();
  string::size_type pos = 0, p;
  string res;
  while ( pos != string::npos ){
    p = src.find( sep, pos );
    if ( p == string::npos ){
      res = src.substr( pos );
      pos = p;
    }
    else {
      res = src.substr( pos, p - pos );
      pos = p + sep.length();
    }
    if ( !res.empty() )
      results.push_back( res );
  }
  return results.size();
}

template< typename T >
T stringTo( const std::string& str ) {
  T result;
  std::stringstream dummy ( str );
  if ( !( dummy >> result ) ) {
    throw( std::runtime_error( "conversion from string '"
			       + str + "' failed" ) );
  }
  return result;
}

int main( int argc, char **argv ){
  cerr << "maximaal aantal bits=" << numeric_limits<bitType>::max() << endl;
  if ( argc < 3 ){
    usage();
    exit(1);
  }
  string anaFile = argv[1];
  string confFile = argv[2];
  string outFile = "outfile";
  if ( argc == 4 ){
    outFile = argv[3];
  }

  ifstream ana( anaFile.c_str() );
  if ( !ana ){
    cerr << "problem opening " << anaFile << endl;
    exit(1);
  }
  cout << "reading hash values" << endl;
  set<bitType> anaSet;
  while ( ana ){
    bitType bit;
    ana >> bit;
    ana.ignore( INT_MAX, '\n' );
    anaSet.insert( bit );
  }
  cout << "read " << anaSet.size() << " anagram values" << endl;

  ifstream conf( confFile.c_str() );
  if ( !conf ){
    cerr << "problem opening confusable file: " << confFile << endl;
    exit(1);
  }

  cout << "reading confusion values" << endl;
  set<bitType> confSet;
  string line;
  size_t step = 0;
  while ( getline( conf, line ) ){
    if ( step % 50000 == 0 ){
      cerr << endl << step << ":" ;
    }
    if ( ++step % 1000 == 0 ){
      cerr << ".";
      cerr.flush();
    }
    vector<string> parts;
    if ( split( line, parts, "#" ) > 0 ){
      bitType bit = stringTo<bitType>( parts[0] );
      confSet.insert(bit);
    }
    else {
      cerr << "problems with line " << line << endl;
      cerr << "bail out " << endl;
      exit(1);
    }
  }
  cout << endl << "read " << confSet.size() << " confusion values" << endl;
  bitType max = *confSet.rbegin();
  cout <<"max value = " << max << endl;

  map<bitType,set<bitType> > result;
  set<bitType>::const_iterator it1 = anaSet.begin();
  size_t size = anaSet.size();
  step = 0;
  while ( it1 != anaSet.end() ){
    if ( step % 5000 == 0 ){
      cerr << step << "/" << size << ":" ;
    }
    if ( ++step % 100 == 0 ){
      cerr << ".";
      cerr.flush();
    }
    set<bitType>::const_iterator it2 = it1;
    ++it2;
    while ( it2 != anaSet.end() ){
      bitType diff = *it2 - *it1;
      if ( diff > max )
	break;
      set<bitType>::const_iterator sit = confSet.find( diff );
      if ( sit != confSet.end() ){
	result[diff].insert(*it1);
      }
      ++it2;
    }
    ++it1;
  }
  cerr << endl;
  ofstream of( outFile.c_str() );

  map<bitType,set<bitType> >::const_iterator rit = result.begin();
  while ( rit != result.end() ){
    of << rit->first << "#";
    set<bitType>::const_iterator it = rit->second.begin();
    while ( it != rit->second.end() ){
      of << *it;
      ++it;
      if ( it != rit->second.end() ){
	of << ",";
      }
    }
    of << endl;
    ++rit;
  }
}
