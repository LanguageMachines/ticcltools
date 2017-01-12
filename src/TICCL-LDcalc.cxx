/*
  Copyright (c) 2006 - 2017
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

#include <unistd.h>
#include <set>
#include <map>
#include <limits>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcl/unicode.h"
#include "config.h"

using namespace std;
typedef signed long int bitType;

bool verbose = false;

void usage( const string& progname ){
  cerr << "usage: " << progname << endl;
  cerr << "\t--index <confuslist> as produced by TICCL-indexer or TICCL-indexerNT." << endl;
  cerr << "\t--hash <anahash>, as produced by TICCl-anahash," << endl;
  cerr << "\t--clean <cleanfile> as produced by TICCL-unk" << endl;
  cerr << "\t--diac <diacritics file> a list of 'diacritical' confusions." << endl;
  cerr << "\t--hist <historicalfile> a list of 'historical' confusions." << endl;
  cerr << "\t--alph <alphabet> an alphabet file (as produced by TICCL-lexstat)" << endl;
  cerr << "\t--nohld ignore --LD for 'historical' confusions." << endl;
  cerr << "\t-o <outputfile>" << endl;
  cerr << "\t-t <threads> Number of threads to run on." << endl;
  cerr << "\t--LD <distance> The Levensthein (or edit) distance to use" << endl;
  cerr << "\t--artifrq <artifreq> " << endl;
  cerr << "\t-h or --help this message " << endl;
  cerr << "\t-v verbose " << endl;
  cerr << "\t-V or --version show version " << endl;
  exit( EXIT_FAILURE );
}

bitType high_five( int val ){
  bitType result = val;
  result *= val;
  result *= val;
  result *= val;
  result *= val;
  return result;
}

unsigned int ldCompare( const UnicodeString& s1, const UnicodeString& s2 ){
  const size_t len1 = s1.length(), len2 = s2.length();
  vector<unsigned int> col(len2+1), prevCol(len2+1);
  for ( unsigned int i = 0; i < prevCol.size(); ++i ){
    prevCol[i] = i;
  }
  for ( unsigned int i = 0; i < len1; ++i ) {
    col[0] = i+1;
    for ( unsigned int j = 0; j < len2; ++j )
      col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]),
		      prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
    col.swap(prevCol);
  }
  unsigned int result = prevCol[len2];
  return result;
}

bool isClean( const UnicodeString& us, const set<UChar>& alfabet ){
  if ( alfabet.empty() )
    return true;
  for ( int i=0; i < us.length(); ++i ){
    if ( alfabet.find( us[i] ) == alfabet.end() )
      return false;
  }
  return true;
}

void handleTranspositions( ostream& os, const set<string>& s,
			   const map<string,size_t>& freqMap,
			   const map<UnicodeString,size_t>& low_freqMap,
			   const set<UChar>& alfabet,
			   size_t freqTreshold,
			   bool isKHC,
			   bool noKHCld,
			   bool isDIAC ){
  set<string>::const_iterator it1 = s.begin();
  while ( it1 != s.end() ) {
    string str1 = *it1;
    map<string,size_t>::const_iterator fit = freqMap.find( str1 );
    if ( fit == freqMap.end() ){
      if ( verbose ){
#pragma omp critical
	{
	  cout << "not found in freq file " << str1 << endl;
	}
      }
      ++it1;
      continue;
    }
    size_t freq1 = fit->second;
    set<string>::const_iterator it2 = it1;
    ++it2;
    while ( it2 != s.end() ) {
      string str2 = *it2;
      map<string,size_t>::const_iterator fit = freqMap.find( str2 );
      if ( fit == freqMap.end() ){
	if ( verbose ){
#pragma omp critical
	  {
	    cout << "not found in freq file " << str2 << endl;
	  }
	}
	++it2;
	continue;
      }
      size_t freq2 = fit->second;
      UnicodeString us1 = UTF8ToUnicode( str1 );
      us1.toLower();
      UnicodeString us2 = UTF8ToUnicode( str2 );
      us2.toLower();

      size_t out_freq1;
      size_t out_low_freq1;
      size_t out_freq2;
      size_t out_low_freq2;
      string out_str1;
      string out_str2;
      size_t low_freq1 = low_freqMap.at(us1);
      size_t low_freq2 = low_freqMap.at(us2);
      if ( low_freq1 >= freqTreshold && low_freq2 >= freqTreshold
	   && !isDIAC ){
	++it2;
	continue;
      }
      if ( low_freq1 >= low_freq2 ){
	if ( low_freq1 < freqTreshold ){
	  ++it2;
	  continue;
	}
      }
      else {
	if ( low_freq2 < freqTreshold ){
	  ++it2;
	  continue;
	}
      }

      size_t canon_freq = 0;
      UnicodeString candidate;
      if ( low_freq1 > low_freq2 ){
	canon_freq = low_freq1;
	out_freq1 = freq2;
	out_low_freq1 = low_freq2;
	out_freq2 = freq1;
	out_low_freq2 = low_freq1;
	out_str1 = str2;
	out_str2 = str1;
	candidate = us1;
      }
      else {
	canon_freq = low_freq2;
	out_freq1 = freq1;
	out_low_freq1 = low_freq1;
	out_freq2 = freq2;
	out_low_freq2 = low_freq2;
	out_str1 = str1;
	out_str2 = str2;
	candidate = us2;
      }
      if ( !isClean( candidate, alfabet ) ){
	if ( verbose ){
#pragma omp critical
	  {
	    cout << "ignore dirty candidate " << candidate << endl;
	  }
	}
	++it2;
	continue;
      }
      unsigned int ld = ldCompare( us1, us2 );
      if ( ld != 2 ){
	if ( !( isKHC && noKHCld ) ){
	  if ( verbose ){
#pragma omp critical
	    {
	      cout << " LD != 2 " << str1 << "," << str2 << endl;
	    }
	  }
	  ++it2;
	  continue;
	}
      }

      int cls = max(us1.length(),us2.length()) - ld;
      string canon = "0";
      if ( canon_freq >= freqTreshold ){
	canon = "1";
      }
      string FLoverlap = "0";
      if ( us1[0] == us2[0] ){
	FLoverlap = "1";
      }
      string LLoverlap = "0";
      if ( us1.length() > 1 && us2.length() > 1
	   && us1[us1.length()-1] == us2[us2.length()-1]
	   && us1[us1.length()-2] == us2[us2.length()-2] ){
	LLoverlap = "1";
      }
      string KHC = "0";
      if ( isKHC ){
	KHC = "1";
      }
      string result = out_str1 + "~" + TiCC::toString(out_freq1) + "~"
	+ TiCC::toString(out_low_freq1) + "~"
	+ out_str2 + "~" + TiCC::toString( out_freq2 ) + "~"
	+ TiCC::toString(out_low_freq2) + "~"
	+ "~0~" + TiCC::toString( ld ) + "~"
	+ TiCC::toString(cls) + "~" + canon + "~"
	+ FLoverlap + "~" + LLoverlap + "~"
	+ KHC;
#pragma omp critical
      {
	os << result << endl;
      }
      ++it2;
    }
    ++it1;
  }
}

void compareSets( ostream& os, unsigned int ldValue,
		  const string& KWC,
		  const set<string>& s1, const set<string>& s2,
		  const map<string,size_t>& freqMap,
		  const map<UnicodeString,size_t>& low_freqMap,
		  set<UChar>& alfabet,
		  size_t freqTreshold,
		  bool isKHC,
		  bool noKHCld,
		  bool isDIAC ){
  // using TiCC::operator<<;
  // cerr << "set 1 " << s1 << endl;
  // cerr << "set 2 " << s2 << endl;
  set<string>::const_iterator it1 = s1.begin();
  while ( it1 != s1.end() ) {
    string str1 = *it1;
    if ( verbose ){
#pragma omp critical
      {
	cout << "string 1 " << str1 << endl;
      }
    }
    map<string,size_t>::const_iterator fit = freqMap.find( str1 );
    if ( fit == freqMap.end() ){
      if ( verbose ){
#pragma omp critical
	{
	  cout << "not found in freq file " << str1 << endl;
	}
      }
      ++it1;
      continue;
    }
    size_t freq1 = fit->second;
    UnicodeString us1 = UTF8ToUnicode( str1 );
    us1.toLower();
    set<string>::const_iterator it2 = s2.begin();
    while ( it2 != s2.end() ) {
      string str2 = *it2;
      if ( verbose ){
#pragma omp critical
	{
	  cout << "string 2 " << str2 << endl;
	}
      }
      fit = freqMap.find( str2 );
      if ( fit == freqMap.end() ){
	if ( verbose ){
#pragma omp critical
	  {
	    cout << "not found in freq file " << str2 << endl;
	  }
	}
	++it2;
	continue;
      }

      size_t freq2 = fit->second;
      UnicodeString us2 = UTF8ToUnicode( str2 );
      us2.toLower();
      unsigned int ld = ldCompare( us1, us2 );
      if ( ld > ldValue ){
	if ( !( isKHC && noKHCld ) ){
	  if ( verbose ){
#pragma omp critical
	    {
	      cout << " LD too high " << str1 << "," << str2 << endl;
	    }
	  }
	  ++it2;
	  continue;
	}
      }

      size_t out_freq1;
      size_t out_low_freq1;
      size_t out_freq2;
      size_t out_low_freq2;
      string out_str1;
      string out_str2;
      size_t low_freq1 = low_freqMap.at(us1);
      size_t low_freq2 = low_freqMap.at(us2);
      size_t canon_freq = 0;
      UnicodeString candidate;
      if ( low_freq1 > low_freq2 ){
	canon_freq = low_freq1;
	out_freq1 = freq2;
	out_low_freq1 = low_freq2;
	out_freq2 = freq1;
	out_low_freq2 = low_freq1;
	out_str1 = str2;
	out_str2 = str1;
	candidate = us1;
      }
      else {
	canon_freq = low_freq2;
	out_freq1 = freq1;
	out_low_freq1 = low_freq1;
	out_freq2 = freq2;
	out_low_freq2 = low_freq2;
	out_str1 = str1;
	out_str2 = str2;
	candidate = us2;
      }
      if ( !isClean( candidate, alfabet ) ){
	if ( verbose ){
#pragma omp critical
	  {
	    cout << "ignore dirty candidate " << candidate << endl;
	  }
	}
	++it2;
	continue;
      }

      if ( out_low_freq1 >= freqTreshold && !isDIAC ){
	if ( verbose ){
#pragma omp critical
	  {
	    cout << "lexical word " << out_str1 << endl;
	  }
	}
	++it2;
	continue;
      }

      int cls = max(us1.length(),us2.length()) - ld;
      string canon = "0";
      if ( canon_freq >= freqTreshold ){
	canon = "1";
      }
      string FLoverlap = "0";
      if ( us1[0] == us2[0] ){
	FLoverlap = "1";
      }
      string LLoverlap = "0";
      if ( us1.length() > 1 && us2.length() > 1
	   && us1[us1.length()-1] == us2[us2.length()-1]
	   && us1[us1.length()-2] == us2[us2.length()-2] ){
	LLoverlap = "1";
      }
      string KHC = "0";
      if ( isKHC ){
	KHC = "1";
      }
      string result = out_str1 + "~" + TiCC::toString(out_freq1) +
	+ "~" + TiCC::toString(out_low_freq1) +
	+ "~" + out_str2 + "~" + TiCC::toString( out_freq2 )
	+ "~" + TiCC::toString( out_low_freq2 )
	+ "~" + KWC + "~" + TiCC::toString( ld ) + "~"
	+ TiCC::toString(cls) + "~" + canon + "~"
	+ FLoverlap + "~" + LLoverlap + "~"
	+ KHC;
#pragma omp critical
      {
	os << result << endl;
      }
      ++it2;
    }
    ++it1;
  }
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "diac:,hist:,nohld,artifrq:,LD:,hash:,clean:,"
			   "alph:,index:,help,version" );
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
  if ( opts.extract('h') || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V') || opts.extract("version") ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  bool verbose = opts.extract( 'v' );

  string indexFile;
  string anahashFile;
  string frequencyFile;
  string histconfFile;
  string diaconfFile;
  string alfabetFile;
  int numThreads=1;
  int LDvalue=2;
  bool backward = false;
  bool noKHCld = opts.extract("nohld");
  if ( !opts.extract( "index", indexFile ) ){
    cerr << "missing --index option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( TiCC::match_back( indexFile, ".index" ) ){
    backward = false;
  }
  else if ( TiCC::match_back( indexFile, ".indexNT" ) ){
    backward = true;
  }
  else {
    cerr << "--index files must have extension: '.index' or '.indexNT' "
	 << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "hash", anahashFile ) ){
    cerr << "missing --hash option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "clean", frequencyFile ) ){
    cerr << "missing --clean option" << endl;
    exit( EXIT_FAILURE );
  }
  opts.extract( "alph", alfabetFile );
  opts.extract( "hist", histconfFile );
  if ( opts.extract( "diac", diaconfFile ) ){
    if ( !TiCC::match_back( diaconfFile, ".diac" ) ){
      cerr << "invalid extension for --diac file '" << diaconfFile
	   << "' (must be .diac) " << endl;
      exit(EXIT_FAILURE);
    }
  }
  string outFile;
  if ( opts.extract( 'o', outFile ) ){
    if ( !TiCC::match_back( outFile, ".ldcalc" ) )
      outFile += ".ldcalc";
  }
  else {
    outFile = indexFile + ".ldcalc";
  }
  size_t artifreq = 0;
  string value;

  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract( 't', value ) ){
#ifdef HAVE_OPENMP
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
#else
    cerr << "You don't have OpenMP supprt. The -t option is useless" << endl;
    exit( EXIT_FAILURE );
#endif
  }
  if ( opts.extract( "LD", value ) ){
    if ( !TiCC::stringTo(value,LDvalue) ) {
      cerr << "illegal value for --LD (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    if ( LDvalue < 1 || LDvalue > 10 ){
      cerr << "invalid LD value: " << LDvalue << " (1-10 is OK)" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

  set<UChar> alfabet;
  if ( !alfabetFile.empty() ){
    ifstream lexicon( alfabetFile );
    if ( !lexicon ){
      cerr << "problem opening alfabet file: " << alfabetFile << endl;
      exit(1);
    }
    cout << "reading alphabet: " << alfabetFile << endl;
    string line;
    while ( getline( lexicon, line ) ){
      if ( line.size() == 0 || line[0] == '#' )
	continue;
      vector<string> vec;
      if ( TiCC::split( line, vec ) != 3 ){
	cerr << "invalid line '" << line << "' in " << alfabetFile << endl;
	exit( EXIT_FAILURE );
      }
      UnicodeString key = UTF8ToUnicode(vec[0]);
      alfabet.insert(key[0]);
    }
  }
  cout << "read " << alfabet.size() << " letters with frequencies" << endl;

  ifstream ff( frequencyFile  );
  if ( !ff ){
    cerr << "problem opening " << frequencyFile << endl;
    exit(1);
  }
  cout << "reading clean file: " << frequencyFile << endl;
  map<string, size_t> freqMap;
  map<UnicodeString, size_t> low_freqMap;
  string line;
  size_t ign = 0;
  while ( getline( ff, line ) ){
    vector<string> v1;
    if ( TiCC::split( line, v1 ) != 2 ){
      ++ign;
      continue;
    }
    else {
      string s = v1[0];
      size_t freq = TiCC::stringTo<size_t>( v1[1] );
      freqMap[s] = freq;
      UnicodeString us = UTF8ToUnicode( s );
      us.toLower();
      low_freqMap[us] +=freq;
    }
  }
  cout << "read " << freqMap.size() << " clean words with frequencies" << endl;
  cout << "skipped " << ign << " n-grams" << endl;

  set<bitType> histMap;
  if ( !histconfFile.empty() ){
    ifstream ff( histconfFile );
    if ( !ff ){
      cerr << "problem opening " << histconfFile << endl;
      exit(1);
    }
    string line;
    while ( getline( ff, line ) ){
      vector<string> v;
      if ( TiCC::split_at( line, v, "#" ) < 2
	   || TiCC::split_at( line, v, "#" ) > 3 ){
	continue;
      }
      bitType val = TiCC::stringTo<bitType>(v[0]);
      histMap.insert(val);
    }
    if ( histMap.size() == 0 ){
      cerr << "the historical confusions file " << histconfFile
	   << " doesn't seem to be in the right format." << endl
	   << " should contain lines like: 10331739614#f~s" << endl;
    }
    else {
      cout << "read " << histMap.size() << " historical confusions." << endl;
    }
  }

  set<bitType> diaMap;
  if ( !diaconfFile.empty() ){
    ifstream ff( diaconfFile );
    if ( !ff ){
      cerr << "problem opening " << diaconfFile << endl;
      exit(1);
    }
    string line;
    while ( getline( ff, line ) ){
      vector<string> v;
      if ( TiCC::split_at( line, v, "#" ) < 2
	   || TiCC::split_at( line, v, "#" ) > 3 ){
	continue;
      }
      bitType val = TiCC::stringTo<bitType>(v[0]);
      diaMap.insert(val);
    }
    if ( diaMap.size() == 0 ){
      cerr << "the diacritical confusions file " << histconfFile
	   << " doesn't seem to be in the right format." << endl
	   << " should contain lines like: 10331739614#e~é" << endl;
    }
    else {
      cout << "read " << diaMap.size() << " diacritical confusions." << endl;
    }
  }

  ifstream indexf( indexFile );
  if ( !indexf ){
    cerr << "problem opening: " << indexFile << endl;
    exit(1);
  }
  ifstream anaf( anahashFile );
  if ( !anaf ){
    cerr << "problem opening anagram hashes file: " << anahashFile << endl;
    exit(1);
  }
  map<bitType,set<string> > hashMap;
  while ( getline( anaf, line ) ){
    vector<string> v1;
    if ( TiCC::split_at( line, v1, "~" ) != 2 )
      continue;
    else {
      vector<string> v2;
      if ( TiCC::split_at( v1[1], v2, "#" ) < 1 ){
	cerr << "strange line: " << line << endl;
      }
      else {
	bitType key = TiCC::stringTo<bitType>( v1[0] );
	for ( size_t i=0; i < v2.size(); ++i )
	  hashMap[key].insert( v2[i] );
      }
    }
  }
  cout << "read " << hashMap.size() << " hash values" << endl;
#ifdef HAVE_OPENMP
  omp_set_num_threads( numThreads );
#endif

  size_t count=0;
  ofstream os( outFile );
  set<bitType> handledTrans;
  while ( getline( indexf, line ) ){
    if ( verbose ){
      cerr << "bekijk " << line << endl;
    }
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "#" ) != 2 ){
      cerr << "ARGL: " << line << endl;
      exit( EXIT_FAILURE);
    }
    else {
      string mainKeyS = parts[0];
      if ( ++count % 1000 == 0 ){
	cout << ".";
	cout.flush();
	if ( count % 50000 == 0 ){
	  cout << endl << count << endl;;
	}
      }
      string rest = parts[1];
      if ( verbose ){
	cerr << "extract parts from " << rest << endl;
      }
      if ( TiCC::split_at( rest, parts, "," ) < 1 ){
	cerr << "arggl, line=" << line << endl;
	exit( EXIT_FAILURE);
      }
      else {
	bitType mainKey = TiCC::stringTo<bitType>(mainKeyS);
	bool isKHC = false;
	if ( histMap.find( mainKey ) != histMap.end() ){
	  isKHC = true;
	}
	bool isDIAC = false;
	if ( diaMap.find( mainKey ) != diaMap.end() ){
	  isDIAC = true;
	}
	string result;
#pragma omp parallel for schedule(dynamic,1)
	for ( size_t i=0; i < parts.size(); ++i ){
	  string keyS = parts[i];
	  bitType key = TiCC::stringTo<bitType>(keyS);
	  if ( verbose ){
#pragma omp critical
	    cout << "bekijk key1 " << key << endl;
	  }
	  map<bitType,set<string> >::const_iterator sit1 = hashMap.find(key);
	  if ( sit1 == hashMap.end() ){
#pragma omp critical
	    cerr << "WARNING: found a key '" << key
		 << "' in the input that isn't present in the hashes." << endl;
	    continue;
	  }
	  if ( sit1->second.size() > 0
	       && LDvalue >= 2 ){
	    bool do_trans = false;
#pragma omp critical
	    {
	      set<bitType>::const_iterator it = handledTrans.find( key );
	      if ( it == handledTrans.end() ){
		handledTrans.insert( key );
		do_trans = true;
	      }
	    }
	    if ( do_trans ){
	      handleTranspositions( os, sit1->second,
				    freqMap, low_freqMap, alfabet,
				    artifreq, isKHC, noKHCld, isDIAC );
	    }
	  }
	  if ( verbose ){
#pragma omp critical
	    cout << "bekijk key2 " << mainKey + key << endl;
	  }
	  map<bitType, set<string> >::const_iterator sit2 = hashMap.find(mainKey+key);
	  if ( sit2 == hashMap.end() ){
	    if ( verbose ){
#pragma omp critical
	      cerr << "WARNING: found a key '" << key
		   << "' in the input that, when added to '" << mainKey
		   << "' isn't present in the hashes." << endl;
	    }
	    continue;
	  }
	  compareSets( os, LDvalue, mainKeyS,
		       sit1->second, sit2->second,
		       freqMap, low_freqMap, alfabet,
		       artifreq, isKHC, noKHCld, isDIAC );
	  if ( backward ){
	    if ( verbose ){
#pragma omp critical
	      cout << "BACKWARD bekijk key2 " << key - mainKey << endl;
	    }
	    map<bitType,set<string> >::const_iterator sit2 = hashMap.find(key-mainKey);
	    if ( sit2 == hashMap.end() ){
	      if ( verbose ){
#pragma omp critical
		cerr << "WARNING: found a key '" << key
		     << "' in the input that, when substracked from '"
		     << mainKey << "' isn't present in the hashes." << endl;
	      }
	      continue;
	    }
	    compareSets( os, LDvalue, mainKeyS,
			 sit1->second, sit2->second,
			 freqMap, low_freqMap, alfabet,
			 artifreq, isKHC, noKHCld, isDIAC );
	  }
	}
      }
    }
  }
  cout << "Done" << endl;

}
