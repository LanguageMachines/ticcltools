/*
  Copyright (c) 2006 - 2021
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
#include <cassert>
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
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"
#include "config.h"

using namespace std;
using namespace icu;

string progname;
int verbose = 0;

void usage( const string& progname ){
  cerr << "usage: " << progname << endl;
  cerr << "\t--index <confuslist>\t inputfile produced by TICCL-indexer or TICCL-indexerNT." << endl;
  cerr << "\t--hash <anahash>\t a file produced by TICCl-anahash," << endl;
  cerr << "\t--clean <cleanfile>\t a file produced by TICCL-unk" << endl;
  cerr << "\t--diac <diacriticsfile>\t a list of 'diacritical' confusions." << endl;
  cerr << "\t--hist <historicalfile>\t a list of 'historical' confusions." << endl;
  cerr << "\t--alph <alphabet>\t alphabet file (as produced by TICCL-lexstat)" << endl;
  cerr << "\t--nohld\t\t ignore --LD for 'historical' confusions." << endl;
  cerr << "\t-o <outputfile>\t the name of an outputfile." << endl;
  cerr << "\t-t <threads> or --threads <threads> Number of threads to run on." << endl;
  cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  cerr << "\t\t\t reasonable value. (OMP_NUM_TREADS - 2)" << endl;
  cerr << "\t--LD <distance>\t The Levensthein (or edit) distance to use" << endl;
  cerr << "\t--artifrq <artifreq>\t (default=0)" << endl;
  cerr << "\t--low=<low>\t skip entries from the anagram file shorter than "
       << endl;
  cerr << "\t\t\t'low' characters. (default=5)" << endl;
  cerr << "\t--high=<high>\t skip entries from the anagram file longer than "
       << endl;
  cerr << "\t\t\t'high' characters. (default=35)" << endl;
  cerr << "\t-v\t\t be verbose, repeat to be more verbose " << endl;
  cerr << "\t-h or --help\t this message " << endl;
  cerr << "\t-V or --version\t show version " << endl;
  exit( EXIT_FAILURE );
}

set<string> follow_words;
map<UChar,bitType> alphabet;

class ld_record {
public:
  ld_record( const string&,
	     const string&,
	     bitType key1,
	     bitType key2,
	     const map<string,size_t>&,
	     const map<UnicodeString,size_t>&,
	     bool, bool, bool,
	     bool );
  void flip(){
    str1.swap(str2);
    ls1.swap(ls2);
    swap( freq1, freq2 );
    swap( low_freq1, low_freq2 );
  }
  bool analyze_ngrams( const map<UnicodeString, size_t>&,
		       size_t, size_t,
		       map<UnicodeString,set<UnicodeString>>&,
		       map<UnicodeString, size_t>&,
		       map<UnicodeString, size_t>& );
  bool handle_the_pair( const UnicodeString&,
			const UnicodeString&,
			const map<UnicodeString, size_t>&,
			size_t,
			size_t,
			map<UnicodeString,set<UnicodeString>>&,
			map<UnicodeString, size_t>&,
			map<UnicodeString, size_t>& );
  bool ld_is( int );
  bool ld_check( int );
  void fill_fields( size_t );
  void sort_high_second();
  bool test_frequency( size_t );
  bool acceptable( size_t, const map<UChar,bitType>& );
  UnicodeString get_key() const;
  string toString() const;
  UnicodeString str1;
  UnicodeString ls1;
  size_t freq1;
  size_t low_freq1;
  UnicodeString str2;
  UnicodeString ls2;
  size_t freq2;
  size_t low_freq2;
  int ld;
  int cls;
  bitType KWC;
  bitType _key1;
  bitType _key2;
  bool canon;
  bool FLoverlap;
  bool LLoverlap;
  UnicodeString KHC;
  int ngram_point;
  bool isKHC;
  bool noKHCld;
  bool is_diac;
  bool follow;
};


ld_record::ld_record( const string& s1, const string& s2,
		      bitType key1, bitType key2,
		      const map<string,size_t>& f_map,
		      const map<UnicodeString,size_t>& low_f_map,
		      bool is_KHC, bool no_KHCld, bool is_diachrone,
		      bool following ):
  ld(-1),
  cls(0),
  KWC(0),
  _key1(key1),
  _key2(key2),
  canon(false),
  FLoverlap(false),
  LLoverlap(false),
  ngram_point(0),
  isKHC(is_KHC),
  noKHCld(no_KHCld),
  is_diac(is_diachrone)
{
  str1 = TiCC::UnicodeFromUTF8(s1);
  auto const it1 = f_map.find( s1 );
  if ( it1 != f_map.end() ){
    freq1 = it1->second;
  }
  else {
    freq1 = 0;
  }
  ls1 = str1;
  ls1.toLower();
  auto const lit1 = low_f_map.find( ls1 );
  if ( lit1 != low_f_map.end() ){
    low_freq1 = lit1->second;
  }
  else {
    low_freq1 = 0;
  }
  auto const it2 = f_map.find( s2 );
  if ( it2 != f_map.end() ){
    freq2 = it2->second;
  }
  else {
    freq2 = 0;
  }
  str2 = TiCC::UnicodeFromUTF8(s2);
  ls2 = str2;
  ls2.toLower();
  auto const lit2 = low_f_map.find( ls2 );
  if ( lit2 != low_f_map.end() ){
    low_freq2 = lit2->second;
  }
  else {
    low_freq2 = 0;
  }
  follow = following;
}

UnicodeString ld_record::get_key() const {
  return str1 + "~" + str2;
}

bool ld_record::handle_the_pair( const UnicodeString& diff_part1,
				 const UnicodeString& diff_part2,
				 const map<UnicodeString, size_t>& low_freqMap,
				 size_t freqThreshold,
				 size_t low_limit,
				 map<UnicodeString,set<UnicodeString>>& dis_map,
				 map<UnicodeString, size_t>& dis_count,
				 map<UnicodeString, size_t>& ngram_count ){
  //
  // Ok, so we have a pair
  //
  if ( follow ){
#pragma omp critical (debugout)
    {
      cerr << "ngram candidate: '" << diff_part1 << "~" << diff_part2
	   << "' in n-grams pair: " << str1 << " # " << str2 << endl;
    }
  }
  if ( diff_part1.isEmpty() ) {
    // can this happen?
    // anyway: nothing to do
    return false; // nothing special
  }

  UnicodeString lp = diff_part1;
  lp.toLower();
  auto const& entry1 = low_freqMap.find( lp );
  lp = diff_part2;
  lp.toLower();
  if ( entry1 != low_freqMap.end()
       && entry1->second >= freqThreshold ){
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "ngram part1: " << diff_part1 << " is high frequent: "
	     << "skipping" << endl;
      }
    }
    return true; // no use to keep this
  }
  // so this IS a potential good correction
  ngram_point = 1;
  UnicodeString disamb_pair = diff_part1 + "~" + diff_part2;
  if ( (size_t)diff_part1.length() < low_limit ){
    // a 'short' word
    // count this short words pair AND store the original n-gram pair
#pragma omp critical (update)
    {
      dis_map[disamb_pair].insert( str1 + "~" + str2 );
      ++dis_count[disamb_pair];
    }
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "stored: short " << disamb_pair << " and forget about "
	     << str1 << "~" << str2 << endl;
      }
    }
  }
  else {
    // count the pair
#pragma omp critical (update)
    {
      ++ngram_count[disamb_pair];
      // keep pair for later
      //      ngram_count.erase(str1 + "~" + str2);
    }
    // signal to discard this ngram (in favor of the unigram within)
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "stored: " << disamb_pair << " and forget about "
	     << str1 << "~" << str2 << endl;
      }
    }
    //    return false;
  }
  return true; // forget the original parents
}

bool ld_record::analyze_ngrams( const map<UnicodeString, size_t>& low_freqMap,
				size_t freqThreshold,
				size_t low_limit,
				map<UnicodeString,set<UnicodeString>>& dis_map,
				map<UnicodeString, size_t>& dis_count,
				map<UnicodeString, size_t>& ngram_count ){
  ngram_point = 0;
  vector<UnicodeString> parts1 = TiCC::split_at( str1, US_SEPARATOR );
  vector<UnicodeString> parts2 = TiCC::split_at( str2, US_SEPARATOR );
  if ( parts1.size() == 1 && parts2.size() == 1 ){
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "ngram candidates: " << parts1[0] << " AND " << parts2[0]
	     << " are UNIGRAMS: nothing to do" << endl;
      }
    }
    return false; // nothing special for unigrams
  }
  UnicodeString diff_part1;
  UnicodeString diff_part2;
  if ( parts1.size() == parts2.size() ){
    //
    // search for a pair of 'uncommon' parts in the 2 ngrams.
    for ( size_t i=0; i < parts1.size(); ++i ){
      UnicodeString left = parts1[i];
      left.toLower();
      UnicodeString right = parts2[i];
      right.toLower();
      if ( left == right ){
	// ok, a common part.
      }
      else if ( diff_part1.isEmpty() ) {
	// not yet an uncommon part found. store it.
	diff_part1 = parts1[i];
	diff_part2 = parts2[i];
      }
      else {
	// another uncommon part. these n-grams are too uncommon
	if ( follow ){
#pragma omp critical (debugout)
	  {
	    cerr << "ngram candidates: " << str1 << " AND " << str2
		 << " are too different. Discard" << endl;
	  }
	}
	return true; // discard
      }
    }
    return handle_the_pair( diff_part1, diff_part2, low_freqMap,
			    freqThreshold,
			    low_limit,
			    dis_map,
			    dis_count,
			    ngram_count );
  }
  else {
    bool uncommon = true;
    using TiCC::operator<<;
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "analyze ngram candidates: " << str1 << " AND " << str2
	     << endl;
      }
    }

    while( !parts1.empty() && !parts2.empty() ){
      UnicodeString left = parts1.back();
      left.toLower();
      UnicodeString right = parts2.back();
      right.toLower();
      if ( left == right ){
	// remove all common parts at the end.
	parts1.pop_back();
	parts2.pop_back();
	uncommon = false; // signal this
      }
      else {
	break;
      }
    }
    while( !parts1.empty() && !parts2.empty() ){
      UnicodeString left = parts1.front();
      left.toLower();
      UnicodeString right = parts2.front();
      right.toLower();
      if ( left == right ){
	// remove all common parts at the begin.
	parts1.erase(parts1.begin());
	parts2.erase(parts2.begin());
	uncommon = false; // signal this
      }
      else {
	break;
      }
    }
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "after reduction, candidates: " << parts1 << " AND " << parts2
	     << endl;
      }
    }
    if ( uncommon ){
      // no common parts at begin or end.
      if ( follow ){
#pragma omp critical (debugout)
	{
	  cerr << "ngram candidates: " << str1 << " AND " << str2
	       << " are too different. Discard" << endl;
	}
      }
      return false;
    }
    if ( parts1.empty() || parts2.empty() ){
      if ( follow ){
#pragma omp critical (debugout)
	{
	  cerr << "ngram candidate is empty after removing commom parts "
	       << endl;
	}
      }
      return true; // discard
    }
    for ( const auto& w1 : parts1 ){
      if ( !diff_part1.isEmpty() ){
	diff_part1 += "_";
      }
      diff_part1 += w1;
    }
    for ( const auto& w2 : parts2 ){
      if ( !diff_part2.isEmpty() ){
	diff_part2 += "_";
      }
      diff_part2 += w2;
    }
    if ( follow ){
#pragma omp critical (debugout)
      {
	cerr << "FOUND 1-2-3 " << diff_part1 << " " << diff_part2 << endl;
      }
    }
    return handle_the_pair( diff_part1, diff_part2, low_freqMap,
			    freqThreshold,
			    low_limit,
			    dis_map,
			    dis_count,
			    ngram_count );
  }
}

bool ld_record::ld_is( int wanted ) {
  ld = ldCompare( ls1, ls2 );
  if ( ld != wanted ){
    if ( !( isKHC && noKHCld ) ){
      if ( follow ){
#pragma omp critical (debugout)
	{
	  cout << "LD " << ld << " ≠ " << wanted << " and rejected, no KHC"
	       << endl;
	}
      }
      return false;
    }
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << "LD " << ld << " ≠ " << wanted << " but kept, KHC"  << endl;
      }
    }
  }
  if ( follow ){
#pragma omp critical (debugout)
    {
      cout << "LD(" << ls1 << "," << ls2 << ")=" << ld << " OK!" << endl;
    }
  }
  return true;
}

bool ld_record::ld_check( int ldvalue ) {
  ld = ldCompare( ls1, ls2 );
  if ( ld <= ldvalue ){
    // LD is ok
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << "LD(" << ls1 << "," << ls2 << ") =" << ld
	     << " OK,  <= " << ldvalue << endl;
      }
    }
    return true;
  }
  else {
    // reject if the LD exeeds ldvalue AND
    //   not a Historical word OR nohld is not specified (terrible logic this)
    if ( !( isKHC && noKHCld ) ){
      if ( follow ){
#pragma omp critical (debugout)
	{
	  cout << "LD " << ld << " > " << ldvalue << " and rejected, no KHC"
	       << endl;
	}
      }
      return false;
    }
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << "LD " << ld << " > " << ldvalue << " but accepted, KHC"
	     << endl;
      }
    }
    return true;
  }
  if ( follow ){
#pragma omp critical (debugout)
    {
      cout << "LD(" << ls1 << "," << ls2 << ") =" << ld
	   << " rejected > " << ldvalue << endl;
    }
  }
  return false;
}

bool ld_record::acceptable( size_t threshold,
			    const map<UChar,bitType>& alphabet ) {
  if ( low_freq1 >= threshold && !is_diac ){
    // reject correction of lexical words, except for diachrone translations
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << str1 << "~" << str2 << " rejected: Lexical, and not diachrone"
	     << endl;
      }
    }
    return false;
  }
  if ( !alphabet.empty() ){
    // reject non lexically clean Corection Candidates
    for ( int i=0; i < ls2.length(); ++i ){
      if ( alphabet.find( ls2[i] ) == alphabet.end() ){
	if ( follow ){
#pragma omp critical (debugout)
	  {
	    cout << str1 << "~" << str2 << " rejected: "
		 << UnicodeString( ls2[i] ) << " not in alphabet" << endl;
	  }
	}
	return false;
      }
    }
  }
  return true;
}

bool ld_record::test_frequency( size_t threshold ){
  // avoid non lexical Correction Candidates
  if ( low_freq2 < threshold ){
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << str1 << "~" << str2 << " rejected: " << str2
	     << " is low frequent: " << low_freq2 << endl;
      }
    }
    return false;
  }
  return true;
}

void ld_record::sort_high_second(){
  // order the record with the highest (most probable) freqency as CC
  if ( low_freq1 == low_freq2 ){
    //    if ( ::hash(str1,alphabet) > ::hash(str2,alphabet) ){
    if ( _key1 < _key2 ){
      flip();
    }
  }
  else if ( low_freq1 > low_freq2 ){
    if ( follow ){
#pragma omp critical (debugout)
      {
	cout << "flip " << str1 << "~" << str2 << endl;
      }
    }
    flip();
  }
}

void ld_record::fill_fields( size_t freqThreshold ) {
  cls = max(ls1.length(),ls2.length()) - ld;
  LLoverlap = false;
  if ( ls1.length() > 1 && ls2.length() > 1
       && ls1[ls1.length()-1] == ls2[ls2.length()-1]
       && ls1[ls1.length()-2] == ls2[ls2.length()-2] ){
    LLoverlap = true;
  }
  FLoverlap = false;
  if ( ls1[0] == ls2[0] ){
    FLoverlap = true;
  }
  canon = false;
  if ( low_freq2 >= freqThreshold ){
    canon = true;
  }
}

string ld_record::toString() const {
  string canon_s = (canon?"1":"0");;
  string FLoverlap_s = (FLoverlap?"1":"0");;
  string LLoverlap_s = (LLoverlap?"1":"0");;
  string KHC = (isKHC?"1":"0");
  stringstream ss;
  ss << str1 << "~" << freq1 << "~" << low_freq1 << "~"
     << str2 << "~" << freq2 << "~" << low_freq2 << "~"
     << KWC << "~" << ld << "~"
     << cls << "~" << canon_s << "~"
     << FLoverlap_s << "~" << LLoverlap_s << "~"
     << KHC << "~" << ngram_point;
  return ss.str();
}

bool transpose_pair( ld_record& record,
		     const map<UnicodeString,size_t>& low_freqMap,
		     map<UnicodeString,set<UnicodeString>>& dis_map,
		     map<UnicodeString, size_t>& dis_count,
		     map<UnicodeString, size_t>& ngram_count,
		     size_t freqThreshold,
		     size_t low_limit,
		     const map<UChar,bitType>& alphabet,
		     bool following ){
  if ( following ){
#pragma omp critical (debugout)
    {
      cout << "TRANSPOSE: string 1 " << record.str1
	   << " string 2 " << record.str2 << endl;
    }
  }
  record.sort_high_second();
  if ( !record.acceptable( freqThreshold, alphabet ) ){
    return false;
  }
  if ( !record.test_frequency( freqThreshold ) ){
    return false;
  }
  if ( record.analyze_ngrams( low_freqMap, freqThreshold, low_limit,
			      dis_map, dis_count, ngram_count ) ){
    return false;
  }
  if ( !record.ld_is( 2 ) ){
    if ( following ){
#pragma omp critical (debugout)
      {
	cout << " LD != 2 " << record.str1 << "," << record.str2 << endl;
      }
    }
    return false;
  }
  record.fill_fields( freqThreshold );
  if ( following ){
    cerr << "Transpose result: " << record.toString() << endl;
  }
  return true;
}

void handleTranspositions( const set<string>& s,
			   bitType key,
			   const map<string,size_t>& freqMap,
			   const map<UnicodeString,size_t>& low_freqMap,
			   const map<UChar,bitType>& alphabet,
			   map<UnicodeString,set<UnicodeString>>& dis_map,
			   map<UnicodeString, size_t>& dis_count,
			   map<UnicodeString, size_t>& ngram_count,
			   size_t freqThreshold,
			   size_t low_limit,
			   bool isKHC,
			   bool noKHCld,
			   bool isDIAC,
			   map<UnicodeString,ld_record>& record_store ){
  auto it1 = s.begin();
  while ( it1 != s.end() ) {
    bool following = false;
    string str1 = *it1;
    if ( follow_words.find( str1 ) != follow_words.end() ){
      following = true;
    }
    auto it2 = it1;
    ++it2;
    while ( it2 != s.end() ) {
      string str2 = *it2;
      if ( follow_words.find( str2 ) != follow_words.end() ){
	following = true;
      }
      ld_record record( str1, str2,
			key, key,
			freqMap, low_freqMap,
			isKHC, noKHCld, isDIAC, following );
      if ( transpose_pair( record, low_freqMap,
			   dis_map, dis_count, ngram_count,
			   freqThreshold, low_limit, alphabet, following ) ){
	UnicodeString key = record.get_key();
#pragma omp critical (output)
	{
	  if ( following ){
	    if ( record_store.find(key) == record_store.end() ){
	      cerr << "1 insert: " << record.toString() << endl;
	    }
	    else {
	      cerr << "1 emplace: " << record_store.find(key)->second.toString() << endl;
	    }
	  }
	  record_store.emplace(key,record);
	  if ( following ){
	    cerr << "1 emplaced result      : " << record.toString() << endl;
	  }
	}
      }
      ++it2;
    }
    ++it1;
  }
}


bool compare_pair( ld_record& record,
		   const map<UnicodeString,size_t>& low_freqMap,
		   int ldValue,
		   bitType KWC,
		   map<UnicodeString,set<UnicodeString>>& dis_map,
		   map<UnicodeString, size_t>& dis_count,
		   map<UnicodeString, size_t>& ngram_count,
		   size_t freqThreshold,
		   size_t low_limit,
		   const map<UChar,bitType>& alphabet,
		   bool following ){
  if ( !record.ld_check( ldValue ) ){
    return false;
  }
  record.sort_high_second();
  if ( !record.acceptable( freqThreshold, alphabet) ){
    return false;
  }
  if ( record.analyze_ngrams( low_freqMap, freqThreshold, low_limit,
			      dis_map, dis_count, ngram_count ) ){
    return false;
  }
  record.fill_fields( freqThreshold );
  record.KWC = KWC;
  return true;
}

void compareSets( int ldValue,
		  bitType KWC,
		  bitType key1,
		  const set<string>& s1, const set<string>& s2,
		  const map<string,size_t>& freqMap,
		  const map<UnicodeString,size_t>& low_freqMap,
		  const map<UChar,bitType>& alphabet,
		  map<UnicodeString,set<UnicodeString>>& dis_map,
		  map<UnicodeString, size_t>& dis_count,
		  map<UnicodeString, size_t>& ngram_count,
		  size_t freqThreshold,
		  size_t low_limit,
		  bool isKHC,
		  bool noKHCld,
		  bool isDIAC,
		  map<UnicodeString,ld_record>& record_store ){
  // using TiCC::operator<<;
  // cerr << "set 1 " << s1 << endl;
  // cerr << "set 2 " << s2 << endl;
  auto it1 = s1.begin();
  while ( it1 != s1.end() ) {
    bool following = false;
    string str1 = *it1;
    if ( follow_words.find( str1 ) != follow_words.end() ){
      following = true;
    }
    if ( following ){
#pragma omp critical (debugout)
      {
	cout << "SET: string 1 " << str1 << endl;
      }
    }
    auto it2 = s2.begin();
    while ( it2 != s2.end() ) {
      string str2 = *it2;
      if ( follow_words.find( str2 ) != follow_words.end() ){
	following = true;
      }
      if ( following ){
#pragma omp critical (debugout)
	{
	  cout << "SET: string 2 " << str2 << endl;
	}
      }
      ld_record record( str1, str2,
			key1, KWC + key1,
			freqMap, low_freqMap,
			isKHC, noKHCld, isDIAC, following );
      if ( compare_pair( record, low_freqMap, ldValue, KWC,
			 dis_map, dis_count, ngram_count,
			 freqThreshold, low_limit, alphabet, following ) ){
	UnicodeString key = record.get_key();
#pragma omp critical (output)
	{
	  if ( following ){
	    if ( record_store.find(key) == record_store.end() ){
	      cerr << "2 insert: " << record.toString() << endl;
	    }
	    else {
	      cerr << "2 emplace: "
		   << record_store.find(key)->second.toString() << endl
		   << " By      : " << record.toString() << endl;
	    }
	  }
	  record_store.emplace(key,record);
	}
      }
      ++it2;
    }
    ++it1;
  }
}

void add_short( ostream& os,
		const map<UnicodeString,size_t>& dis_count,
		const map<string,size_t>& freqMap,
		const map<UnicodeString,size_t>& low_freqMap,
		int max_ld, size_t threshold ){
  for ( const auto& entry : dis_count ){
    vector<UnicodeString> parts = TiCC::split_at( entry.first, "~" );
    ld_record rec( TiCC::UnicodeToUTF8(parts[0]),
		   TiCC::UnicodeToUTF8(parts[1]),
		   0, 0,
		   freqMap, low_freqMap,
		   false, false, false, false );
    if ( !rec.ld_check( max_ld ) ){
      continue;
    }
    rec.fill_fields( threshold );
    rec.ngram_point = entry.second;
    os << rec.toString() << endl;
  }
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "diac:,hist:,nohld,artifrq:,LD:,hash:,clean:,alph:,"
			   "index:,help,version,threads:,follow:,low:,high:" );
    opts.init( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    progname = opts.prog_name();
    cerr << e.what() << endl;
    usage( progname );
    exit( EXIT_FAILURE );
  }
  progname = opts.prog_name();
  if ( argc < 2	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  if ( opts.extract('h') || opts.extract("help") ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V') || opts.extract("version") ){
    cerr << progname << ": " << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  while ( opts.extract( 'v' ) ){
    ++verbose;
  }
  string value;
  while ( opts.extract( "follow", value ) ){
    vector<string> parts = TiCC::split_at( value, "," );
    for ( const auto& p : parts ){
      follow_words.insert( p );
    }
  }

  string indexFile;
  string anahashFile;
  string frequencyFile;
  string histconfFile;
  string diaconfFile;
  string alfabetFile;
  int LDvalue=2;
  bool noKHCld = opts.extract("nohld");
  if ( !opts.extract( "index", indexFile ) ){
    cerr << progname << ": missing --index option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !TiCC::match_back( indexFile, ".index" )
       && !TiCC::match_back( indexFile, ".indexNT" ) ){
    cerr << progname << ": --index files must have extension: '.index' or '.indexNT' "
	 << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "hash", anahashFile ) ){
    cerr << progname << ": missing --hash option" << endl;
    exit( EXIT_FAILURE );
  }
  if ( !opts.extract( "clean", frequencyFile ) ){
    cerr << progname << ": missing --clean option" << endl;
    exit( EXIT_FAILURE );
  }
  opts.extract( "alph", alfabetFile );
  opts.extract( "hist", histconfFile );
  if ( opts.extract( "diac", diaconfFile ) ){
    if ( !TiCC::match_back( diaconfFile, ".diac" ) ){
      cerr << progname << ": invalid extension for --diac file '" << diaconfFile
	   << "' (must be .diac) " << endl;
      exit(EXIT_FAILURE);
    }
  }
  string outFile;
  string shortFile;
  if ( opts.extract( 'o', outFile ) ){
    if ( !TiCC::match_back( outFile, ".ldcalc" ) ){
      shortFile = outFile + ".short.ldcalc";
      outFile += ".ldcalc";
    }
    else {
      shortFile = outFile;
      shortFile.insert( shortFile.length() - 7, ".short" );
    }
  }
  else {
    outFile = indexFile + ".ldcalc";
    shortFile = indexFile + ".short.ldcalc";
  }
  string ambiFile = outFile + ".ambi";
  size_t artifreq = 0;

  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << progname << ": illegal value for --artifrq (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  int low_limit = 5;
  if ( opts.extract( "low", value ) ){
    if ( !TiCC::stringTo(value,low_limit) ){
      cerr << progname << ": illegal value for --low (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }

  int high_limit = 35;
  if ( opts.extract( "high", value ) ){
    if ( !TiCC::stringTo(value,high_limit) ) {
      cerr << progname << ": illegal value for --high (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
#ifdef HAVE_OPENMP
  int numThreads;
  if ( TiCC::lowercase(value) == "max" ){
    numThreads = omp_get_max_threads() - 2;
    omp_set_num_threads( numThreads );
    cout << "running on " << numThreads << " threads." << endl;
  }
  else {
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    omp_set_num_threads( numThreads );
    cout << "running on " << numThreads << " threads." << endl;
  }
#else
  if ( value != "1" ){
    cerr << "unable to set number of threads!.\nNo OpenMP support available!"
	 <<endl;
    exit(EXIT_FAILURE);
  }
#endif
  if ( opts.extract( "LD", value ) ){
    if ( !TiCC::stringTo(value,LDvalue) ) {
      cerr << progname << ": illegal value for --LD (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
    if ( LDvalue < 1 || LDvalue > 10 ){
      cerr << progname << ": invalid LD value: " << LDvalue << " (1-10 is OK)" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( !opts.empty() ){
    cerr << progname << ": unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

  if ( !alfabetFile.empty() ){
    ifstream lexicon( alfabetFile );
    if ( !lexicon ){
      cerr << progname << ": problem opening alfabet file: " << alfabetFile << endl;
      exit(EXIT_FAILURE);
    }
    cout << progname << ": reading alphabet: " << alfabetFile << endl;
    fillAlphabet( lexicon, alphabet );
    cout << progname << ": read " << alphabet.size() << " letters with frequencies" << endl;
  }
  ifstream ff( frequencyFile  );
  if ( !ff ){
    cerr << progname << ": problem opening " << frequencyFile << endl;
    exit(EXIT_FAILURE);
  }
  cout << progname << ": reading clean file: " << frequencyFile << endl;
  map<string, size_t> freqMap;
  map<UnicodeString, size_t> low_freqMap;
  string line;
  size_t ign = 0;
  size_t skipped = 0;
  while ( getline( ff, line ) ){
    vector<string> v1;
    if ( TiCC::split( line, v1 ) != 2 ){
      ++ign;
      continue;
    }
    else {
      string s = v1[0];
      UnicodeString ls = TiCC::UnicodeFromUTF8(s);
      if ( low_limit > 0 && ls.length() < low_limit ){
	++skipped;
	continue;
      }
      if ( high_limit > 0 && ls.length() > high_limit ){
	++skipped;
	continue;
      }
      size_t freq = TiCC::stringTo<size_t>( v1[1] );
      freqMap[s] = freq;
      ls.toLower();
      if ( freq >= artifreq ){
	// make sure that the artifrq is counted only once!
	if ( low_freqMap[ls] == 0 ){
	  low_freqMap[ls] = freq;
	}
	else {
	  low_freqMap[ls] += freq-artifreq;
	}
      }
      else {
	low_freqMap[ls] +=freq;
      }
    }
  }
  cout << progname << ": read " << freqMap.size()
       << " clean words with frequencies." << endl;
  if ( skipped > 0 ){
    cout << progname << ": skipped " << skipped << " out-of-band words."
	 << endl;
  }
  if ( ign > 0 ){
    cout << progname << ": skipped " << ign << " spaced words in the clean file" << endl;
  }
  set<bitType> histMap;
  if ( !histconfFile.empty() ){
    ifstream ff( histconfFile );
    if ( !ff ){
      cerr << "problem opening " << histconfFile << endl;
      exit(EXIT_FAILURE);
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
    if ( histMap.empty() ){
      cerr << progname << ": the historical confusions file " << histconfFile
	   << " doesn't seem to be in the right format." << endl
	   << " should contain lines like: 10331739614#f~s" << endl;
    }
    else {
      cout << progname << ": read " << histMap.size() << " historical confusions." << endl;
    }
  }

  set<bitType> diaMap;
  if ( !diaconfFile.empty() ){
    ifstream ff( diaconfFile );
    if ( !ff ){
      cerr << progname << ": problem opening " << diaconfFile << endl;
      exit(EXIT_FAILURE);
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
    if ( diaMap.empty() ){
      cerr << progname << ": the diacritical confusions file " << histconfFile
	   << " doesn't seem to be in the right format." << endl
	   << " should contain lines like: 10331739614#e~é" << endl;
      exit(EXIT_FAILURE);
    }
    else {
      cout << progname << ": read " << diaMap.size() << " diacritical confusions." << endl;
    }
  }

  ifstream indexf( indexFile );
  if ( !indexf ){
    cerr << progname << ": problem opening: " << indexFile << endl;
    exit(EXIT_FAILURE);
  }
  ifstream anaf( anahashFile );
  if ( !anaf ){
    cerr << progname << ": problem opening anagram hashes file: " << anahashFile << endl;
    exit(EXIT_FAILURE);
  }
  map<bitType,set<string> > hashMap;
  while ( getline( anaf, line ) ){
    vector<string> v1;
    if ( TiCC::split_at( line, v1, "~" ) != 2 )
      continue;
    else {
      vector<string> v2;
      if ( TiCC::split_at( v1[1], v2, "#" ) < 1 ){
	cerr << progname << ": strange line: " << line << endl
	     << " in anagram hashes file" << endl;
	exit(EXIT_FAILURE);
      }
      else {
	bitType key = TiCC::stringTo<bitType>( v1[0] );
	for ( size_t i=0; i < v2.size(); ++i ){
	  auto it = freqMap.find( v2[i] );
	  if ( it != freqMap.end() ){
	    // only store words from the .clean lexicon
	    hashMap[key].insert( v2[i] );
	  }
	  else {
	    if ( verbose > 1 ){
	      cerr << "skip hash for " << v2[i] << " (not in lexicon)" << endl;
	    }
	  }
	}
      }
    }
  }
  cout << progname << ": read " << hashMap.size() << " hash values" << endl;

  size_t count=0;
  set<bitType> handledTrans;
  map<UnicodeString,set<UnicodeString>> dis_map;
  map<UnicodeString,size_t> dis_count;
  map<UnicodeString,size_t> ngram_count;
  map<UnicodeString,ld_record> record_store;
  size_t line_nr = 0;
  int err_cnt = 0;

  size_t file_lines = 0;
  while ( getline( indexf, line ) ){
    ++file_lines;
  }
  if ( file_lines == 0 ){
    cerr << "the indexfile: '" << indexFile
	 << "' is empty! No further processing possible." << endl;
    exit( EXIT_FAILURE );
  }
  cout << progname << ": " << file_lines << " character confusion values to be read.\n\t\tWe indicate progress by printing a dot for every 1000 confusion values processed" << endl;
  indexf.clear();
  indexf.seekg( 0 );
  while ( getline( indexf, line ) ){
    if ( err_cnt > 9 ){
      cerr << progname << ": FATAL ERROR: too many problems in indexfile: " << indexFile
	   << " terminated" << endl;
      exit( EXIT_FAILURE);
    }
    ++line_nr;
    if ( verbose > 1 ){
      cerr << "examine " << line << endl;
    }
    line = TiCC::trim(line);
    if ( line.empty() ){
      continue;
    }
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "#" ) != 2 ){
      cerr << progname << ": ERROR in line " << line_nr
	   << " of the indexfile: unable to split in 2 parts at #"
	   << endl << "line was" << endl << line << endl;
      ++err_cnt;
    }
    else {
      string key_s = parts[0];
      if ( ++count % 1000 == 0 ){
	cout << ".";
	cout.flush();
	if ( count % 50000 == 0 ){
	  cout << endl << count << endl;;
	}
      }
      string rest = parts[1];
      if ( verbose > 1 ){
	cerr << "extract parts from " << rest << endl;
      }
      if ( TiCC::split_at( rest, parts, "," ) < 1 ){
	cerr << progname << ": ERROR in line " << line_nr
	     << " of indexfile: unable to split in parts separated by ','"
	     << endl << "line was" << endl << line << endl;
	++err_cnt;
      }
      else {
	bitType mainKey = TiCC::stringTo<bitType>(key_s);
	bool isKHC = false;
	if ( histMap.find( mainKey ) != histMap.end() ){
	  isKHC = true;
	}
	bool isDIAC = false;
	if ( diaMap.find( mainKey ) != diaMap.end() ){
	  isDIAC = true;
	}
#pragma omp parallel for schedule(dynamic,1)
	for ( size_t i=0; i < parts.size(); ++i ){
	  string keyS = parts[i];
	  bitType key = TiCC::stringTo<bitType>(keyS);
	  map<bitType,set<string> >::const_iterator sit1 = hashMap.find(key);
	  if ( sit1 == hashMap.end() ){
	    if ( verbose > 1 ){
#pragma omp critical (debugout)
	      cerr << progname << ": WARNING: found a key '" << key
		   << "' in the input that isn't present in the hashes." << endl;
	    }
	    continue;
	  }
	  if ( verbose > 1 ){
#pragma omp critical (debugout)
	    cout << "bekijk key1 " << key << endl;
	  }
	  if ( sit1->second.size() > 0
	       && LDvalue >= 2 ){
	    bool do_trans = false;
#pragma omp critical (debugout)
	    {
	      set<bitType>::const_iterator it = handledTrans.find( key );
	      if ( it == handledTrans.end() ){
		handledTrans.insert( key );
		do_trans = true;
	      }
	    }
	    if ( do_trans ){
	      handleTranspositions( sit1->second,
				    key,
				    freqMap, low_freqMap, alphabet,
				    dis_map, dis_count, ngram_count,
				    artifreq, low_limit, isKHC, noKHCld, isDIAC,
				    record_store );
	    }
	  }
	  map<bitType, set<string> >::const_iterator sit2 = hashMap.find(mainKey+key);
	  if ( sit2 == hashMap.end() ){
	    if ( verbose > 4 ){
#pragma omp critical (debugout)
	      cerr << progname << ": WARNING: found a key '" << key
		   << "' in the input that, when added to '" << mainKey
		   << "' isn't present in the hashes." << endl;
	    }
	    continue;
	  }
	  if ( verbose > 1 ){
#pragma omp critical (debugout)
	    cout << "bekijk key2 " << mainKey + key << endl;
	  }
	  compareSets( LDvalue, mainKey, key,
		       sit1->second, sit2->second,
		       freqMap, low_freqMap, alphabet,
		       dis_map, dis_count, ngram_count,
		       artifreq, low_limit, isKHC, noKHCld, isDIAC,
		       record_store );
	}
      }
    }
  }
  cout << endl << "creating .short file: " << shortFile << endl;
  ofstream shortf( shortFile );
  add_short( shortf, dis_count, freqMap, low_freqMap, LDvalue, artifreq );
  cout << endl << "creating .ambi file: " << ambiFile << endl;
  ofstream amb( ambiFile );
  for ( const auto& ambi : dis_map ){
    amb << ambi.first << "#";
    for ( const auto& val : ambi.second ){
      amb << val << "#";
    }
    amb << endl;
  }
  map<UnicodeString,unsigned int> low_ngramcount;
  for ( const auto& ng : ngram_count ){
    UnicodeString lv = ng.first;
    lv.toLower();
    low_ngramcount[lv] += ng.second;
  }
  for ( const auto& it : ngram_count ){
    if ( record_store.find( it.first ) != record_store.end() ){
      UnicodeString lv = it.first;
      lv.toLower();
      assert( low_ngramcount.find( lv ) != low_ngramcount.end() );
      record_store.find(it.first)->second.ngram_point += low_ngramcount[lv];
    }
    else {
      // Ok, our data seems to be incomplete
      // that is not our problem, so ignore
      if ( verbose > 1 ){
	cerr << "ignoring " << it.first << endl;
      }
    }
  }
  ofstream os( outFile );
  for ( const auto& r : record_store ){
    os << r.second.toString() << endl;
  }
  cout << progname << ": Done" << endl;
}
