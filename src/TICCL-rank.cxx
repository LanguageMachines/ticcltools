/*
  Copyright (c) 2006 - 2018
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
#include <algorithm>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>
#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/Unicode.h"
#include "ticcl/unicode.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace icu;
using TiCC::operator<<;

typedef signed long int bitType;

const int RANK_COUNT=14;
const string SEPARATOR = "_";
set<string> follow_words;

bool verbose = false;

void usage( const string& name ){
  cerr << "usage: " << name << " --alph <alphabetfile> --charconf <lexstat file>[--wordvec <wordvectorfile>] [-o <outputfile>] [-t threads] [--clip <clip>] [--debugfile <debugfile>] [--artifrq art] [--skipcols <skip>] infile" << endl;
  cerr << "\t'infile'\t is a file in TICCL-LDcalc format" << endl;
  cerr << "\t--alph 'alpha'\t an alphabet file in TICCL-lexstat format." << endl;
  cerr << "\t--charconf 'charconfus'\t a character confusion file in TICCL-lexstat format." << endl;
  cerr << "\t--charconfreq 'name'\t Extract a character confusion frequency file" << endl;
  cerr << "\t--wordvec<wordvecfile> read in a google word2vec file." << endl;
  cerr << "\t-o 'outfile'\t name of the output file." << endl;
  cerr << "\t-t <threads>\n\t--threads <threads> Number of threads to run on." << endl;
  cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  cerr << "\t\t\t reasonable value. (OMP_NUM_TREADS - 2)" << endl;
  cerr << "\t--clip 'clip'\t limit the number of variants per word to 'clip'." << endl;
  cerr << "\t--debugfile 'debug'\t produce a more verbose outputfile in parallel." << endl;
  cerr << "\t\t\t (for debugging.)" << endl;
  cerr << "\t--subtractartifrqfeature1 'arti'\t decrease the frequencies for feature 1 with value 'arti'." << endl;
  cerr << "\t\t\t (which should match the artifreq used in TICCL-LDcalc)" << endl;
  cerr << "\t--subtractartifrqfeature2 'arti'\t decrease the frequencies for feature 2 with value 'arti'." << endl;
  cerr << "\t\t\t (which should match the artifreq used in TICCL-LDcalc)" << endl;
  cerr << "\t--artifrq 'arti'\t OBSOLETE. use --subtractartifrqfeature2." << endl;
  cerr << "\t--skipcols=arglist\t skip the named columns in the ranking." << endl;
  cerr << "\t\t\t e.g. if arglist=3,9, then the columns 3 and 9 are not used." << endl;
  cerr << "\t-v\t\t run (very) verbose" << endl;
  exit( EXIT_FAILURE );
}

class record {
public:
  record( const string &, size_t, size_t, const vector<word_dist>& );
  string extractResults() const;
  string extractLong( const vector<bool>& skip ) const;
  string variant;
  string candidate;
  string lower_candidate;
  double variant_count;
  double variant_rank;
  size_t variant_freq;
  size_t low_variant_freq;
  size_t candidate_freq;
  size_t reduced_candidate_freq;
  size_t low_candidate_freq;
  double freq_rank;
  bitType kwc;
  size_t f2len;
  size_t f2len_rank;
  int ld;
  double ld_rank;
  int cls;
  double cls_rank;
  int canon;
  double canon_rank;
  size_t pairs1;
  double pairs1_rank;
  size_t pairs2;
  double pairs2_rank;
  size_t median;
  double median_rank;
  int fl;
  double fl_rank;
  int ll;
  double ll_rank;
  int khc;
  double khc_rank;
  double cosine;
  double cosine_rank;
  int ngram_points;
  double ngram_rank;
  double rank;
};

float lookup( const vector<word_dist>& vec,
	      const string& word ){
  for( size_t i=0; i < vec.size(); ++i ){
    if ( vec[i].w == word ){
      //	cerr << "JA! " << vec[i].d << endl;
      return vec[i].d;
    }
  }
  return 0.0;
}

record::record( const string& line,
		size_t sub_artifreq,
		size_t sub_artifreq_f2,
		const vector<word_dist>& WV ):
  variant_count(-1),
  f2len_rank(-1),
  ld(-1),
  pairs1(0),
  pairs1_rank(-1),
  pairs2(0),
  pairs2_rank(-1),
  median(0),
  median_rank(-1),
  rank(-10000)
{
  vector<string> parts = TiCC::split_at( line, "~" );
  // file a record with the RANK_COUNT parts of one line from a LDcalc output file
  if ( parts.size() == RANK_COUNT ){
    variant = parts[0];
    variant_freq = TiCC::stringTo<size_t>(parts[1]);
    low_variant_freq = TiCC::stringTo<size_t>(parts[2]);
    candidate = parts[3];
    UnicodeString us = TiCC::UnicodeFromUTF8( candidate );
    us.toLower();
    lower_candidate = TiCC::UnicodeToUTF8( us );
    variant_rank = -2000;  // bogus value, is set later
    string f2_string = parts[4];
    candidate_freq = TiCC::stringTo<size_t>( f2_string );
    reduced_candidate_freq = candidate_freq;
    if ( sub_artifreq > 0 && reduced_candidate_freq >= sub_artifreq ){
      reduced_candidate_freq -= sub_artifreq;
    }
    if ( sub_artifreq_f2 > 0 && candidate_freq >= sub_artifreq_f2 ){
      size_t rf2 = candidate_freq - sub_artifreq_f2;
      string rf2_string = TiCC::toString( rf2 );
      f2len = rf2_string.length();
      //      cerr << "f2len is nu: " << f2len << " (" << rf2_string << ") was " <<  f2_string.length() << " (" << f2_string << ")" << endl;
    }
    else {
      f2len = f2_string.length();
    }
    low_candidate_freq = TiCC::stringTo<size_t>(parts[5]);
    freq_rank = -20;  // bogus value, is set later
    kwc   = TiCC::stringTo<bitType>(parts[6]);
    ld    = TiCC::stringTo<int>(parts[7]);
    ld_rank = -4.5;  // bogus value, is set later
    cls   = TiCC::stringTo<int>(parts[8]);
    cls_rank = -5.6; // bogus value, is set later
    canon = TiCC::stringTo<int>(parts[9]);
    if ( canon == 0 )
      canon_rank = 10;
    else
      canon_rank = 1;
    fl = TiCC::stringTo<int>(parts[10]);
    if ( fl == 0 )
      fl_rank = 2;
    else
      fl_rank = 1;
    ll = TiCC::stringTo<int>(parts[11]);
    if ( ll == 0 )
      ll_rank = 2;
    else
      ll_rank = 1;
    khc = TiCC::stringTo<int>(parts[12]);
    if ( khc == 0 )
      khc_rank = 2;
    else
      khc_rank = 1;
    ngram_points = TiCC::stringTo<int>(parts[13]);
    ngram_rank = -6.7;  // bogus value, is set later
    cosine = lookup( WV, candidate );
    if ( cosine <= 0.001 )
      cosine_rank = 1;
    else
      cosine_rank = 10;
  }
}

string record::extractLong( const vector<bool>& skip ) const {
  string result = variant + "#";
  result += TiCC::toString(variant_freq) + "#";
  result += TiCC::toString(low_variant_freq) + "#";
  result += candidate + "#";
  result += TiCC::toString(candidate_freq) + "#";
  result += TiCC::toString(low_candidate_freq) + "#";
  result += TiCC::toString(kwc) + "#";
  result += TiCC::toString(f2len) + "~";
  double the_rank = 0;
  if ( skip[0] ){
    result += "N#";  }
  else {
    the_rank += f2len_rank;
    result += TiCC::toString(f2len_rank) + "#";
  }
  result += TiCC::toString(reduced_candidate_freq) + "~";
  if ( skip[1] ){
    result += "N#";
  }
  else {
    the_rank += freq_rank;
    result += TiCC::toString(freq_rank) + "#";
  }
  result += TiCC::toString(ld) + "~";
  if ( skip[2] ){
    result += "N#";
  }
  else {
    the_rank += ld_rank;
    result += TiCC::toString(ld_rank) + "#";
  }
  result += TiCC::toString(cls) + "~";
  if ( skip[3] ){
    result += "N#";
  }
  else {
    the_rank += cls_rank;
    result += TiCC::toString(cls_rank) + "#";
  }
  result += TiCC::toString(canon) + "~";
  if ( skip[4] ){
    result += "N#";
  }
  else {
    the_rank += canon_rank;
    result += TiCC::toString(canon_rank) + "#";
  }
  result += TiCC::toString(fl) + "~";
  if ( skip[5] ){
    result += "N#";
  }
  else {
    the_rank += fl_rank;
    result += TiCC::toString(fl_rank) + "#";
  }
  result += TiCC::toString(ll) + "~";
  if ( skip[6] ){
    result += "N#";
  }
  else {
    the_rank += ll_rank;
    result += TiCC::toString(ll_rank) + "#";
  }
  result += TiCC::toString(khc) + "~";
  if ( skip[7] ){
    result += "N#";
  }
  else {
    the_rank += khc_rank;
    result += TiCC::toString(khc_rank) + "#";
  }
  result += TiCC::toString(pairs1) + "~";
  if ( skip[8] ){
    result += "N#";
  }
  else {
    the_rank += pairs1_rank;
    result += TiCC::toString(pairs1_rank) + "#";
  }
  result += TiCC::toString(pairs2) + "~";
  if ( skip[9] ){
    result += "N#";
  }
  else {
    the_rank += pairs2_rank;
    result += TiCC::toString(pairs2_rank) + "#";
  }
  result += TiCC::toString(median) + "~";
  if ( skip[10] ){
    result += "N#";
  }
  else {
    the_rank += median_rank;
    result += TiCC::toString(median_rank) + "#";
  }
  result += TiCC::toString(variant_count) + "~";
  if ( skip[11] ){
    result += "N#";
  }
  else {
    the_rank += variant_rank;
    result += TiCC::toString(variant_rank) + "#";
  }
  result += TiCC::toString(cosine) + "~";
  if ( skip[12] ){
    result += "N#";
  }
  else {
    the_rank += cosine_rank;
    result += TiCC::toString(cosine_rank) + "#";
  }
  result += TiCC::toString(ngram_points) + "~";
  if ( skip[13] ){
    result += "N#";
  }
  else {
    the_rank += ngram_rank;
    result += TiCC::toString(ngram_rank) + "#";
  }
  result += TiCC::toString(the_rank) + "#";
  result += TiCC::toString(rank);
  return result;
}

string record::extractResults() const {
  string result = variant + "#";
  result += TiCC::toString(variant_freq) + "#";
  result += candidate + "#";
  result += TiCC::toString(candidate_freq) + "#";
  result += TiCC::toString(ld) + "#";
  result += TiCC::toString(rank);
  return result;
}

template <typename TObject, typename TMember, typename TValue>
void set_val( TObject& object, TMember member, TValue value )
{
    ( object ).*member = value;
}

template< class Tmap, typename TMember >
void rank_desc_map( const Tmap& desc_map,
		    vector<record>& recs,
		    TMember member ){
  // the map is a (multi-)map, sorted descending on the first value
  // whichs is an integer value.
  // the second value of the map is an index in the vector of records.
  // so the map contains records, sorted descending on specific values
  // e.g:
  //  clsmap holds the CommonSubstringLengths of all vectors, longest first
  if ( desc_map.empty() ){
    return;
  }
  size_t last = desc_map.begin()->first; // start with the longest
  int ranking = 1;                    // it will be ranked 1
  for ( const auto& rit : desc_map ){
    if ( rit.first < last ){
      // we find a shorter. so ranking is incremented (meaning LOWER ranking)
      last = rit.first;
      ++ranking;
    }
    // set the currect ranking for the record at hand
    set_val( recs[rit.second], member, ranking );
  }
}

void rank( vector<record>& recs,
	   map<string,multimap<double,record,std::greater<double>>>& results,
	   int clip,
	   const map<bitType,size_t>& kwc_counts,
	   const map<bitType,size_t>& kwc2_counts,
	   const map<bitType,size_t>& kwc_medians,
	   ostream* db, vector<bool>& skip, int factor ){
  bool follow = follow_words.find(recs.begin()->variant) != follow_words.end();
  if ( follow||verbose ){
#pragma omp critical (log)
    {
#ifdef HAVE_OPENMP
      int numt = omp_get_thread_num();
      cerr << numt << "-";
#endif
      cerr << "RANK " << recs[0].variant
	   << " with " << recs.size() << " variants" << endl;
    }
  }
  multimap<size_t,size_t,std::greater<size_t>> freqmap;  // freqs sorted descending
  multimap<size_t,size_t,std::greater<size_t>> f2lenmap; // f2 lenghts sorted descending
  multimap<size_t,size_t> ldmap;
  multimap<size_t,size_t, std::greater<size_t>> clsmap; // Common substring lengths descending
  multimap<size_t,size_t,std::greater<size_t>> pairmap1;
  multimap<size_t,size_t,std::greater<size_t>> pairmap2;
  multimap<size_t,size_t,std::greater<size_t>> median_map;
  multimap<size_t,size_t,std::greater<size_t>> ngram_map;
  map<string,int> lowvarmap;
  size_t count = 0;

  for ( auto& it : recs ){
    // for every record, we store information in descending multimaps
    // So in freqmap, the (index of) the records with highest freq
    //   are stored in front
    //same for f2len, ld, cls and ngram points
    freqmap.insert( make_pair(it.reduced_candidate_freq, count ) ); // freqs descending
    f2lenmap.insert( make_pair(it.f2len, count ) ); // f2lengths descending
    ldmap.insert( make_pair(it.ld,count) ); // lds sorted ASCENDING
    clsmap.insert( make_pair(it.cls,count) ); // cls sorted descending
    ngram_map.insert( make_pair(it.ngram_points,count) ); // ngrampoints sorted descending
    size_t var1_cnt = kwc_counts.at(it.kwc);
    it.pairs1 = var1_cnt;
    pairmap1.insert( make_pair(var1_cnt,count )); // #variants descending
    size_t var2_cnt = 0;
    try {
      var2_cnt += kwc2_counts.at(it.kwc);
    }
    catch(...){
    }
    it.pairs2 = var2_cnt;
    pairmap2.insert( make_pair(var2_cnt,count )); // #variants decending
    it.median = kwc_medians.at(it.kwc);
    median_map.insert( make_pair(it.median,count )); // #medians decending
    ++lowvarmap[it.lower_candidate]; // count frequency of variants
    ++count;
  }
  multimap<int,size_t,std::greater<int>> lower_variantmap; // descending map
  count = 0;
  for ( const auto& it : recs ){
    lower_variantmap.insert( make_pair( lowvarmap[it.lower_candidate], count ) );
    ++count;
  }
  if ( follow ){
    cout << "1 f2lenmap = " << f2lenmap << endl;
    cout << "2 freqmap = " << freqmap << endl;
    cout << "3 ldmap = " << ldmap << endl;
    cout << "4 clsmap = " << clsmap << endl;
    cout << "9 pairmap1 = " << pairmap1 << endl;
    cout << "10 pairmap2 = " << pairmap2 << endl;
    cout << "11 medianmap = " << median_map << endl;
    //    cout << "12-a lowvarmap = " << lowvarmap << endl;
    cout << "12 lower_variantmap = " << lower_variantmap << endl;
    cout << "14 ngram_map = " << ngram_map << endl;
  }
  rank_desc_map( f2lenmap, recs, &record::f2len_rank );
  if ( follow ){
    cout << "step 1: f2len_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.f2len_rank << endl;
    }
  }

  rank_desc_map( freqmap, recs, &record::freq_rank );
  if ( follow ){
    cout << "step 2: freq_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.freq_rank << endl;
    }
  }

  if ( !ldmap.empty() ){
    int ranking = 1;
    size_t last = ldmap.begin()->first;
    for ( const auto& sit : ldmap ){
      if ( sit.first > last ){
   	last = sit.first;
   	++ranking;
      }
      recs[sit.second].ld_rank = ranking;
    }
  }
  if ( follow ){
    cout << "step 3: ld_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ld_rank << endl;
    }
  }

  rank_desc_map( clsmap, recs, &record::cls_rank );
  if ( follow ){
    cout << "step 4: cls_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.cls_rank << endl;
    }
  }

  if ( follow ){
    cout << "step 5: canon_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.canon_rank << endl;
    }
  }

  if ( follow ){
    cout << "step 6: fl_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.fl_rank << endl;
    }
  }

  if ( follow ){
    cout << "step 7: ll_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ll_rank << endl;
    }
  }

  if ( follow ){
    cout << "step 8: khc_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.khc_rank << endl;
    }
  }

  rank_desc_map( pairmap1, recs, &record::pairs1_rank );
  if ( follow ){
    cout << "step 9: pairs1_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.pairs1_rank << endl;
    }
  }

  rank_desc_map( pairmap2, recs, &record::pairs2_rank );
  if ( follow ){
    cout << "step 10: pairs2_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.pairs2_rank << endl;
    }
  }

  rank_desc_map( median_map, recs, &record::median_rank );
  if ( follow ){
    cout << "step 11: median_rank: for " << recs.begin()->variant << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.median_rank << endl;
    }
  }

  if ( !lower_variantmap.empty() ){
    int ranking = 1;
    int last = lower_variantmap.begin()->first;
    for ( const auto& it1 : lower_variantmap ){
      if ( it1.first < last ){
	last = it1.first;
	++ranking;
      }
      recs[it1.second].variant_count = it1.first;
      recs[it1.second].variant_rank = ranking;
    }
  }
  if ( follow ){
    cout << "step 12: lower_variant_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " count=" << r.variant_count << " rank= " << r.variant_rank << endl;
    }
  }

  if ( follow ){
    cout << "step 13: cosine_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.cosine_rank << endl;
    }
  }

  rank_desc_map( ngram_map, recs, &record::ngram_rank );
  if ( follow ){
    cout << "step 14: ngram_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ngram_rank << endl;
    }
  }

  double sum = 0.0;
  vector<record>::iterator vit = recs.begin();
  while ( vit != recs.end() ){
    double rank =
      (skip[0]?0:(*vit).f2len_rank) +  // number of characters in the frequency
      (skip[1]?0:(*vit).freq_rank) +   // frequency of the CC
      (skip[2]?0:(*vit).ld_rank) +     // levenshtein distance
      (skip[3]?0:(*vit).cls_rank) +    // common longest substring
      (skip[4]?0:(*vit).canon_rank) +  // is it a validated word form
      (skip[5]?0:(*vit).fl_rank) +     // first character equality
      (skip[6]?0:(*vit).ll_rank) +     // last 2 characters equality
      (skip[7]?0:(*vit).khc_rank) +    // known historical confusion
      (skip[8]?0:(*vit).pairs1_rank) + //
      (skip[9]?0:(*vit).pairs2_rank) + //
      (skip[10]?0:(*vit).median_rank) + //
      (skip[11]?0:(*vit).variant_rank) + // # of decapped versions of the CC
      (skip[12]?0:(*vit).cosine_rank) + // WordVector rank
      (skip[13]?0:(*vit).ngram_rank);
    if ( follow ){
      cerr << "Rank=" << rank << endl;
    }
    rank = rank/factor;
    if ( follow ){
      cerr << "Rank/" << factor << " = " << rank << endl;
    }
    sum += rank;
    if ( follow ){
      cerr << "Sum =" << sum << endl;
    }
    (*vit).rank = rank;
    ++vit;
  }

  if ( recs.size() == 1 ){
    recs[0].rank = 1.0;
  }
  else {
    for ( auto& it : recs ){
      if ( follow ){
	cerr << "Rank=(1-" << it.rank << "/" << sum << ") = ";
      }
      it.rank = 1 - it.rank/sum;
      if ( follow ){
	cerr << it.rank << endl;
      }
    }
  }

  // sort records on alphabeticaly on variant and descending on rank
  map<string,multimap<double,record*,std::greater<double>>> output;
  for ( auto& it : recs ){
    string variant = it.variant;
    auto p = output.find(variant);
    if ( p != output.end()  ){
      p->second.insert( make_pair( it.rank, &it) );
      // rank sorted descending per variant.
    }
    else {
      // new variant.
      multimap<double,record*,std::greater<double>> tmp;
      tmp.insert( make_pair( it.rank, &it ) );
      output.insert( make_pair(variant,tmp) );
    }
  }

  // now extract the first 'clip' records for every variant, (best ranked)
  for ( const auto& it : output ){
    const auto& mm = it.second;
    int cnt = 0;
    multimap<double,record,std::greater<double>> tmp;
    for ( const auto& mit : mm ){
      tmp.insert( make_pair( mit.first, *mit.second ) );
      if ( ++cnt >= clip ){
	break;
      }
    }
    // store the result vector
#pragma omp critical (store)
    {
      results.insert( make_pair(it.first,tmp) );
    }
  }

  if ( db ){
    vector<record>::iterator vit = recs.begin();
    multimap<double,string,greater<double>> outv;
    while ( vit != recs.end() ){
      outv.insert( make_pair( (*vit).rank, vit->extractLong(skip) ) );
      ++vit;
    }
#pragma omp critical (debugoutput)
    for ( const auto& oit : outv ){
      *db << oit.second << endl;
    }
  }
}

void collect_ngrams( const vector<record>& records, set<string>& variants_set ){
  vector<record> sorted = records;
  //  cerr << "\nCollecting NEW variant " << records[0].variant << endl;
  // for ( auto const& it : records ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  sort( sorted.begin(), sorted.end(),
	[]( const record& lhs, const record& rhs ){
	  return lhs.ngram_points > rhs.ngram_points;} );
  // cerr << "OUT: " << endl;
  // for ( auto const& it : sorted ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  set<string> variants;
  auto it = sorted.begin();
  while ( it != sorted.end() ){
    if ( verbose ){
#pragma omp critical (log)
      {
	cerr << "NEXT it: " << it->variant << "~" << it->candidate
	     << "::" << it->ngram_points << endl;
      }
    }
    if ( it->ngram_points > 0 ){
      if ( verbose ){
#pragma omp critical (log)
	{
	  cerr << "Remember: " << it->variant << endl;
	}
      }
      variants.insert( it->variant );
    }
    ++it;
  }
#pragma omp critical (update)
  {
    variants_set.insert( variants.begin(), variants.end() );
  }
}

vector<record> filter_ngrams( const vector<record>& records,
			      const set<string>& variants_set ){
  //  cerr << "\nexamining NEW variant " << records[0].variant << endl;
  // for ( auto const& it : records ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  vector<record> result;
  auto it = records.begin();
  while ( it != records.end() ){
    if ( verbose ){
#pragma omp critical (log)
      {
	cerr << "NEXT it: " << it->variant << "~" << it->candidate
	     << "::" << it->ngram_points << endl;
      }
    }
    if ( it->ngram_points == 0 ){
      vector<string> parts = TiCC::split_at(it->variant,SEPARATOR);
      if ( verbose ){
#pragma omp critical (log)
	{
	  cerr << "bekijk variant: " << it->variant << "~" << it->candidate
	       << "::" << it->ngram_points << endl;
	}
      }
      bool forget = false;
      for ( const auto& p: parts ){
	if ( variants_set.find(p) != variants_set.end() ){
	  if ( verbose ){
#pragma omp critical (log)
	    {
	      cerr << "ERASE: " << it->variant << "~" << it->candidate
		   << "::" << it->ngram_points << endl;
	    }
	  }
	  forget = true;
	  break;
	}
      }
      if ( !forget ){
	result.push_back( *it );
      }
    }
    else {
      result.push_back( *it );
    }
    ++it;
  }
  return result;
}

struct wid {
  wid( const string& s, const set<streamsize>& st ): _s(s), _st(st) {};
  string _s;
  set<streamsize> _st;
};

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "alph:,debugfile:,skipcols:,charconf:,charconfreq:,"
			   "artifrq:,subtractartifrqfeature1:,subtractartifrqfeature2:,"
			   "wordvec:,clip:,numvec:,threads:,verbose,follow:,ALTERNATIVE" );
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
  if ( opts.extract('h' ) ){
    usage( progname );
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  bool verbose = opts.extract( 'v' ) || opts.extract("verbose");
  bool ALTERNATIVE = opts.extract( "ALTERNATIVE" );
  string alfabetFile;
  string lexstatFile;
  string freqOutFile;
  string wordvecFile;
  string outFile;
  string debugFile;
  int clip = 0;
  string skipC;
  size_t sub_artifreq = 0;
  size_t sub_artifreq_f1 = 0;
  if ( !opts.extract("charconf",lexstatFile) ){
    cerr << "missing --charconf option" << endl;
    exit(EXIT_FAILURE);
  }
  opts.extract("charconfreq",freqOutFile);
  if ( !opts.extract("alph",alfabetFile) ){
    cerr << "missing --alph option" << endl;
    exit(EXIT_FAILURE);
  }
  opts.extract( "wordvec", wordvecFile );
  opts.extract( 'o', outFile );
  opts.extract( "debugfile", debugFile );
  opts.extract( "skipcols", skipC );
  string value;
  if ( opts.extract( "clip", value ) ){
    if ( !TiCC::stringTo(value,clip) ) {
      cerr << "illegal value for --clip (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract( "subtractartifrqfeature2", value ) ){
    if ( !TiCC::stringTo(value,sub_artifreq) ) {
      cerr << "illegal value for --subtractartifrqfeature2 (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( sub_artifreq == 0 ){
    if ( opts.extract( "artifrq", value ) ){
      cerr << "WARNING: Obsolete option 'artifrq'. Use 'subtractartifrqfeature2' instead." << endl;
      if ( !TiCC::stringTo(value,sub_artifreq) ) {
	cerr << "illegal value for --artifrq (" << value << ")" << endl;
	exit( EXIT_FAILURE );
      }
    }
  }
  if ( opts.extract( "subtractartifrqfeature1", value ) ){
    if ( !TiCC::stringTo(value,sub_artifreq_f1) ) {
      cerr << "illegal value for --subtractartifrqfeature1 (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  int numThreads=1;
  value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
  value = "1";
#ifdef HAVE_OPENMP
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

  while ( opts.extract( "follow", value ) ){
    vector<string> parts = TiCC::split_at( value, "," );
    for ( const auto& p : parts ){
      follow_words.insert( p );
    }
  }

#ifdef HAVE_OPENMP
  if ( !follow_words.empty() && numThreads > 1 ){
    cerr << "FORCING # threads to 1 because of --follow option!" << endl;
    numThreads = 1;
    omp_set_num_threads( numThreads );
  }
#endif

  //#define TESTWV
#ifdef TESTWV
  size_t num_vec = 20;
  if ( opts.extract( "numvec", value ) ){
    if ( !TiCC::stringTo(value,num_vec) ) {
      cerr << "illegal value for --numvec (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
#endif
  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }
  vector<string> fileNames = opts.getMassOpts();
  if ( fileNames.empty() ){
    cerr << "missing an inputfile" << endl;
    exit(EXIT_FAILURE);
  }
  if ( fileNames.size() > 1 ){
    cerr << "only one inputfile may be provided." << endl;
    exit(EXIT_FAILURE);
  }
  string inFile = fileNames[0];
  if ( !TiCC::match_back( inFile, ".ldcalc" ) ){
    cerr << "inputfile must have extension .ldcalc" << endl;
    exit(EXIT_FAILURE);
  }
  if ( !outFile.empty() ){
    if ( !TiCC::match_back( outFile, ".ranked" ) )
      outFile += ".ranked";
  }
  else {
    outFile = inFile + ".ranked";
  }
  if ( outFile == inFile ){
    cerr << "same filename for input and output!" << endl;
    exit(EXIT_FAILURE);
  }
  else if ( outFile == debugFile ){
    cerr << "same filename for output and debug!" << endl;
    exit(EXIT_FAILURE);
  }
  else if ( debugFile == inFile ){
    cerr << "same filename for input and debug!" << endl;
    exit(EXIT_FAILURE);
  }

  ifstream input( inFile );
  if ( !input ){
    cerr << "problem opening confusie file: " << inFile << endl;
    exit(1);
  }

  ifstream lexicon( alfabetFile );
  if ( !lexicon ){
    cerr << "problem opening alfabet file: " << alfabetFile << endl;
    exit(1);
  }

  ifstream lexstats( lexstatFile );
  if ( !lexstats ){
    cerr << "problem opening lexstat file: " << lexstatFile << endl;
    exit(1);
  }

  wordvec_tester WV;
  if ( !wordvecFile.empty() ){
    cerr << "loading word vectors" << endl;
    bool res = WV.fill( wordvecFile );
    if ( !res ){
      cerr << "problem opening wordvec file: " << wordvecFile << endl;
      exit(1);
    }
    cerr << "loaded " << WV.size() << " word vectors" << endl;
#ifdef TESTWV
    vector<word_dist> wv_result;
    if ( !WV.lookup( "dofter", num_vec, wv_result ) ){
      cerr << "faal!" << endl;
    }
    for( size_t i=0; i < wv_result.size(); ++i ){
      cerr << wv_result[i].w << "\t\t" << wv_result[i].d << endl;
    }
    if ( !WV.lookup( "van", num_vec, wv_result ) ){
      cerr << "faal!" << endl;
    }
    for( size_t i=0; i < wv_result.size(); ++i ){
      cerr << wv_result[i].w << "\t\t" << wv_result[i].d << endl;
    }
    if ( !WV.lookup( "dofter van", num_vec, wv_result ) ){
      cerr << "faal!" << endl;
    }
    for( size_t i=0; i < wv_result.size(); ++i ){
      cerr << wv_result[i].w << "\t\t" << wv_result[i].d << endl;
    }
    exit(EXIT_SUCCESS);
#endif
  }


  size_t count=0;
  ofstream os( outFile );
  ofstream *db = 0;
  if ( !debugFile.empty() ){
    db = new ofstream( debugFile );
  }

  set<int> skip_cols;
  vector<bool> skip;
  int skip_factor = RANK_COUNT;
  skip.resize(RANK_COUNT,false);
  if ( !skipC.empty() ){
    vector<string> vec;
    if ( TiCC::split_at( skipC, vec, "," ) == 0 ){
      cerr << "unable te retrieve column numbers from " << skipC << endl;
      exit(EXIT_FAILURE);
    }
    for( size_t i=0; i < vec.size(); ++i ){
      int kol = TiCC::stringTo<int>(vec[i]);
      if ( kol <= 0 || kol > RANK_COUNT ){
	cerr << "invalid skip value in --skipcols. All values must be between 1 and "
	     << RANK_COUNT << endl;
	exit(EXIT_FAILURE);
      }
      skip_cols.insert(kol);
    }
    if ( skip_cols.size() == (unsigned)RANK_COUNT ){
      cerr << "you may not skip all value using --skipcols." << endl;
      exit(EXIT_FAILURE);
    }
    cerr << "skips = " << skip_cols << endl;
    // now stuff it in a bool vector
    for ( int i=0; i < RANK_COUNT; ++i ){
      if ( skip_cols.find(i+1) != skip_cols.end() ){
	skip[i] = true;
	--skip_factor;
      }
    }
  }

  map<UChar,bitType> alfabet;
  cout << "reading alphabet." << endl;
  string line;
  while ( getline( lexicon, line ) ){
    if ( line.size() == 0 || line[0] == '#' )
      continue;
    vector<string> vec;
    if ( TiCC::split( line, vec ) != 3 ){
      cerr << "invalid line '" << line << "' in " << alfabetFile << endl;
      exit( EXIT_FAILURE );
    }
    UnicodeString key = TiCC::UnicodeFromUTF8(vec[0]);
    bitType value = TiCC::stringTo<bitType>( vec[2] );
    alfabet[key[0]] = value;
  }

  map<string,set<streamsize> > fileIds;
  map<bitType,size_t> kwc_counts;
  map<bitType,vector<size_t>> cc_freqs;
  cout << "start indexing input and determining KWC counts AND CC freq per KWC" << endl;
  int failures = 0;
  streamsize pos = input.tellg();
  while ( getline( input, line ) ){
    if ( verbose ){
      cerr << "bekijk " << line << endl;
    }
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "~" ) != RANK_COUNT ){
      cerr << "invalid line: " << line << endl;
      cerr << "expected " << RANK_COUNT << " ~ separated values." << endl;
      if ( ++failures > 50 ){
	cerr << "too many invalid lines" << endl;
	exit(EXIT_FAILURE);
      }
    }
    else {
      string variant = parts[0];
      fileIds[variant].insert( pos );
      bitType kwc = TiCC::stringTo<bitType>(parts[6]);
      ++kwc_counts[kwc];
      size_t ccf = TiCC::stringTo<size_t>(parts[4]);
      cc_freqs[kwc].push_back(ccf);
      if ( ++count % 10000 == 0 ){
	cout << ".";
	cout.flush();
	if ( count % 500000 == 0 ){
	  cout << endl << count << endl;
	}
      }
    }
    pos = input.tellg();
  }
  cout << endl << "Done indexing" << endl;

  map<bitType,size_t> kwc_medians;
  for ( auto& it :  cc_freqs ){
    sort( it.second.begin(), it.second.end() );
    //    cerr << "vector: " << it.second << endl;
    size_t size = it.second.size();
    size_t median =0;
    if ( size %2 == 0 ){
      // even
      median = ( it.second[size/2 -1] + it.second[size/2] ) / 2;
    }
    else {
      median = it.second[size/2];
    }
    //    cerr << "median " << it.first << " = " << median << endl;
    kwc_medians[it.first] = median;
  }
  map<bitType,size_t> kwc2_counts;
  map<bitType,string> kwc_string;

  cout << "reading lexstat file " << lexstatFile
       << " and extracting pairs." << endl;
  while ( getline( lexstats, line ) ){
    vector<string> vec;
    if ( TiCC::split_at( line, vec, "#" ) < 2 ){
      cerr << "invalid line '" << line << "' in " << lexstatFile << endl;
      exit( EXIT_FAILURE );
    }
    bitType key = TiCC::stringTo<bitType>( vec[0] );
    kwc_string[key] = vec[1];
    if ( kwc_counts[key] > 0 ){
      UnicodeString value = TiCC::UnicodeFromUTF8( vec[1] );
      if ( value.length() == 5 && value[2] == '~' ){
	if ( verbose ){
	  cerr << "bekijk tweetal: " << value << " met freq=" << kwc_counts[key]
	       << endl;
	}
	// look up diffs for value[0] - value[3], value[0] - value[4],
	// value[1] - value[3] and value [1] - value[4];
	bitType b1 = alfabet[value[0]];
	bitType b2 = alfabet[value[1]];
	bitType b3 = alfabet[value[3]];
	bitType b4 = alfabet[value[4]];
	if ( b1 != 0 && b2 != 0 && b3 != 0 && b4 != 0 ){
	  vector<size_t> counts(4);
	  bitType diff = abs(b1 - b3);
	  counts[0] = kwc_counts[diff]; // may be 0!
	  diff = abs(b1 - b4);
	  counts[1] = kwc_counts[diff];
	  diff = abs(b2 - b3);
	  counts[2] = kwc_counts[diff];
	  diff = abs(b2 - b4);
	  counts[3] = kwc_counts[diff];
	  size_t max = 0;
	  size_t maxPos = 0;
	  for ( size_t i=0; i < 3; ++i ){
	    if ( counts[i] > max ){
	      max = counts[i];
	      maxPos = i;
	    }
	  }
	  if ( verbose ){
	    cerr << "paarsgewijs max= " << max << " op positie " << maxPos << endl;
	  }
	  if ( max != 0 ){
	    if ( maxPos == 0 ){
	      if ( verbose ){
		cerr << "dus paar = 1-3: " << UnicodeString(value[0]) << "-"
		     << UnicodeString(value[3]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[1]) << "-"
		     << UnicodeString(value[4]) << " met freq=" << counts[3] << endl;
	      }
	      kwc2_counts[key] = max + counts[3];
	    }
	    else if ( maxPos == 1 ){
	      if ( verbose ){
		cerr << "dus paar = 1-4: " << UnicodeString(value[0]) << "-"
		     <<  UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[1]) << "-"
		     << UnicodeString(value[3]) << " met freq=" << counts[2] << endl;
	      }
	      kwc2_counts[key] = max + counts[2];
	    }
	    else if ( maxPos == 2 ){
	      if ( verbose ){
		cerr << "dus paar = 2-3: " << UnicodeString(value[1]) << "-"
		     <<  UnicodeString(value[3]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[0]) << "-"
		    << UnicodeString(value[4]) << " met freq=" << counts[1] << endl;
	      }
	      kwc2_counts[key] = max + counts[1];
	    }
	    else if ( maxPos == 3 ){
	      if ( verbose ){
		cerr << "dus paar = 2-4: " << UnicodeString(value[1]) << "-"
		     <<  UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[0]) << "-"
		    << UnicodeString(value[3]) << " met freq=" << counts[0] << endl;
	      }
	      kwc2_counts[key] = max + counts[0];
	    }
	  }
	}
	else {
	  if ( value.indexOf( "*" ) == -1 &&
	       value.indexOf( "$" ) == -1 ){
	    cerr << "one of the characters from " << value
		 << " isn't in the alphabet" << endl;
	    continue;
	  }
	}
      }
    }
  }

  if ( !freqOutFile.empty() ){
    ofstream os( freqOutFile );
    if ( os.good() ){
      cout << "dumping character confusions into " << freqOutFile << endl;
      multimap<size_t,bitType, std::greater<int> > sorted;
      for ( const auto& it: kwc_counts ){
	if ( it.second != 0 ){
	  // only store non-0 frequencies
	  sorted.insert( make_pair(it.second,it.first) );
	}
      }
      for ( const auto& it : sorted ){
	string tr = kwc_string[it.second];
	if ( tr.empty() ){
	  if ( it.second == 0 ){
	    os << it.second << "\ttransposition\t" << it.first << endl;
	  }
	  else {
	    cerr << "no translation for kwc: " << it.second << endl;
	    os << it.second << "\tmissing\t" << it.first << endl;
	  }
	}
	else {
	  os << it.second << "\t" << tr << "\t" << it.first << endl;
	}
      }
    }
    else {
      cerr << "unable to open " << freqOutFile << endl;
    }
  }

  vector<wid> work;
  for ( const auto& it : fileIds ){
    work.push_back( wid( it.first, it.second ) );
  }
  count = 0;

  cout << "Start searching for ngram proof, with " << work.size()
       << " iterations on " << numThreads << " thread(s)." << endl;
  set<string> variants_set;
#pragma omp parallel for schedule(dynamic,1) shared(variants_set,verbose)
  for( size_t i=0; i < work.size(); ++i ){
    const set<streamsize>& ids = work[i]._st;
    ifstream in( inFile );
    vector<record> records;
    set<streamsize>::const_iterator it = ids.begin();
    while ( it != ids.end() ){
      vector<word_dist> vec;
      in.seekg( *it );
      ++it;
      string line;
      getline( in, line );
      record rec( line, sub_artifreq, sub_artifreq_f1, vec );
      records.push_back( rec );
      if ( verbose ){
	int tmp = 0;
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
	int numt = 0;
#ifdef HAVE_OPENMP
	numt = omp_get_thread_num();
#endif
	if ( numt == 0 && tmp % 10000 == 0 ){
	  cout << ".";
	  cout.flush();
	  if ( tmp % 500000 == 0 ){
	    cout << endl << tmp << endl;
	  }
	}
      }
    }
    collect_ngrams( records, variants_set );
  }

  map<string,multimap<double,record,std::greater<double>>> results;
  cout << "Start the REAL work, with " << work.size()
       << " iterations on " << numThreads << " thread(s)." << endl;
#pragma omp parallel for schedule(dynamic,1) shared(verbose,db)
  for( size_t i=0; i < work.size(); ++i ){
    const set<streamsize>& ids = work[i]._st;
    vector<word_dist> vec;
    if ( WV.size() > 0 ){
      WV.lookup( work[i]._s, 20, vec );
      if ( verbose ){
#pragma omp critical (log)
	{
	  cerr << "looked up: " << work[i]._s << endl;
	}
      }
    }
    ifstream in( inFile );
    vector<record> records;
    set<streamsize>::const_iterator it = ids.begin();
    while ( it != ids.end() ){
      in.seekg( *it );
      ++it;
      string line;
      getline( in, line );
      record rec( line, sub_artifreq, sub_artifreq_f1, vec );
      records.push_back( rec );
      if ( verbose ){
	int tmp = 0;
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
	int numt = 0;
#ifdef HAVE_OPENMP
	numt = omp_get_thread_num();
#endif
	if ( numt == 0 && tmp % 10000 == 0 ){
	  cout << ".";
	  cout.flush();
	  if ( tmp % 500000 == 0 ){
	    cout << endl << tmp << endl;
	  }
	}
      }
    }
    records = filter_ngrams( records, variants_set );
    if ( !records.empty() ){
      if ( ALTERNATIVE ){
	map<bitType,vector<size_t>> local_cc_freqs;
	for ( const auto& it : records ){
	  local_cc_freqs[it.kwc].push_back( it.candidate_freq );
	}
	map<bitType,size_t> local_kwc_medians;
	for ( auto& it : local_cc_freqs ){
	  sort( it.second.begin(), it.second.end() );
	  //    cerr << "vector: " << it.second << endl;
	  size_t size = it.second.size();
	  size_t median =0;
	  if ( size %2 == 0 ){
	    // even
	    median = ( it.second[size/2 -1] + it.second[size/2] ) / 2;
	  }
	  else {
	    median = it.second[size/2];
	  }
	  //    cerr << "median " << it.first << " = " << median << endl;
	  local_kwc_medians[it.first] = median;
	}
	::rank( records, results, clip, kwc_counts, kwc2_counts,
		local_kwc_medians,
		db, skip, skip_factor );
      }
      else {
	::rank( records, results, clip, kwc_counts, kwc2_counts,
		kwc_medians,
		db, skip, skip_factor );
      }
    }
  }

  if ( clip == 1 ){
    // we re-sort the output on descending frequency AND descending on rank,
    // needed for chaining
    // map<string,multimap<double,record,std::greater<double>>> results;
    // but we know that every multimap has only 1 entry for clip = 1
    multimap< size_t, multimap<double, string, std::greater<double>>, std::greater<size_t> > o_vec;
    for ( auto& it : results ){
      const record *rec = &it.second.begin()->second;
      auto oit = o_vec.find( rec->candidate_freq );
      if ( oit != o_vec.end() ){
	oit->second.insert( make_pair( rec->rank, rec->extractResults() ) );
      }
      else {
	multimap<double, string, std::greater<double>> tmp;
	tmp.insert( make_pair( rec->rank, rec->extractResults() ) );
	o_vec.insert( make_pair( rec->candidate_freq, tmp ) );
      }
    }
    // output the results
    for ( const auto& oit : o_vec ){
      for ( const auto& mit : oit.second ){
	os << mit.second << endl;
      }
    }
  }
  else {
    // output the result
    // map<string,multimap<double,record,std::greater<double>>> results;
    for ( const auto& it : results ){
      for( const auto& mit : it.second ){
	os << mit.second.extractResults() << endl;
      }
    }
  }
  cout << "results in " << outFile << endl;
}
