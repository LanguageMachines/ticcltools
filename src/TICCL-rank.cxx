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
#include "config.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"
#include "ticcl/word2vec.h"

using namespace std;
using namespace icu;
using ticcl::bitType;
using TiCC::operator<<;

const int RANK_COUNT=14;

set<UnicodeString> follow_words;

bool verbose = false;

void usage( const string& name ){
  cerr << "usage: " << name << " --alph <alphabetfile> --charconf <lexstat file> [--wordvec <wordvectorfile>] [-o <outputfile>] [-t threads] [--clip <clip>] [--debugfile <debugfile>] [--artifrq art] [--skipcols <skip>] infile" << endl;
  cerr << "\t'infile'\t is a file in TICCL-LDcalc format" << endl;
  cerr << "\t--alph 'alpha'\t an alphabet file in TICCL-lexstat format." << endl;
  cerr << "\t--charconf 'charconfus'\t a character confusion file in TICCL-lexstat format." << endl;
  cerr << "\t--charconfreq 'name'\t Extract a character confusion frequency file" << endl;
  cerr << "\t--wordvec<wordvecfile> read in a google word2vec file." << endl;
  cerr << "\t-o 'outfile'\t name of the output file." << endl;
  cerr << "\t-t <threads> or --threads <threads> Number of threads to run on." << endl;
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

class rank_record {
public:
  rank_record( const UnicodeString &,
	       size_t,
	       size_t,
	       const vector<word_dist>& );
  UnicodeString extractResults() const;
  UnicodeString extractLong( const vector<bool>& skip ) const;
  UnicodeString variant;
  UnicodeString candidate;
  UnicodeString lower_candidate;
  double variant_count;
  double variant_rank;
  size_t variant_freq;
  size_t low_variant_freq;
  size_t candidate_freq;
  size_t reduced_candidate_freq;
  size_t low_candidate_freq;
  double freq_rank;
  bitType char_conf_val;
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
	      const UnicodeString& word ){
  for( size_t i=0; i < vec.size(); ++i ){
    if ( vec[i].w == TiCC::UnicodeToUTF8(word) ){
      //	cerr << "JA! " << vec[i].d << endl;
      return vec[i].d;
    }
  }
  return 0.0;
}

rank_record::rank_record( const UnicodeString& line,
			  size_t sub_artifreq_f1,
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
  vector<UnicodeString> parts = TiCC::split_at( line, "~" );
  // file a rank_record with the RANK_COUNT parts of one line from a LDcalc output file
  if ( parts.size() == RANK_COUNT ){
    variant = parts[0];
    variant_freq = TiCC::stringTo<size_t>(parts[1]);
    low_variant_freq = TiCC::stringTo<size_t>(parts[2]);
    candidate = parts[3];
    UnicodeString ls = candidate;
    ls.toLower();
    lower_candidate = ls;
    variant_rank = -2000;  // bogus value, is set later
    UnicodeString f2_string = parts[4];
    candidate_freq = TiCC::stringTo<size_t>( f2_string );
    reduced_candidate_freq = candidate_freq;
    if ( sub_artifreq_f1 > 0 && reduced_candidate_freq >= sub_artifreq_f1 ){
      reduced_candidate_freq -= sub_artifreq_f1;
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
    char_conf_val   = TiCC::stringTo<bitType>(parts[6]);
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

UnicodeString rank_record::extractLong( const vector<bool>& skip ) const {
  ostringstream ss;
  ss << variant << "#"
     << variant_freq << "#"
     << low_variant_freq << "#"
     << candidate << "#"
     << candidate_freq << "#"
     << low_candidate_freq << "#"
     << char_conf_val << "#"
     << f2len << "~";
  double the_rank = 0;
  if ( skip[0] ){
    ss << "N#";
  }
  else {
    the_rank += f2len_rank;
    ss << f2len_rank << "#";
  }
  ss << reduced_candidate_freq << "~";
  if ( skip[1] ){
    ss << "N#";
  }
  else {
    the_rank += freq_rank;
    ss << freq_rank << "#";
  }
  ss << ld << "~";
  if ( skip[2] ){
    ss << "N#";
  }
  else {
    the_rank += ld_rank;
    ss << ld_rank << "#";
  }
  ss << cls << "~";
  if ( skip[3] ){
    ss << "N#";
  }
  else {
    the_rank += cls_rank;
    ss << cls_rank << "#";
  }
  ss << canon << "~";
  if ( skip[4] ){
    ss << "N#";
  }
  else {
    the_rank += canon_rank;
    ss << canon_rank << "#";
  }
  ss << fl << "~";
  if ( skip[5] ){
    ss << "N#";
  }
  else {
    the_rank += fl_rank;
    ss << fl_rank << "#";
  }
  ss << ll << "~";
  if ( skip[6] ){
    ss << "N#";
  }
  else {
    the_rank += ll_rank;
    ss << ll_rank << "#";
  }
  ss << khc << "~";
  if ( skip[7] ){
    ss << "N#";
  }
  else {
    the_rank += khc_rank;
    ss << khc_rank << "#";
  }
  ss << pairs1 << "~";
  if ( skip[8] ){
    ss << "N#";
  }
  else {
    the_rank += pairs1_rank;
    ss << pairs1_rank << "#";
  }
  ss <<pairs2 << "~";
  if ( skip[9] ){
    ss << "N#";
  }
  else {
    the_rank += pairs2_rank;
    ss << pairs2_rank << "#";
  }
  ss << median << "~";
  if ( skip[10] ){
    ss << "N#";
  }
  else {
    the_rank += median_rank;
    ss << median_rank << "#";
  }
  ss << variant_count << "~";
  if ( skip[11] ){
    ss << "N#";
  }
  else {
    the_rank += variant_rank;
    ss << variant_rank << "#";
  }
  ss <<cosine << "~";
  if ( skip[12] ){
    ss << "N#";
  }
  else {
    the_rank += cosine_rank;
    ss << cosine_rank << "#";
  }
  ss << ngram_points << "~";
  if ( skip[13] ){
    ss << "N#";
  }
  else {
    the_rank += ngram_rank;
    ss << ngram_rank << "#";
  }
  ss << the_rank << "#"
     << rank;
  return TiCC::UnicodeFromUTF8(ss.str());
}

UnicodeString rank_record::extractResults() const {
  ostringstream ss;
  ss << variant << "#"
     << variant_freq << "#"
     << candidate << "#"
     << candidate_freq << "#"
     << char_conf_val << "#"
     << ld << "#"
     << rank;
  return TiCC::UnicodeFromUTF8(ss.str());
}

template <typename TObject, typename TMember, typename TValue>
void set_val( TObject& object, TMember member, TValue value )
{
    ( object ).*member = value;
}

template< class Tmap, typename TMember >
void rank_desc_map( const Tmap& desc_map,
		    vector<rank_record>& recs,
		    TMember member ){
  // the map is a (multi-)map, sorted descending on the first value
  // whichs is an integer value.
  // the second value of the map is an index in the vector of rank_records.
  // so the map contains rank_records, sorted descending on specific values
  // e.g:
  //  clsmap holds the CommonSubstringLengths of all vectors, longest first
  if ( desc_map.empty() ){
    return;
  }
  size_t last = desc_map.begin()->first; // start with the longest
  int ranking = 1;                    // it will be ranked 1
  for ( const auto& [val,index] : desc_map ){
    if ( val < last ){
      // we find a shorter. so ranking is incremented (meaning LOWER ranking)
      last = val;
      ++ranking;
    }
    // set the currect ranking for the rank_record at hand
    set_val( recs[index], member, ranking );
  }
}

void rank( vector<rank_record>& recs,
	   map<UnicodeString,multimap<double,rank_record,std::greater<double>>>& results,
	   int clip,
	   const map<bitType,size_t>& char_conf_val_counts,
	   const map<bitType,size_t>& char_conf_val2_counts,
	   const map<bitType,size_t>& char_conf_val_medians,
	   ostream* db, const vector<bool>& skip, int factor ){
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
  map<UnicodeString,int> lowvarmap;
  size_t count = 0;

  for ( auto& it : recs ){
    // for every rank_record, we store information in descending multimaps
    // So in freqmap, the (index of) the rank_records with highest freq
    //   are stored in front
    //same for f2len, ld, cls and ngram points
    f2lenmap.insert( make_pair(it.f2len, count ) ); // f2lengths descending
    freqmap.insert( make_pair(it.reduced_candidate_freq, count ) ); // freqs descending
    ldmap.insert( make_pair(it.ld,count) ); // lds sorted ASCENDING
    clsmap.insert( make_pair(it.cls,count) ); // cls sorted descending
    size_t var1_cnt = char_conf_val_counts.at(it.char_conf_val);
    it.pairs1 = var1_cnt;
    pairmap1.insert( make_pair(var1_cnt,count )); // #variants descending
    size_t var2_cnt = 0;
    try {
      var2_cnt += char_conf_val2_counts.at(it.char_conf_val);
    }
    catch(...){
    }
    it.pairs2 = var2_cnt;
    pairmap2.insert( make_pair(var2_cnt,count )); // #variants decending
    it.median = char_conf_val_medians.at(it.char_conf_val);
    median_map.insert( make_pair(it.median,count )); // #medians decending
    ngram_map.insert( make_pair(it.ngram_points,count) ); // ngrampoints sorted descending
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
    cout << "8 pairmap1 = " << pairmap1 << endl;
    cout << "9 pairmap2 = " << pairmap2 << endl;
    cout << "10 medianmap = " << median_map << endl;
    //    cout << "12-a lowvarmap = " << lowvarmap << endl;
    cout << "11 lower_variantmap = " << lower_variantmap << endl;
    cout << "14 ngram_map = " << ngram_map << endl;
  }
  rank_desc_map( f2lenmap, recs, &rank_record::f2len_rank );
  if ( follow ){
    cout << "step 1: f2len_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.f2len_rank << endl;
    }
  }

  rank_desc_map( freqmap, recs, &rank_record::freq_rank );
  if ( follow ){
    cout << "step 2: freq_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.freq_rank << endl;
    }
  }

  if ( !ldmap.empty() ){
    int ranking = 1;
    size_t last = ldmap.begin()->first;
    for ( const auto& [val,index] : ldmap ){
      if ( val > last ){
   	last = val;
   	++ranking;
      }
      recs[index].ld_rank = ranking;
    }
  }
  if ( follow ){
    cout << "step 3: ld_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ld_rank << endl;
    }
  }

  rank_desc_map( clsmap, recs, &rank_record::cls_rank );
  if ( follow ){
    cout << "step 4: cls_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.cls_rank << endl;
    }
    cout << "step 5: canon_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.canon_rank << endl;
    }
    cout << "step 6: fl_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.fl_rank << endl;
    }
    cout << "step 7: ll_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ll_rank << endl;
    }
    cout << "step 8: khc_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.khc_rank << endl;
    }
  }

  rank_desc_map( pairmap1, recs, &rank_record::pairs1_rank );
  if ( follow ){
    cout << "step 9: pairs1_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.pairs1_rank << endl;
    }
  }

  rank_desc_map( pairmap2, recs, &rank_record::pairs2_rank );
  if ( follow ){
    cout << "step 10: pairs2_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.pairs2_rank << endl;
    }
  }

  rank_desc_map( median_map, recs, &rank_record::median_rank );
  if ( follow ){
    cout << "step 11: median_rank: for " << recs.begin()->variant << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.median_rank << endl;
    }
  }

  if ( !lower_variantmap.empty() ){
    int ranking = 1;
    int last = lower_variantmap.begin()->first;
    for ( const auto& [cnt,index] : lower_variantmap ){
      if ( cnt < last ){
	last = cnt;
	++ranking;
      }
      recs[index].variant_count = cnt;
      recs[index].variant_rank = ranking;
    }
  }
  if ( follow ){
    cout << "step 11: lower_variant_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " count=" << r.variant_count << " rank= " << r.variant_rank << endl;
    }
    cout << "step 13: cosine_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.cosine_rank << endl;
    }
  }

  rank_desc_map( ngram_map, recs, &rank_record::ngram_rank );
  if ( follow ){
    cout << "step 14: ngram_rank: " << endl;
    for ( const auto& r : recs ){
      cout << "\t" << r.candidate << " rank= " << r.ngram_rank << endl;
    }
  }

  double sum = 0.0;
  for( auto& vit : recs ){
    double rank_v =
      (skip[0]?0:vit.f2len_rank) +  // number of characters in the frequency
      (skip[1]?0:vit.freq_rank) +   // frequency of the CC
      (skip[2]?0:vit.ld_rank) +     // levenshtein distance
      (skip[3]?0:vit.cls_rank) +    // common longest substring
      (skip[4]?0:vit.canon_rank) +  // is it a validated word form
      (skip[5]?0:vit.fl_rank) +     // first character equality
      (skip[6]?0:vit.ll_rank) +     // last 2 characters equality
      (skip[7]?0:vit.khc_rank) +    // known historical confusion
      (skip[8]?0:vit.pairs1_rank) + //
      (skip[9]?0:vit.pairs2_rank) + //
      (skip[10]?0:vit.median_rank) + //
      (skip[11]?0:vit.variant_rank) + // # of decapped versions of the CC
      (skip[12]?0:vit.cosine_rank) + // WordVector rank
      (skip[13]?0:vit.ngram_rank);
    if ( follow ){
      cerr << "Rank=" << rank_v << endl;
    }
    rank_v = rank_v/factor;
    if ( follow ){
      cerr << "Rank/" << factor << " = " << rank_v << endl;
    }
    sum += rank_v;
    if ( follow ){
      cerr << "Sum =" << sum << endl;
    }
    vit.rank = rank_v;
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

  // sort rank_records on alphabeticaly on variant and descending on rank
  map<UnicodeString,multimap<double,rank_record*,std::greater<double>>> output;
  for ( auto& it : recs ){
    UnicodeString variant = it.variant;
    auto p = output.find(variant);
    if ( p != output.end()  ){
      p->second.insert( make_pair( it.rank, &it) );
      // rank sorted descending per variant.
    }
    else {
      // new variant.
      multimap<double,rank_record*,std::greater<double>> tmp;
      tmp.insert( make_pair( it.rank, &it ) );
      output.insert( make_pair(variant,tmp) );
    }
  }

  // now extract the first 'clip' rank_records for every variant, (best ranked)
  for ( const auto& [word,rec_map] : output ){
    int cnt = 0;
    multimap<double,rank_record,std::greater<double>> tmp;
    for ( const auto& [val,rec] : rec_map ){
      tmp.insert( make_pair( val, *rec ) );
      if ( ++cnt >= clip ){
	break;
      }
    }
    // store the result vector
#pragma omp critical (store)
    {
      results.insert( make_pair(word,tmp) );
    }
  }

  if ( db ){
    multimap<double,UnicodeString,greater<double>> outv;
    for ( const auto& vit : recs ){
      outv.insert( make_pair( vit.rank, vit.extractLong(skip) ) );
    }
#pragma omp critical (debugoutput)
    for ( const auto& oit : outv ){
      *db << oit.second << endl;
    }
  }
}

void collect_ngrams( const vector<rank_record>& rank_records,
		     set<UnicodeString>& variants_set ){
  vector<rank_record> sorted = rank_records;
  //  cerr << "\nCollecting NEW variant " << rank_records[0].variant << endl;
  // for ( auto const& it : rank_records ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  sort( sorted.begin(), sorted.end(),
	[]( const rank_record& lhs, const rank_record& rhs ){
	  return lhs.ngram_points > rhs.ngram_points;} );
  // cerr << "OUT: " << endl;
  // for ( auto const& it : sorted ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  set<UnicodeString> variants;
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

vector<rank_record> filter_ngrams( const vector<rank_record>& rank_records,
				   const set<UnicodeString>& variants_set ){
  //  cerr << "\nexamining NEW variant " << rank_records[0].variant << endl;
  // for ( auto const& it : rank_records ){
  //   cerr << it.variant << "~" << it.candidate << "::" << it.ngram_points << endl;
  // }
  // cerr << endl;
  vector<rank_record> result;
  auto it = rank_records.begin();
  while ( it != rank_records.end() ){
    if ( verbose ){
#pragma omp critical (log)
      {
	cerr << "NEXT it: " << it->variant << "~" << it->candidate
	     << "::" << it->ngram_points << endl;
      }
    }
    if ( it->ngram_points == 0 ){
      vector<UnicodeString> parts = TiCC::split_at(it->variant,
						   ticcl::US_SEPARATOR);
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
  wid( const UnicodeString& s, const set<streamsize>& st ): _s(s), _st(st) {};
  UnicodeString _s;
  set<streamsize> _st;
};

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.add_short_options( "vVho:t:" );
    opts.add_long_options( "alph:,debugfile:,skipcols:,charconf:,charconfreq:,"
			   "artifrq:,"
			   "subtractartifrqfeature1:,subtractartifrqfeature2:,"
			   "wordvec:,clip:,numvec:,threads:,verbose,follow:,"
			   "help,version,ALTERNATIVE" );
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
  verbose = opts.extract( 'v' ) || opts.extract("verbose");
  bool ALTERNATIVE = opts.extract( "ALTERNATIVE" );
  string alphabetFile;
  string lexstatFile;
  string freqOutFile;
  string wordvecFile;
  string outFile;
  string debugFile;
  int clip = 0;
  string skipC;
  size_t sub_artifreq_f1 = 0;
  size_t sub_artifreq_f2 = 0;
  if ( !opts.extract("charconf",lexstatFile) ){
    cerr << "missing --charconf option" << endl;
    exit(EXIT_FAILURE);
  }
  opts.extract("charconfreq",freqOutFile);
  if ( !opts.extract("alph",alphabetFile) ){
    cerr << "missing --alph option" << endl;
    exit(EXIT_FAILURE);
  }
  opts.extract( "wordvec", wordvecFile );
  opts.extract( 'o', outFile );
  opts.extract( "debugfile", debugFile );
  opts.extract( "skipcols", skipC );
  string arg_val;
  if ( opts.extract( "clip", arg_val ) ){
    if ( !TiCC::stringTo(arg_val,clip) ) {
      cerr << "illegal value for --clip (" << arg_val << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( opts.extract( "subtractartifrqfeature2", arg_val ) ){
    if ( !TiCC::stringTo(arg_val,sub_artifreq_f2) ) {
      cerr << "illegal value for --subtractartifrqfeature2 (" << arg_val
	   << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  if ( sub_artifreq_f2 == 0 ){
    if ( opts.extract( "artifrq", arg_val ) ){
      cerr << "WARNING: Obsolete option 'artifrq'. Use 'subtractartifrqfeature2' instead." << endl;
      if ( !TiCC::stringTo( arg_val, sub_artifreq_f2) ) {
	cerr << "illegal value for --artifrq (" << arg_val << ")" << endl;
	exit( EXIT_FAILURE );
      }
    }
  }
  if ( opts.extract( "subtractartifrqfeature1", arg_val ) ){
    if ( !TiCC::stringTo(arg_val,sub_artifreq_f1) ) {
      cerr << "illegal value for --subtractartifrqfeature1 (" << arg_val
	   << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  int numThreads=1;
  arg_val = "1";
  if ( !opts.extract( 't', arg_val ) ){
    opts.extract( "threads", arg_val );
  }
#ifdef HAVE_OPENMP
  if ( TiCC::lowercase(arg_val) == "max" ){
    numThreads = omp_get_max_threads() - 2;
    omp_set_num_threads( numThreads );
    cout << "running on " << numThreads << " threads." << endl;
  }
  else {
    if ( !TiCC::stringTo(arg_val,numThreads) ) {
      cerr << "illegal value for -t (" << arg_val << ")" << endl;
      exit( EXIT_FAILURE );
    }
    omp_set_num_threads( numThreads );
    cout << "running on " << numThreads << " threads." << endl;
  }
#else
  if ( arg_val != "1" ){
    cerr << "unable to set number of threads!.\nNo OpenMP support available!"
	 <<endl;
    exit(EXIT_FAILURE);
  }
#endif

  while ( opts.extract( "follow", arg_val ) ){
    vector<UnicodeString> parts = TiCC::split_at( TiCC::UnicodeFromUTF8(arg_val), "," );
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
  if ( opts.extract( "numvec", arg_val ) ){
    if ( !TiCC::stringTo(value,num_vec) ) {
      cerr << "illegal value for --numvec (" << arg_val << ")" << endl;
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

  if ( !TiCC::isFile( inFile ) ){
    cerr << "problem opening confusie file: " << inFile << endl;
    exit(1);
  }

  if ( !TiCC::isFile( alphabetFile ) ){
    cerr << "problem opening alphabet file: " << alphabetFile << endl;
    exit(1);
  }

  if ( !TiCC::isFile( lexstatFile ) ){
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
    vector<string> vec = TiCC::split_at( skipC,  "," );
    if ( vec.empty() ){
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

  cout << "reading alphabet." << endl;
  map<UChar,bitType> alphabet;
  ifstream is( alphabetFile );
  ticcl::fillAlphabet( is, alphabet );
  map<UnicodeString,set<streamsize> > fileIds;
  map<bitType,size_t> char_conf_val_counts;
  map<bitType,vector<size_t>> cc_freqs;
  cout << "start indexing input and determining CHAR_CONF_VAL counts AND CC freq per CHAR_CONF_VAL" << endl;
  int failures = 0;
  ifstream input( inFile );
  streamsize pos = input.tellg();
  UnicodeString input_line;
  while ( TiCC::getline( input, input_line ) ){
    if ( verbose ){
      cerr << "bekijk " << input_line << endl;
    }
    vector<UnicodeString> parts = TiCC::split_at( input_line, "~" );
    if ( parts.size() != RANK_COUNT ){
      cerr << "invalid line: " << input_line << endl;
      cerr << "expected " << RANK_COUNT << " ~ separated values." << endl;
      if ( ++failures > 50 ){
	cerr << "too many invalid lines" << endl;
	exit(EXIT_FAILURE);
      }
    }
    else {
      UnicodeString variant = parts[0];
      fileIds[variant].insert( pos );
      bitType char_conf_val = TiCC::stringTo<bitType>(parts[6]);
      ++char_conf_val_counts[char_conf_val];
      size_t ccf = TiCC::stringTo<size_t>(parts[4]);
      cc_freqs[char_conf_val].push_back(ccf);
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

  map<bitType,size_t> char_conf_val_medians;
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
    char_conf_val_medians[it.first] = median;
  }
  map<bitType,size_t> char_conf_val2_counts;
  map<bitType,UnicodeString> char_conf_val_string;

  cout << "reading lexstat file " << lexstatFile
       << " and extracting pairs." << endl;
  ifstream lexstats( lexstatFile );
  UnicodeString stats_line;
  while ( TiCC::getline( lexstats, stats_line ) ){
    vector<UnicodeString> vec = TiCC::split_at( stats_line, "#" );
    if ( vec.size() < 2 ){
      cerr << "invalid line '" << stats_line << "' in " << lexstatFile << endl;
      exit( EXIT_FAILURE );
    }
    bitType key = TiCC::stringTo<bitType>( vec[0] );
    char_conf_val_string[key] = vec[1];
    if ( char_conf_val_counts[key] > 0 ){
      UnicodeString value = vec[1];
      if ( value.length() == 5 && value[2] == '~' ){
	if ( verbose ){
	  cerr << "bekijk tweetal: " << value << " met freq=" << char_conf_val_counts[key]
	       << endl;
	}
	// look up diffs for value[0] - value[3], value[0] - value[4],
	// value[1] - value[3] and value [1] - value[4];
	bitType b1 = alphabet[value[0]];
	bitType b2 = alphabet[value[1]];
	bitType b3 = alphabet[value[3]];
	bitType b4 = alphabet[value[4]];
	if ( b1 != 0 && b2 != 0 && b3 != 0 && b4 != 0 ){
	  vector<size_t> counts(4);
	  bitType diff;
	  if ( b1 > b3 ){
	    diff = b1 - b3;
	  }
	  else {
	    diff = b3 - b1;
	  }
	  counts[0] = char_conf_val_counts[diff]; // may be 0!
	  if ( b1 > b4 ){
	    diff = b1 - b4;
	  }
	  else {
	    diff = b4 - b1;
	  }
	  counts[1] = char_conf_val_counts[diff];
	  if ( b2 > b3 ){
	    diff = b2 - b3;
	  }
	  else {
	    diff = b3 - b2;
	  }
	  counts[2] = char_conf_val_counts[diff];
	  if ( b2 > b4 ){
	    diff = b2 - b4;
	  }
	  else {
	    diff = b4 - b2;
	  }
	  counts[3] = char_conf_val_counts[diff];
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
	      char_conf_val2_counts[key] = max + counts[3];
	    }
	    else if ( maxPos == 1 ){
	      if ( verbose ){
		cerr << "dus paar = 1-4: " << UnicodeString(value[0]) << "-"
		     <<  UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[1]) << "-"
		     << UnicodeString(value[3]) << " met freq=" << counts[2] << endl;
	      }
	      char_conf_val2_counts[key] = max + counts[2];
	    }
	    else if ( maxPos == 2 ){
	      if ( verbose ){
		cerr << "dus paar = 2-3: " << UnicodeString(value[1]) << "-"
		     <<  UnicodeString(value[3]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[0]) << "-"
		    << UnicodeString(value[4]) << " met freq=" << counts[1] << endl;
	      }
	      char_conf_val2_counts[key] = max + counts[1];
	    }
	    else if ( maxPos == 3 ){
	      if ( verbose ){
		cerr << "dus paar = 2-4: " << UnicodeString(value[1]) << "-"
		     <<  UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << UnicodeString(value[0]) << "-"
		    << UnicodeString(value[3]) << " met freq=" << counts[0] << endl;
	      }
	      char_conf_val2_counts[key] = max + counts[0];
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
    ofstream fs( freqOutFile );
    if ( fs.good() ){
      cout << "dumping character confusions into " << freqOutFile << endl;
      multimap<size_t,bitType, std::greater<int> > sorted;
      for ( const auto& [val,freq]: char_conf_val_counts ){
	if ( freq != 0 ){
	  // only store non-0 frequencies
	  sorted.insert( make_pair(freq,val) );
	}
      }
      for ( const auto& [val,index] : sorted ){
	UnicodeString tr = char_conf_val_string[index];
	if ( tr.isEmpty() ){
	  if ( index == 0 ){
	    fs << index << "\ttransposition\t" << val << endl;
	  }
	  else {
	    cerr << "no translation for char_conf_val: " << index << endl;
	    fs << index << "\tmissing\t" << val << endl;
	  }
	}
	else {
	  fs << index << "\t" << tr << "\t" << val << endl;
	}
      }
    }
    else {
      cerr << "unable to open " << freqOutFile << endl;
    }
  }

  vector<wid> work;
  for ( const auto& [file,id] : fileIds ){
    work.push_back( wid( file, id ) );
  }
  count = 0;

  cout << "Start searching for ngram proof, with " << work.size()
       << " iterations on " << numThreads << " thread(s)." << endl;
  set<UnicodeString> variants_set;
#pragma omp parallel for schedule(dynamic,1) shared(variants_set,verbose)
  for( size_t i=0; i < work.size(); ++i ){
    const set<streamsize>& ids = work[i]._st;
    ifstream in( inFile );
    vector<rank_record> rank_records;
    auto id_iter = ids.begin();
    while ( id_iter != ids.end() ){
      vector<word_dist> vec;
      in.seekg( *id_iter );
      ++id_iter;
      UnicodeString rec_line;
      TiCC::getline( in, rec_line );
      rank_record rec( rec_line, sub_artifreq_f1, sub_artifreq_f2, vec );
      rank_records.push_back( rec );
      if ( verbose ){
	int tmp = 0;
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
#ifdef HAVE_OPENMP
	int numt = omp_get_thread_num();
	if ( numt == 0 ){
#endif
	  if ( tmp % 10000 == 0 ){
	    cout << ".";
	    cout.flush();
	    if ( tmp % 500000 == 0 ){
	      cout << endl << tmp << endl;
	    }
	  }
#ifdef HAVE_OPENMP
	}
#endif
      }
    }
    collect_ngrams( rank_records, variants_set );
  }

  map<UnicodeString,multimap<double,rank_record,std::greater<double>>> results;
  cout << "Start the REAL work, with " << work.size()
       << " iterations on " << numThreads << " thread(s)." << endl;
#pragma omp parallel for schedule(dynamic,1) shared(verbose,db)
  for( size_t i=0; i < work.size(); ++i ){
    const set<streamsize>& ids = work[i]._st;
    vector<word_dist> vec;
    if ( WV.size() > 0 ){
      WV.lookup( TiCC::UnicodeToUTF8(work[i]._s), 20, vec );
      if ( verbose ){
#pragma omp critical (log)
	{
	  cerr << "looked up: " << work[i]._s << endl;
	}
      }
    }
    ifstream in( inFile );
    vector<rank_record> rank_records;
    auto id_iter = ids.begin();
    while ( id_iter != ids.end() ){
      in.seekg( *id_iter );
      ++id_iter;
      UnicodeString line;
      TiCC::getline( in, line );
      rank_record rec( line, sub_artifreq_f1, sub_artifreq_f2, vec );
      rank_records.push_back( rec );
      if ( verbose ){
	int tmp = 0;
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
#ifdef HAVE_OPENMP
	int numt = omp_get_thread_num();
	if ( numt == 0 ){
#endif
	  if ( tmp % 10000 == 0 ){
	    cout << ".";
	    cout.flush();
	    if ( tmp % 500000 == 0 ){
	      cout << endl << tmp << endl;
	    }
	  }
#ifdef HAVE_OPENMP
	}
#endif
      }
    }
    rank_records = filter_ngrams( rank_records, variants_set );
    if ( !rank_records.empty() ){
      if ( ALTERNATIVE ){
	map<bitType,vector<size_t>> local_cc_freqs;
	for ( const auto& r_it : rank_records ){
	  local_cc_freqs[r_it.char_conf_val].push_back( r_it.candidate_freq );
	}
	map<bitType,size_t> local_char_conf_val_medians;
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
	  local_char_conf_val_medians[it.first] = median;
	}
	::rank( rank_records, results, clip, char_conf_val_counts, char_conf_val2_counts,
		local_char_conf_val_medians,
		db, skip, skip_factor );
      }
      else {
	::rank( rank_records, results, clip, char_conf_val_counts, char_conf_val2_counts,
		char_conf_val_medians,
		db, skip, skip_factor );
      }
    }
  }

  if ( clip == 1 ){
    // we re-sort the output on descending frequency AND descending on rank,
    // needed for chaining
    // map<string,multimap<double,rank_record,std::greater<double>>> results;
    // but we know that every multimap has only 1 entry for clip = 1
    multimap< size_t, multimap<double, UnicodeString, std::greater<double>>, std::greater<size_t> > o_vec;
    for ( const auto& it : results ){
      const rank_record *rec = &it.second.begin()->second;
      auto oit = o_vec.find( rec->candidate_freq );
      if ( oit != o_vec.end() ){
	oit->second.insert( make_pair( rec->rank, rec->extractResults() ) );
      }
      else {
	multimap<double, UnicodeString, std::greater<double>> tmp;
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
    for ( const auto& it : results ){
      for( const auto& mit : it.second ){
	os << mit.second.extractResults() << endl;
      }
    }
  }
  cout << "results in " << outFile << endl;
}
