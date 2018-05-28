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
typedef signed long int bitType;

const int RANK_COUNT=14;

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
  cerr << "\t--artifrq 'arti'\t decrease frequencies with value 'arti'." << endl;
  cerr << "\t\t\t (which should match the artifreq used in TICCL-LDcalc)" << endl;
  cerr << "\t--skipcols=arglist\t skip the named columns in the ranking." << endl;
  cerr << "\t\t\t e.g. if arglist=3,9, then the columns 3 and 9 are not used." << endl;
  cerr << "\t-v\t\t run (very) verbose" << endl;
  exit( EXIT_FAILURE );
}

class record {
public:
  record( const string &, size_t, const vector<word_dist>& );
  string variant1;
  string variant2;
  string lowervariant2;
  double variant_count;
  double variant_rank;
  size_t freq1;
  size_t low_freq1;
  size_t freq2;
  size_t reduced_freq2;
  size_t low_freq2;
  double freq_rank;
  bitType kwc;
  size_t f2len;
  size_t f2len_rank;
  int ld;
  double ld_rank;
  int csl;
  double csl_rank;
  int canon;
  double canon_rank;
  size_t pairs1;
  double pairs1_rank;
  size_t pairs2;
  double pairs2_rank;
  size_t pairs_combined;
  double pairs_combined_rank;
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
		size_t artifreq,
		const vector<word_dist>& WV ):
  variant_count(-1),
  f2len_rank(-1),
  ld(-1),
  pairs1(0),
  pairs1_rank(-1),
  pairs2(0),
  pairs2_rank(-1),
  pairs_combined(0),
  pairs_combined_rank(-1),
  rank(-10000)
{
  vector<string> parts = TiCC::split_at( line, "~" );
  if ( parts.size() == 14 ){
    variant1 = parts[0];
    freq1 = TiCC::stringTo<size_t>(parts[1]);
    low_freq1 = TiCC::stringTo<size_t>(parts[2]);
    variant2 = parts[3];
    icu::UnicodeString us = TiCC::UnicodeFromUTF8( variant2 );
    us.toLower();
    lowervariant2 = TiCC::UnicodeToUTF8( us );
    variant_rank = -2000;
    freq2 = TiCC::stringTo<size_t>(parts[4]);
    f2len = parts[4].length();
    reduced_freq2 = freq2;
    if ( artifreq > 0 && reduced_freq2 >= artifreq ){
      reduced_freq2 -= artifreq;
    }
    low_freq2 = TiCC::stringTo<size_t>(parts[5]);
    freq_rank = -20;
    kwc   = TiCC::stringTo<bitType>(parts[6]);
    ld    = TiCC::stringTo<int>(parts[7]);
    ld_rank = -4.5;
    csl   = TiCC::stringTo<int>(parts[8]);
    csl_rank = -5.6;
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
    ngram_rank = -6.7;
    cosine = lookup( WV, variant2 );
    if ( cosine <= 0.001 )
      cosine_rank = 1;
    else
      cosine_rank = 10;
  }
}

string extractLong( const record& rec, const vector<bool>& skip ){
  string result = rec.variant1 + "#";
  result += TiCC::toString(rec.freq1) + "#";
  result += TiCC::toString(rec.low_freq1) + "#";
  result += rec.variant2 + "#";
  result += TiCC::toString(rec.freq2) + "#";
  result += TiCC::toString(rec.low_freq2) + "#";
  result += TiCC::toString(rec.kwc) + "#";
  result += TiCC::toString(rec.f2len) + "~";
  double rank = 0;
  if ( skip[0] ){
    result += "N#";  }
  else {
    rank += rec.f2len_rank;
    result += TiCC::toString(rec.f2len_rank) + "#";
  }
  result += TiCC::toString(rec.reduced_freq2) + "~";
  if ( skip[1] ){
    result += "N#";
  }
  else {
    rank += rec.freq_rank;
    result += TiCC::toString(rec.freq_rank) + "#";
  }
  result += TiCC::toString(rec.ld) + "~";
  if ( skip[2] ){
    result += "N#";
  }
  else {
    rank += rec.ld_rank;
    result += TiCC::toString(rec.ld_rank) + "#";
  }
  result += TiCC::toString(rec.csl) + "~";
  if ( skip[3] ){
    result += "N#";
  }
  else {
    rank += rec.csl_rank;
    result += TiCC::toString(rec.csl_rank) + "#";
  }
  result += TiCC::toString(rec.canon) + "~";
  if ( skip[4] ){
    result += "N#";
  }
  else {
    rank += rec.canon_rank;
    result += TiCC::toString(rec.canon_rank) + "#";
  }
  result += TiCC::toString(rec.fl) + "~";
  if ( skip[5] ){
    result += "N#";
  }
  else {
    rank += rec.fl_rank;
    result += TiCC::toString(rec.fl_rank) + "#";
  }
  result += TiCC::toString(rec.ll) + "~";
  if ( skip[6] ){
    result += "N#";
  }
  else {
    rank += rec.ll_rank;
    result += TiCC::toString(rec.ll_rank) + "#";
  }
  result += TiCC::toString(rec.khc) + "~";
  if ( skip[7] ){
    result += "N#";
  }
  else {
    rank += rec.khc_rank;
    result += TiCC::toString(rec.khc_rank) + "#";
  }
  result += TiCC::toString(rec.pairs1) + "~";
  if ( skip[8] ){
    result += "N#";
  }
  else {
    rank += rec.pairs1_rank;
    result += TiCC::toString(rec.pairs1_rank) + "#";
  }
  result += TiCC::toString(rec.pairs2) + "~";
  if ( skip[9] ){
    result += "N#";
  }
  else {
    rank += rec.pairs2_rank;
    result += TiCC::toString(rec.pairs2_rank) + "#";
  }
  result += TiCC::toString(rec.pairs_combined) + "~";
  if ( skip[10] ){
    result += "N#";
  }
  else {
    rank += rec.pairs_combined_rank;
    result += TiCC::toString(rec.pairs_combined_rank) + "#";
  }
  result += TiCC::toString(rec.variant_count) + "~";
  if ( skip[11] ){
    result += "N#";
  }
  else {
    rank += rec.variant_rank;
    result += TiCC::toString(rec.variant_rank) + "#";
  }
  result += TiCC::toString(rec.cosine) + "~";
  if ( skip[12] ){
    result += "N#";
  }
  else {
    rank += rec.cosine_rank;
    result += TiCC::toString(rec.cosine_rank) + "#";
  }
  result += TiCC::toString(rec.ngram_points) + "~";
  if ( skip[13] ){
    result += "N#";
  }
  else {
    rank += rec.ngram_rank;
    result += TiCC::toString(rec.ngram_rank) + "#";
  }
  result += TiCC::toString(rank) + "#";
  result += TiCC::toString(rec.rank);
  return result;
}

string extractResults( const record& rec ){
  string result = rec.variant1 + "#";
  result += TiCC::toString(rec.freq1) + "#";
  result += rec.variant2 + "#";
  result += TiCC::toString(rec.freq2) + "#";
  result += TiCC::toString(rec.ld) + "#";
  result += TiCC::toString(rec.rank);
  return result;
}

template <typename TObject, typename TMember, typename TValue>
void set_val( TObject* object, TMember member, TValue value )
{
    ( *object ).*member = value;
}

template< class Tmap, typename TMember > void rank_map( const Tmap& f_map,
							vector<record*>& recs,
							TMember member ){
  if ( f_map.empty() ){
    return;
  }
  int ranking = 1;
  size_t last = f_map.begin()->first;
  for ( const auto& rit : f_map ){
    if ( rit.first < last ){
      last = rit.first;
      ++ranking;
    }
    set_val( recs[rit.second], member, ranking );
  }
}

void rank( vector<record>& records,
	   multimap<size_t, record, std::greater<size_t>>& results,
	   int clip,
	   const map<bitType,size_t>& kwc_counts,
	   const map<bitType,size_t>& kwc2_counts,
	   ostream* db, vector<bool>& skip, int factor ){
  if (verbose ){
#pragma omp critical (log)
    {
#ifdef HAVE_OPENMP
      int numt = omp_get_thread_num();
      cerr << numt << "-RANK " << records[0].variant1
	   << " " << records.size() << endl;
#else
      cerr << "RANK " << records[0].variant1
	   << " " << records.size() << endl;
#endif
    }
  }
  multimap<size_t,size_t,std::greater<size_t>> freqmap;  // freqs sorted descending
  multimap<size_t,size_t,std::greater<size_t>> f2lenmap; // f2 lenghts sorted descending
  multimap<size_t,size_t> ldmap;
  multimap<size_t,size_t, std::greater<size_t>> cslmap; // Common substring lengths deescending
  multimap<size_t,size_t,std::greater<size_t>> pairmap1;
  multimap<size_t,size_t,std::greater<size_t>> pairmap2;
  multimap<size_t,size_t,std::greater<size_t>> pairmap_combined;
  multimap<size_t,size_t,std::greater<size_t>> ngram_map;
  map<string,int> lowvarmap;
  size_t count = 0;
  vector<record*> recs;
  for ( vector<record>::iterator it = records.begin();
	it != records.end();
	++it ){
    recs.push_back( &*it );
    freqmap.insert( make_pair(it->reduced_freq2, count ) ); // freqs descending
    f2lenmap.insert( make_pair(it->f2len, count ) ); // f2lengths descending
    ldmap.insert( make_pair(it->ld,count) ); // lds sorted from low to high
    cslmap.insert( make_pair(it->csl,count) ); // csl sorted descending
    ngram_map.insert( make_pair(it->ngram_points,count) ); // ngrampoints sorted descending
    size_t var1_cnt = kwc_counts.at(it->kwc);
    it->pairs1 = var1_cnt;
    pairmap1.insert( make_pair(var1_cnt,count )); // #variants descending
    size_t var2_cnt = 0;
    try {
      var2_cnt += kwc2_counts.at(it->kwc);
    }
    catch(...){
    }
    it->pairs2 = var2_cnt;
    pairmap2.insert( make_pair(var2_cnt,count )); // #variants decending
    size_t var_combined_cnt = var1_cnt + var2_cnt;
    it->pairs_combined = var_combined_cnt;
    pairmap_combined.insert( make_pair(var_combined_cnt,count )); // #variants descending

    ++lowvarmap[it->lowervariant2]; // count frequency of variants
    ++count;
  }
  using TiCC::operator<<;
  multimap<int,size_t,std::greater<int>> lower_variantmap; // descending map
  count = 0;
  for ( const auto& it : records ){
    lower_variantmap.insert( make_pair( lowvarmap[it.lowervariant2], count ) );
    ++count;
  }
  if ( !lower_variantmap.empty() ){
    int ranking = 1;
    int last = lower_variantmap.begin()->first;
    for ( const auto& it1 : lower_variantmap ){
      if ( it1.first < last ){
	last = it1.first;
	++ranking;
      }
      recs[it1.second]->variant_count = it1.first;
      recs[it1.second]->variant_rank = ranking;
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
      recs[sit.second]->ld_rank = ranking;
    }
  }

  rank_map( freqmap, recs, &record::freq_rank );
  rank_map( f2lenmap, recs, &record::f2len_rank );
  rank_map( cslmap, recs, &record::csl_rank );
  rank_map( pairmap1, recs, &record::pairs1_rank );
  rank_map( pairmap2, recs, &record::pairs2_rank );
  rank_map( pairmap_combined, recs, &record::pairs_combined_rank );
  rank_map( ngram_map, recs, &record::ngram_rank );

  double sum = 0.0;
  vector<record*>::iterator vit = recs.begin();
  while ( vit != recs.end() ){
    double rank =
      (skip[0]?0:(*vit)->f2len_rank) +  // number of characters in the frequency
      (skip[1]?0:(*vit)->freq_rank) +   // frequency of the CC
      (skip[2]?0:(*vit)->ld_rank) +     // levenshtein distance
      (skip[3]?0:(*vit)->csl_rank) +    // common longest substring
      (skip[4]?0:(*vit)->canon_rank) +  // is it a validated word form
      (skip[5]?0:(*vit)->fl_rank) +     // first character equality
      (skip[6]?0:(*vit)->ll_rank) +     // last 2 characters equality
      (skip[7]?0:(*vit)->khc_rank) +    // known historical confusion
      (skip[8]?0:(*vit)->pairs1_rank) + //
      (skip[9]?0:(*vit)->pairs2_rank) + //
      (skip[10]?0:(*vit)->pairs_combined_rank) + //
      (skip[11]?0:(*vit)->variant_rank) + // # of decapped versions of the CC
      (skip[12]?0:(*vit)->cosine_rank) + // WordVector rank
      (skip[13]?0:(*vit)->ngram_rank);
    rank = rank/factor;
    sum += rank;
    (*vit)->rank = rank;
    ++vit;
  }

  if ( recs.size() == 1 ){
    recs[0]->rank = 1.0;
  }
  else {
    for ( auto& it : recs ){
      it->rank = 1 - it->rank/sum;
    }
  }

  multimap<size_t, record*, std::greater<size_t>> ccf_sort; // sort descending records on CC frequency;
  vit = recs.begin();
  while ( vit != recs.end() ){
    ccf_sort.insert( make_pair( (*vit)->freq2, *vit ) );
    ++vit;
  }

  int cnt = 0;
  auto it = ccf_sort.begin();
  while ( it != ccf_sort.end() ){
    ++cnt;
    if ( clip > 0 && cnt > clip )
      break;
    // store the result vector
#pragma omp critical (store)
    {
      results.insert( make_pair(it->first,*it->second) );
    }
    ++it;
  }

  if ( db ){
    vector<record*>::iterator vit = recs.begin();
    multimap<double,string,greater<double>> outv;
    while ( vit != recs.end() ){
      outv.insert( make_pair( (*vit)->rank, extractLong(**vit, skip) ) );
      ++vit;
    }
    stringstream outstr;
    for ( const auto& oit : outv ){
      outstr << oit.second << endl;
    }
#pragma omp critical (debugoutput)
    {
      *db << outstr.rdbuf();
    }
  }
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
    opts.set_long_options( "alph:,debugfile:,skipcols:,charconf:,charconfreq:,artifrq:,wordvec:,clip:,numvec:,threads:" );
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
  bool verbose = opts.extract( 'v' );
  string alfabetFile;
  string lexstatFile;
  string freqOutFile;
  string wordvecFile;
  string outFile;
  string debugFile;
  int clip = 0;
  string skipC;
  size_t artifreq = 0;
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
  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
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
	cerr << "invalid skip value in -R. All values must be between 1 and "
	     << RANK_COUNT << endl;
	exit(EXIT_FAILURE);
      }
      skip_cols.insert(kol);
    }
    if ( skip_cols.size() == (unsigned)RANK_COUNT ){
      cerr << "you may not skip all value using -R." << endl;
      exit(EXIT_FAILURE);
    }
    using TiCC::operator<<;
    cerr << "skips = " << skip_cols << endl;
    // no stuf it in a bool vector
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
    icu::UnicodeString key = TiCC::UnicodeFromUTF8(vec[0]);
    bitType value = TiCC::stringTo<bitType>( vec[2] );
    alfabet[key[0]] = value;
  }

  map<string,set<streamsize> > fileIds;
  map<bitType,size_t> kwc_counts;
  cout << "start indexing input and determining KWC counts." << endl;
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
      string variant1 = parts[0];
      fileIds[variant1].insert( pos );
      bitType kwc = TiCC::stringTo<bitType>(parts[6]);
      ++kwc_counts[kwc];
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
      icu::UnicodeString value = TiCC::UnicodeFromUTF8( vec[1] );
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
		cerr << "dus paar = 1-3: " << icu::UnicodeString(value[0]) << "-"
		     << icu::UnicodeString(value[3]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << icu::UnicodeString(value[1]) << "-"
		     << icu::UnicodeString(value[4]) << " met freq=" << counts[3] << endl;
	      }
	      kwc2_counts[key] = max + counts[3];
	    }
	    else if ( maxPos == 1 ){
	      if ( verbose ){
		cerr << "dus paar = 1-4: " << icu::UnicodeString(value[0]) << "-"
		     <<  icu::UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << icu::UnicodeString(value[1]) << "-"
		     << icu::UnicodeString(value[3]) << " met freq=" << counts[2] << endl;
	      }
	      kwc2_counts[key] = max + counts[2];
	    }
	    else if ( maxPos == 2 ){
	      if ( verbose ){
		cerr << "dus paar = 2-3: " << icu::UnicodeString(value[1]) << "-"
		     <<  icu::UnicodeString(value[3]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << icu::UnicodeString(value[0]) << "-"
		    << icu::UnicodeString(value[4]) << " met freq=" << counts[1] << endl;
	      }
	      kwc2_counts[key] = max + counts[1];
	    }
	    else if ( maxPos == 3 ){
	      if ( verbose ){
		cerr << "dus paar = 2-4: " << icu::UnicodeString(value[1]) << "-"
		     <<  icu::UnicodeString(value[4]) << " met freq=" << max << endl;
		cerr << "het anders paar: " << icu::UnicodeString(value[0]) << "-"
		    << icu::UnicodeString(value[3]) << " met freq=" << counts[0] << endl;
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

  multimap<size_t,record, std::greater<size_t>> results;
  cout << "Start the work, with " << work.size()
       << " iterations on " << numThreads << " thread(s)." << endl;
#pragma omp parallel for schedule(dynamic,1)
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
      string line;
      in.seekg( *it );
      getline( in, line );
      records.push_back( record( line, artifreq, vec ) );
      if ( verbose ){
	int tmp = 0;
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
#ifdef HAVE_OPENMP
	int numt = omp_get_thread_num();
	if ( numt == 0 && tmp % 10000 == 0 ){
#else
        if ( tmp % 10000 == 0 ){
#endif
	  cout << ".";
	  cout.flush();
	  if ( tmp % 500000 == 0 ){
	    cout << endl << tmp << endl;
	  }
	}
      }
      ++it;
    }
    ::rank( records, results, clip, kwc_counts, kwc2_counts, db, skip, skip_factor );
  }

  multimap< double, multimap<string,record*>, std::greater<double> > o_vec;
  // we sort the output of one CC frequency descending on rank, and alphabeticaly on first word
  size_t last = 0;
  for ( auto& it : results ){
    if ( last == 0 ){
      last = it.first;
    }
    else if ( last != it.first ){
      // a new key
      // output the vector;
      for ( const auto& oit : o_vec ){
        for ( const auto& vit: oit.second ){
          os << extractResults(*vit.second) << endl;
        }
      }
      o_vec.clear();
      last = it.first;
    }
    // add to the vector
    auto my_map_it = o_vec.find( it.second.rank );
    if ( my_map_it == o_vec.end () ) {
      multimap<string,record*> new_map;
      new_map.insert( make_pair( it.second.variant1, &it.second ) );
      o_vec.insert( make_pair( it.second.rank, new_map ) );
    }
    else {
      my_map_it->second.insert( make_pair( it.second.variant1, &it.second ) );
    }
  }
  if ( !o_vec.empty() ){
    // output the last items
    for ( const auto& oit : o_vec ){
      for ( const auto& vit: oit.second ){
        os << extractResults(*vit.second) << endl;
      }
    }
  }
  cout << "results in " << outFile << endl;
}
