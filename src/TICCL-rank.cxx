/*
  Copyright (c) 2006 - 2016
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
#include "ticcl/word2vec.h"
#include "ticcl/unicode.h"

using namespace std;
typedef signed long int bitType;

const int RANK_COUNT=13;

bool verbose = false;

void usage( const string& name ){
  cerr << "usage: " << name << " --alph <alphabetfile> --charconf <lexstat file>[--wordvec <wordvectorfile>] [-o <outputfile>] [-t threads] [--clip <clip>] [--debugfile <debugfile>] [--artifrq art] [--skipcols <skip>] infile" << endl;
  cerr << "\t'infile'\t is a file in TICCL-LDcalc format" << endl;
  cerr << "\t--alph 'alpha'\t an alphabet file in TICCL-lexstat format." << endl;
  cerr << "\t--charconf 'charconfus'\t a character confusion file in TICCL-lexstat format." << endl;
  cerr << "\t--charconfreq 'name'\t Extract a character confusion frequency file" << endl;
  cerr << "\t--wordvec<wordvecfile> read in a google word2vec file." << endl;
  cerr << "\t-o 'outfile'\t name of the output file." << endl;
  cerr << "\t-t 'threads'\t number of parallel threads to execute." << endl;
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
  int cls;
  double cls_rank;
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
		size_t stripM,
		const vector<word_dist>& WV ){
  vector<string> parts;
  if ( TiCC::split_at( line, parts, "~" ) != 13 ){
    ld = -1;
  }
  else {
    variant1 = parts[0];
    freq1 = TiCC::stringTo<size_t>(parts[1]);
    low_freq1 = TiCC::stringTo<size_t>(parts[2]);
    variant2 = parts[3];
    UnicodeString us = UTF8ToUnicode( variant2 );
    us.toLower();
    lowervariant2 = UnicodeToUTF8( us );
    variant_rank = -2000;
    freq2 = TiCC::stringTo<size_t>(parts[4]);
    f2len = parts[4].length();
    reduced_freq2 = freq2;
    if ( stripM > 0 && reduced_freq2 >= stripM ){
      reduced_freq2 -= stripM;
    }
    low_freq2 = TiCC::stringTo<size_t>(parts[5]);
    freq_rank = -20;
    kwc   = TiCC::stringTo<bitType>(parts[6]);
    ld    = TiCC::stringTo<int>(parts[7]);
    ld_rank = -4.5;
    cls   = TiCC::stringTo<int>(parts[8]);
    cls_rank = -5.6;
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
  result += TiCC::toString(rec.cls) + "~";
  if ( skip[3] ){
    result += "N#";
  }
  else {
    rank += rec.cls_rank;
    result += TiCC::toString(rec.cls_rank) + "#";
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


void rank( ostream& os, vector<record>& records,
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
  multimap<size_t,size_t> freqmap;
  multimap<size_t,size_t> f2lenmap;
  multimap<int,size_t> ldmap;
  multimap<int,size_t> clsmap;
  multimap<size_t,size_t> pairmap1;
  multimap<size_t,size_t> pairmap2;
  multimap<size_t,size_t> pairmap_combined;
  map<string,int> lowvarmap;
  size_t count = 0;
  vector<record*> recs;
  for ( vector<record>::iterator it = records.begin();
	it != records.end();
	++it ){
    recs.push_back( &*it );
    freqmap.insert( make_pair(it->reduced_freq2, count ) ); // freqs sorted from low to high
    f2lenmap.insert( make_pair(it->f2len, count ) ); // f2lengths sorted from low to high
    ldmap.insert( make_pair(it->ld,count) ); // lds sorted from low to high
    clsmap.insert( make_pair(it->cls,count) ); // cls sorted from low to high
    size_t var1_cnt = kwc_counts.at(it->kwc);
    it->pairs1 = var1_cnt;
    pairmap1.insert( make_pair(var1_cnt,count )); // #variants sorted from low to high
    size_t var2_cnt = 0;
    try {
      var2_cnt += kwc2_counts.at(it->kwc);
    }
    catch(...){
    }
    it->pairs2 = var2_cnt;
    pairmap2.insert( make_pair(var2_cnt,count )); // #variants sorted from low to high
    size_t var_combined_cnt = var1_cnt + var2_cnt;
    it->pairs_combined = var_combined_cnt;
    pairmap_combined.insert( make_pair(var_combined_cnt,count )); // #variants sorted from low to high

    ++lowvarmap[it->lowervariant2]; // count frequency of variants
    ++count;
  }
  using TiCC::operator<<;
  multimap<int,size_t> lower_variantmap;
  count = 0;
  for ( vector<record>::const_iterator it = records.begin();
	it != records.end();
	++it ){
    lower_variantmap.insert( make_pair( lowvarmap[it->lowervariant2], count ) );
    ++count;
  }
  multimap<int,size_t>::const_reverse_iterator it1 = lower_variantmap.rbegin();
  if ( it1 != lower_variantmap.rend() ){
    int ranking = 1;
    int last = it1->first;
    while ( it1 != lower_variantmap.rend() ){
      if ( it1->first < last ){
	last = it1->first;
	++ranking;
      }
      recs[it1->second]->variant_count = it1->first;
      recs[it1->second]->variant_rank = ranking;
      ++it1;
    }
  }

  multimap<size_t,size_t>::const_reverse_iterator rit = freqmap.rbegin();
  if ( rit != freqmap.rend() ){
    int ranking = 1;
    size_t last = rit->first;
    while ( rit != freqmap.rend() ){
      if ( rit->first < last ){
	last = rit->first;
	++ranking;
      }
      recs[rit->second]->freq_rank = ranking;
      ++rit;
    }
  }
  multimap<size_t,size_t>::const_reverse_iterator rit2 = f2lenmap.rbegin();
  if ( rit2 != f2lenmap.rend() ){
    int ranking = 1;
    size_t last = rit2->first;
    while ( rit2 != f2lenmap.rend() ){
      if ( rit2->first < last ){
	last = rit2->first;
	++ranking;
      }
      recs[rit2->second]->f2len_rank = ranking;
      ++rit2;
    }
  }

  multimap<int,size_t>::const_iterator sit = ldmap.begin();
  if ( sit != ldmap.end() ){
    int ranking = 1;
    int last = sit->first;
    while ( sit != ldmap.end() ){
      if ( sit->first > last ){
	last = sit->first;
	++ranking;
      }
      recs[sit->second]->ld_rank = ranking;
      ++sit;
    }
  }
  map<int,size_t>::const_reverse_iterator rit1 = clsmap.rbegin();
  if ( rit1 != clsmap.rend() ){
    int ranking = 1;
    int last = rit1->first;
    while ( rit1 != clsmap.rend() ){
      if ( rit1->first < last ){
	last = rit1->first;
	++ranking;
      }
      recs[rit1->second]->cls_rank = ranking;
      ++rit1;
    }
  }

  map<size_t,size_t>::const_reverse_iterator rrit = pairmap1.rbegin();
  if ( rrit != pairmap1.rend() ){
    int ranking = 1;
    size_t last = rrit->first;
    while ( rrit != pairmap1.rend() ){
      if ( rrit->first < last ){
	last = rrit->first;
	++ranking;
      }
      recs[rrit->second]->pairs1_rank = ranking;
      ++rrit;
    }
  }

  rrit = pairmap2.rbegin();
  if ( rrit != pairmap2.rend() ){
    int ranking = 1;
    size_t last = rrit->first;
    while ( rrit != pairmap2.rend() ){
      if ( rrit->first < last ){
	last = rrit->first;
	++ranking;
      }
      recs[rrit->second]->pairs2_rank = ranking;
      ++rrit;
    }
  }

  rrit = pairmap_combined.rbegin();
  if ( rrit != pairmap_combined.rend() ){
    int ranking = 1;
    size_t last = rrit->first;
    while ( rrit != pairmap_combined.rend() ){
      if ( rrit->first < last ){
	last = rrit->first;
	++ranking;
      }
      recs[rrit->second]->pairs_combined_rank = ranking;
      ++rrit;
    }
  }

  double sum = 0.0;
  vector<record*>::iterator vit = recs.begin();
  while ( vit != recs.end() ){
    double rank =
      (skip[0]?0:(*vit)->f2len_rank) +
      (skip[1]?0:(*vit)->freq_rank) +
      (skip[2]?0:(*vit)->ld_rank) +
      (skip[3]?0:(*vit)->cls_rank) +
      (skip[4]?0:(*vit)->canon_rank) +
      (skip[5]?0:(*vit)->fl_rank) +
      (skip[6]?0:(*vit)->ll_rank) +
      (skip[7]?0:(*vit)->khc_rank) +
      (skip[8]?0:(*vit)->pairs1_rank) +
      (skip[9]?0:(*vit)->pairs2_rank) +
      (skip[10]?0:(*vit)->pairs_combined_rank) +
      (skip[11]?0:(*vit)->variant_rank) +
      (skip[12]?0:(*vit)->cosine_rank);
    rank = rank/factor;
    sum += rank;
    (*vit)->rank = rank;
    ++vit;
  }
  multimap<double,string> outv;
  vit = recs.begin();
  if ( recs.size() == 1 ){
    (*vit)->rank = 1.0;
    outv.insert( make_pair( 1.0 , extractResults(**vit) ) );
  }
  else {
    while ( vit != recs.end() ){
      (*vit)->rank = 1 -  (*vit)->rank/sum;
      outv.insert( make_pair( (*vit)->rank, extractResults(**vit) ) );
      ++vit;
    }
  }
  multimap<double,string>::const_reverse_iterator oit = outv.rbegin();
  int cnt = 1;
  stringstream outstr;
  while ( oit != outv.rend() ){
    outstr << oit->second << endl;
    if ( clip != 0 && ++cnt > clip ){
      break;
    }
    ++oit;
  }
#pragma omp critical (output)
  {
    os << outstr.rdbuf();
  }
  if ( db ){
    vector<record*>::iterator vit = recs.begin();
    multimap<double,string> outv;
    while ( vit != recs.end() ){
      outv.insert( make_pair( (*vit)->rank, extractLong(**vit, skip) ) );
      ++vit;
    }
    multimap<double,string>::const_reverse_iterator oit = outv.rbegin();
    stringstream outstr;
    while ( oit != outv.rend() ){
      outstr << oit->second << endl;
      ++oit;
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
    opts.set_long_options( "alph:,debugfile:,skipcols:,charconf:,charconfreq:,artifrq:,wordvec:,clip:,numvec:" );
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
  int numThreads=1;
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
  if ( opts.extract( 't', value ) ){
    if ( !TiCC::stringTo(value,numThreads) ) {
      cerr << "illegal value for -t (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  size_t num_vec = 20;
  if ( opts.extract( "numvec", value ) ){
    if ( !TiCC::stringTo(value,num_vec) ) {
      cerr << "illegal value for --numvec (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
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
    //#define TESTWV
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
    UnicodeString key = UTF8ToUnicode(vec[0]);
    bitType value = TiCC::stringTo<bitType>( vec[2] );
    alfabet[key[0]] = value;
  }

  map<string,set<streamsize> > fileIds;
  multimap<string,record> records;
  map<bitType,size_t> kwc_counts;
  cout << "start indexing input and determining KWC counts." << endl;
  int failures = 0;
  streamsize pos = input.tellg();
  while ( getline( input, line ) ){
    if ( verbose ){
      cerr << "bekijk " << line << endl;
    }
    vector<string> parts;
    if ( TiCC::split_at( line, parts, "~" ) != 13 ){
      cerr << "invalid line: " << line << endl;
      cerr << "expected 13 ~ separated values." << endl;
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
      UnicodeString value = UTF8ToUnicode( vec[1] );
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

#ifdef HAVE_OPENMP
  omp_set_num_threads( numThreads );
#endif
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
      int tmp = 0;
      if ( verbose ){
#pragma omp critical (count)
	tmp = ++count;
	//
	// omp single isn't allowed here. trick!
#ifdef HAVE_OPENMP
	int numt = omp_get_thread_num();
#else
	int numt = 0;
#endif
	if ( numt == 0 && tmp % 10000 == 0 ){
	  cout << ".";
	  cout.flush();
	  if ( tmp % 500000 == 0 ){
	    cout << endl << tmp << endl;
	  }
	}
      }
      ++it;
    }
    ::rank( os, records, clip, kwc_counts, kwc2_counts, db, skip, skip_factor );
  }
  cout << "results in " << outFile << endl;
}
