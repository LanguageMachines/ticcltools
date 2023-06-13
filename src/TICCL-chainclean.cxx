/*
  Copyright (c) 2006 - 2023
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
#include <algorithm>
#include <functional>
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
#include "ticcutils/StringOps.h"
#include "ticcutils/CommandLine.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"
#include "ticcl/word2vec.h"

using namespace std;
using TiCC::operator<<;

set<string> follow_words;

void usage( const string& name ){
  cerr << "usage: " << name << "[options] chainfile " << endl;
  cerr << "\t\t The chainfile is an outputfile from TICCL-chain." << endl;
  cerr << "\t--lexicon A validated lexicon." << endl;
  cerr << "\t--artifrq The artifreq. Default 100000000 ." << endl;
  cerr << "\t--low=<low>\t delete records with ngrams shorter than 'low' "
       << endl;
  cerr << "\t--caseless=<yes|no> Perform caseless string comparison" << endl
       << "\t\t  (default=YES)" << endl;
  cerr << "\t\t characters. (default = 5)" << endl;
  cerr << "\t-o <outputfile> name of the outputfile." << endl;
  cerr << "\t-h or --help this message." << endl;
  cerr << "\t-v be verbose, repeat to be more verbose. " << endl;
  cerr << "\t-V or --version show version. " << endl;
  exit( EXIT_FAILURE );
}

class chain_record {
public:
  chain_record():deleted(false){};
  string variant;
  vector<string> v_parts;
  vector<string> v_dh_parts;
  string v_freq;
  string cc;
  vector<string> cc_parts;
  vector<string> cc_dh_parts;
  string cc_freq;
  string ccv;
  string ld;
  bool deleted;
};

ostream& operator<<( ostream& os, const chain_record& rec ){
  os << rec.variant << "#" << rec.v_freq << "#" << rec.cc << "#"
     << rec.cc_freq << "#";
  if ( rec.ccv != "none" ){
    os << rec.ccv << "#";
  }
  os << rec.ld << (rec.deleted?"#D":"#C");
  return os;
}

ostream& operator<<( ostream& os, const chain_record *rec ){
  os << *rec;
  return os;
}

vector<string> sort( const vector<string>& in,
		     const map<int,string>& cc_order){
  vector<string> uit;
  for ( const auto& it : cc_order ){
    auto vit = std::find( in.begin(), in.end(), it.second );
    if ( vit != in.end() ){
      uit.push_back(*vit);
    }
  }
  return uit;
}

int main( int argc, char **argv ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "lexicon:,artifrq:,follow:,low:,caseless:" );
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
  int verbosity = 0;
  while( opts.extract( 'v' ) ){
    ++verbosity;
  }
  unsigned int artifreq = 100000000;
  string value;
  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
  bool caseless = true;
  if ( opts.extract( "caseless", value ) ){
    if ( !TiCC::stringTo(value,caseless) ) {
      cerr << "illegal value for --caseless (" << value << ")" << endl;
      exit(EXIT_FAILURE);
    }
  }
  size_t low_limit = 5;
  if ( opts.extract( "low", value ) ){
    if ( !TiCC::stringTo(value,low_limit) ){
      cerr << progname << ": illegal value for --low (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  string lex_name;
  opts.extract( "lexicon", lex_name );
  if ( lex_name.empty() ){
    cerr << "missing --lexcion options"<< endl;
    exit(EXIT_FAILURE);
  }

  while ( opts.extract( "follow", value ) ){
    vector<string> parts = TiCC::split_at( value, "," );
    for ( const auto& p : parts ){
      follow_words.insert( p );
    }
  }

  string out_name;
  opts.extract( 'o', out_name );

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
  string in_name = fileNames[0];
  if ( !out_name.empty() ){
    if ( out_name == in_name ){
      cerr << "same filename for input and output!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    out_name = in_name + ".cleaned";
  }
  ifstream input( in_name );
  if ( !input ){
    cerr << "problem opening input file: " << in_name << endl;
    exit(1);
  }
  set<string> valid_words;
  ifstream lexicon( lex_name );
  string line;
  while ( getline( lexicon, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> vec = TiCC::split( line );
    if ( vec.size() < 2 ){
      cerr << progname << ": invalid line '" << line << "' in "
	   << lex_name << endl;
      exit( EXIT_FAILURE );
    }
    unsigned int freq = 0;
    if ( !TiCC::stringTo(vec[1], freq) ) {
      cerr << progname << ": invalid frequency in '" << line << "' in "
	   << lex_name << endl;
      exit( EXIT_FAILURE );
    }
    if ( freq >= artifreq ){
      if ( caseless ){
	valid_words.insert( TiCC::utf8_lowercase( vec[0] ) );
      }
      else {
	valid_words.insert( vec[0] );
      }
    }
    else {
      // the lexicon is sorted on freq. so we can bail out now
      break;
    }
  }
  cout << "read " << valid_words.size() << " validated words from "
       << lex_name << endl;
  cout << "start reading chained results" << endl;
  list<chain_record> chain_records;
  while ( getline( input, line ) ){
    vector<string> vec;
    TiCC::split_exact_at( line, vec, "#" );
    bool no_ccv = false;
    if ( vec.size() == 6 ){
      no_ccv = true;
    }
    if ( !no_ccv && vec.size() != 7 ){
      cerr << progname << ": chained file should have 6 or 7 items per line: '" << line << "' in " << lex_name << endl;
      cerr << "\t found " << vec.size() << endl;
      exit( EXIT_FAILURE );
    }
    chain_record rec;
    rec.variant = vec[0];
    rec.v_freq = vec[1];
    rec.cc = vec[2];
    rec.cc_freq = vec[3];
    if ( no_ccv ){
      rec.ccv = "none";
      rec.ld = vec[4];
    }
    else {
      rec.ccv = vec[4];
      if ( rec.ccv.empty() ){
	cerr << "YES: " << rec.variant << endl;
      }
      rec.ld = vec[5];
    }
    chain_records.push_back( rec );
  }
  cout << "start processing " << chain_records.size() << " chained results" << endl;
  map<string,int> parts_freq;
  for ( auto& rec : chain_records ){
    rec.v_parts = TiCC::split_at( rec.variant, ticcl::S_SEPARATOR );
    rec.v_dh_parts = TiCC::split_at_first_of( rec.variant, ticcl::S_SEPARATOR+"-" );
    rec.cc_parts = TiCC::split_at( rec.cc, ticcl::S_SEPARATOR );
    rec.cc_dh_parts = TiCC::split_at_first_of( rec.cc, ticcl::S_SEPARATOR+"-" );
    if ( rec.v_parts.size() == 1 ){
      continue;
    }
    for ( const auto& p : rec.v_parts ){
      string key;
      if ( caseless ){
	key = TiCC::utf8_lowercase( p );
      }
      else {
	key = p;
      }
      if ( valid_words.find( key ) == valid_words.end() ){
	++parts_freq[key];
      }
    }
  }
  cout << "found " << parts_freq.size() << " unknown parts" << endl;
  //
  // Maybe (but why) sorting parts_freq on freqency is needed?
  //
  multimap<int,string,std::greater<int>> desc_parts_freq;
  // sort on highest frequency first.
  // DOES IT REALLY MATTER???
  for ( const auto& cc : parts_freq ){
    desc_parts_freq.insert( make_pair(cc.second,cc.first) );
  }
  if ( verbosity > 0 ){
    cerr << "The unknown parts:" << endl;
    for ( const auto& it : desc_parts_freq ){
      cerr << it.first << "\t" << it.second << endl;
    }
  }

  list<chain_record*> copy_chain_records;
  for ( auto& rec : chain_records ){
    if ( rec.v_parts.size() > 1 ){
      string tmp;
      for ( const auto& p : rec.v_parts ){
	tmp += p;
      }
      if ( tmp.size() <= low_limit ){
	rec.deleted = true;
      }
    }
    copy_chain_records.push_back( &rec );
  }
  set<chain_record*> done_chain_records;
  map<string,string> done;
  size_t counter = 0;
  for ( const auto& part : desc_parts_freq ) {
    if ( ++counter % 10 == 0 ){
      cout << ".";
      cout.flush();
      if ( counter % 500 == 0 ){
	cout << endl << counter << endl;
      }
    }
    string unk_part;
    if ( caseless ){
      unk_part = TiCC::utf8_lowercase(part.second);
    }
    else {
      unk_part = part.second;
    }
    bool show = (verbosity>0)
      || follow_words.find( unk_part ) != follow_words.end();
    if ( show ){
      cerr << "\n  Loop for part: " << part.second << "/" << unk_part << endl;
    }
    map<string,int> cc_freqs;
    map<int,string> cc_order;
    int oc = 0;
    for ( const auto& it : chain_records ){
      bool match = false;
      for ( const auto& p : it.v_dh_parts ){
	string v_part;
	if ( caseless ){
	  v_part = TiCC::utf8_lowercase(p);
	}
	else {
	  v_part = p;
	}
	if ( verbosity>1 ){
	  cerr << "ZOEK: " << v_part << endl;
	}
	if ( v_part == unk_part ){
	  if ( show ){
	    cerr << "found: " << unk_part << " in: " << it << endl;
	  }
	  match = true;
	  break;
	}
      }
      if ( match ){
	for ( const auto& cp : it.cc_dh_parts ){
	  string c_part;
	  if ( caseless ){
	    c_part = TiCC::utf8_lowercase(cp);
	  }
	  else {
	    c_part = cp;
	  }
	  if ( cc_freqs.find(c_part) == cc_freqs.end() ){
	    // first encounter
	    cc_order[oc++] = c_part;
	  }
	  ++cc_freqs[c_part];
	  if ( show ){
	    cerr << "for: " << unk_part << " increment " << c_part << endl;
	  }
	}
      }
    }
    multimap<int,string,std::greater<int>> desc_cc;
    set<int> keys;
    // sort on highest frequency first.
    // DOES IT REALLY MATTER???
    for ( const auto& cc : cc_freqs ){
      keys.insert(cc.second);
      desc_cc.insert( make_pair(cc.second,cc.first) );
    }
    if ( show ){
      cerr << "found " << desc_cc.size() << " CC's for: " << unk_part << endl;
      for ( const auto& it : desc_cc ){
	cerr << it.first << "\t" << it.second << endl;
      }
    }
    map<int,vector<string>,std::greater<int>> desc_cc_vec_map;
    for ( const auto& key : keys ){
      auto const& pr = desc_cc.equal_range( key );
      vector<string> in;
      for ( auto it = pr.first; it != pr.second; ++it ){
	in.push_back( it->second );
      }
      vector<string> uit = sort(in,cc_order);
      desc_cc_vec_map[key] = uit;
    }
    if ( show ){
      cerr << "found " << cc_order.size() << " CC's for: " << unk_part << endl;
      for ( const auto& it : desc_cc_vec_map ){
	cerr << it.first << "\t" << it.second << endl;
      }
    }
    for ( const auto& dvm_it : desc_cc_vec_map ){
      if ( show ){
	cerr << "With frequency = " << dvm_it.first << endl;
      }
      for ( const auto& dcc : dvm_it.second ){
	string cand_cor;
	if ( caseless ){
	  cand_cor = TiCC::utf8_lowercase( dcc );
	}
	else {
	  cand_cor = dcc;
	}
	if ( show ){
	  cerr << "BEKIJK: " << cand_cor << "[" << dvm_it.first << "]" << endl;
	}
	map<string,int> uniq;
	auto it = copy_chain_records.begin();
	while ( it != copy_chain_records.end() ){
	  chain_record* rec = *it;
	  if ( rec->deleted ){
	    ++it;
	    continue;
	  }
	  if ( done_chain_records.find( rec ) != done_chain_records.end() ){
	    if ( show && rec->variant.find( unk_part) != string::npos ) {
	      cerr << "skip already done " << rec << endl;
	    }
	    ++it;
	    continue;
	  }
	  if ( rec->v_parts.size() == 1 ){
	    string vari;
	    string corr;
	    if ( caseless ){
	      vari = TiCC::utf8_lowercase( rec->variant );
	      corr = TiCC::utf8_lowercase( rec->cc );
	    }
	    else {
	      vari = rec->variant;
	      corr = rec->cc;
	    }
	    if ( vari == unk_part
		 && corr.find(cand_cor) != string::npos ){
	      // this is (might be) THE desired CC
	      if ( show ){
		cerr << "UNI gram: both " << unk_part << " and " << cand_cor
		     << " matched in: " << rec << endl;
		cerr << "KEEP: " << rec << endl;
	      }
	      done[corr] = vari;
	      done_chain_records.insert(rec);
	      if ( rec->cc_parts.size() == 1 ){
		// so this is a unigram CC
		++uniq[vari];
	      }
	    }
	  }
	  else {
	    bool local_show = verbosity > 0;
	    for ( const auto& p : rec->v_parts ){
	      local_show |= follow_words.find( p ) != follow_words.end();
	    }
	    if ( local_show ){
	      cerr << "bekijk met " << cand_cor << ":" << rec << endl;
	    }
	    for ( const auto& vp : rec->v_parts ){
	      if ( uniq.find(vp) != uniq.end() ){
		// a ngram part equals an already resolved unigram
		// discard!
		rec->deleted = true;
		break;
	      }
	    }
	    if ( rec->deleted ){
	      if ( local_show ){
		cerr << "REMOVE uni: " << rec << endl;
	      }
	      ++it;
	      continue;
	    }
	    bool match = false;
	    for( const auto& cp : rec->cc_parts ){
	      string cor_part;
	      if ( caseless ){
		cor_part = TiCC::utf8_lowercase( cp );
	      }
	      else {
		cor_part = cp;
	      }
	      if ( cand_cor == cor_part ){
		// CC match
		for ( const auto& p : rec->v_parts ){
		  string p_part;
		  if ( caseless ){
		    p_part = TiCC::utf8_lowercase( p );
		  }
		  else {
		    p_part = p;
		  }
		  if ( p_part == unk_part ){
		    // variant match too
		    match = true;
		    break;
		  }
		}
		if ( match ){
		  if ( local_show ){
		    cerr << "both " << cor_part << " and " << unk_part
			 << " matched in: " << rec << endl;
		  }
		  string lvar;
		  if ( caseless ){
		    lvar = TiCC::utf8_lowercase(rec->variant);
		  }
		  else {
		    lvar = rec->variant;
		  }
		  if ( done.find( cor_part ) != done.end() ){
		    string v = done[cor_part];
		    if ( uniq.find( unk_part ) != uniq.end() ){
		      if ( local_show ){
			cerr << "REMOVE uni: " << rec << endl;
		      }
		      (*it)->deleted = true;
		    }
		    else if ( lvar.find( v ) != string::npos ){
		      if ( local_show ){
			cerr << "REMOVE match: " << rec << endl;
		      }
		      (*it)->deleted = true;
		    }
		    else {
		      if ( local_show ){
			cerr << "KEEP 1: " << rec << endl;
		      }
		      done[cor_part] = lvar;
		      done_chain_records.insert(rec);
		    }
		  }
		  else {
		    if ( local_show ){
		      cerr << "KEEP 2: " << rec << endl;
		    }
		    done[cor_part] = lvar;
		    done_chain_records.insert(rec);
		  }
		  break;
		}
	      }
	    }
	  }
	  ++it;
	}
      }
    }
  }
  ofstream os( out_name );
  int count = 0;
  for ( const auto it : copy_chain_records ){
    if ( it != 0 ){
      if ( !it->deleted ){
	++count;
	os << it << endl;
      }
    }
  }
  cerr << endl << "wrote " << count << " chain_records to " << out_name << endl;
  ofstream osd( out_name + ".deleted" );
  count = 0;
  for ( const auto it : copy_chain_records ){
    if ( it != 0 ){
      if ( it->deleted ){
	++count;
	osd << it << endl;
      }
    }
  }
  cerr << "wrote " << count << " DELETED records to " << out_name
       << ".deleted" << endl;
}
