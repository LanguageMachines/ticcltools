/*
  Copyright (c) 2019 - 2026
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
#include <cassert>
#include <cmath>
#include "ticcutils/StringOps.h"
#include "ticcl/word2vec.h"

using namespace std;

bool wordvec_tester::fill( const string& name ){
  FILE *f = fopen( name.c_str(), "rb");
  if (f == NULL) {
    cerr << "unable to open " << name << endl;
    return false;
  }
  unsigned long words = 0;
  if ( fscanf(f, "%lud", &words) != 1 ){
    cerr << "reading #words failed" << endl;
    fclose(f);
    return false;
  }
  unsigned long dim = 0;
  if ( fscanf(f, "%lud", &dim) != 1 ){
    cerr << "reading #dimension failed" << endl;
    fclose(f);
    return false;
  }
  _dim = dim;
  cout << "start reading " << words << " vectors, dim=" << dim << endl;
  for ( unsigned b = 0; b < words; b++) {
    string word;
    while (1) {
      int kar = fgetc(f);
      if ( feof(f) || kar == EOF || kar == ' ' )
	break;
      if ( kar != '\n' ){
	word += kar;
      }
    }
    vector<float> vec( _dim );
    for ( size_t a = 0; a < _dim; a++) {
      float val;
      if ( fread(&val, sizeof(float), 1, f ) != 1 ){
	cerr << "reading float failed" << endl;
	exit(1);
      }
      vec[a] = val;
    }
    // normalize the vector
    float len = 0;
    for ( size_t i = 0; i < _dim; ++i )
      len += vec[i] * vec[i];
    len = sqrt(len);
    for ( size_t i = 0; i < _dim; ++i ){
      vec[i] /= len;
    }
    // and insert in the vocabulary
    vocab.insert( make_pair( word, vec ) );
  }
  fclose(f);
  return true;
}

bool wordvec_tester::lookup( const string& sentence, size_t num_vec,
			     vector<word_dist>& result ) const {
  result.clear();
  //  cerr << "looking up: '" << sentence << "'" << endl;
  vector<string> words = TiCC::split( sentence );
  if ( words.empty() ){
    cerr << "empty searchterm" << endl;
    return false;
  }

  // create an aggregated vector of all the words
  vector<float> vec( _dim, 0 );
  for ( size_t b = 0; b < words.size(); ++b ) {
    auto const it = vocab.find( words[b] );
    if ( it == vocab.end() ){
      //      cerr << "couldn't find " << words[a] << endl;
      return false;
    }
    for ( size_t a = 0; a < _dim; ++a ){
      vec[a] += it->second[a];
    }
  }

  // normalize the created vector
  float len = 0;
  for ( size_t a = 0; a < _dim; ++a ) {
    len += vec[a] * vec[a];
  }
  len = sqrt(len);
  for ( size_t a = 0; a < _dim; ++a ) {
    vec[a] /= len;
  }

  // now compare with ALL the vectors in de vocabulary
  // keep de 'num_vec' largest

  result.resize( num_vec, {"", 0.0 } );
  for ( const auto& [word,p_vec]: vocab ) {
    bool hit = false;
    for ( const auto& w_it : words ){
      if ( w_it == word ){
	hit = true;
	break;
      }
    }
    if ( hit ) continue;
    float dist = 0;
    for ( size_t a = 0; a < _dim; ++a ){
      dist += vec[a] * p_vec[a];
    }
    for ( size_t a = 0; a < num_vec; ++a ) {
      if (dist > result[a].d ) {
	// larger so shift the rest to the back
	// and insert
	for ( size_t d = num_vec - 1; d > a; d--) {
	  result[d] = result[d-1];
	}
	result[a].w = word;
	result[a].d = dist;
	break;
      }
    }
  }
  return true;
}

bool wordvec_tester::analogy( const vector<string>& words,
			      size_t num_vec,
			      vector<word_dist>& result ){
  result.clear();
  if ( words.size() != 3 ){
    cerr << "normalize needs 3 words, not " << words.size() << endl;
    return false;
  }
    // create an aggregated vector of all the words
  vector<float> vec( _dim, 0 );

  auto const& it0 = vocab.find( words[0] );
  if ( it0 == vocab.end() ){
    //      cerr << "couldn't find " << words[0] << endl;
    return false;
  }
  auto const& it1 = vocab.find( words[1] );
  if ( it1 == vocab.end() ){
    //      cerr << "couldn't find " << words[1] << endl;
    return false;
  }
  auto const& it2 = vocab.find( words[2] );
  if ( it2 == vocab.end() ){
    //      cerr << "couldn't find " << words[2] << endl;
    return false;
  }

  for ( size_t a = 0; a < _dim; ++a ){
    vec[a] += it1->second[a] - it0->second[a] + it2->second[a];
  }

  // normalize the created vector
  float len = 0;
  for ( size_t a = 0; a < _dim; ++a ) {
    len += vec[a] * vec[a];
  }
  len = sqrt(len);
  for ( size_t a = 0; a < _dim; ++a ) {
    vec[a] /= len;
  }

  // now compare with ALL the vectors in de vocabulary
  // keep de 'num_vec' largest
  result.resize( num_vec, {"", 0.0 } );
  for ( const auto& [word,p_vec]: vocab ) {
    if ( word == it0->first
	 || word == it1->first
	 || word == it2->first ){
      //      cerr << "skip " << word << endl;
      continue;
    }
    float dist = 0;
    for ( size_t a = 0; a < _dim; ++a ){
      dist += vec[a] * p_vec[a];
    }
    for ( size_t a = 0; a < num_vec; ++a ) {
      if (dist > result[a].d ) {
	// larger so shift the rest to the back
	// and insert
	for ( size_t d = num_vec - 1; d > a; d--) {
	  result[d] = result[d-1];
	}
	result[a].w = word;
	result[a].d = dist;
	break;
      }
    }
  }
  return true;
}

double wordvec_tester::distance( const string& s1, const string& s2 ) const {
  //  cerr << "looking up: '" << word << "'" << endl;
  vector<string> words1 = TiCC::split( s1 );
  if ( words1.empty() ){
    cerr << "empty searchterm" << endl;
    return false;
  }
  vector<string> words2 = TiCC::split( s2 );
  if ( words2.empty() ){
    cerr << "empty searchterm" << endl;
    return false;
  }

  // create an aggregated vector of all the words
  vector<float> vec1( _dim, 0 );
  for ( auto const& w : words1 ) {
    auto const& it = vocab.find( w );
    if ( it == vocab.end() ){
      throw "unknown word '" + w + "'";
    }
    for ( size_t a = 0; a < _dim; ++a ){
      vec1[a] += it->second[a];
    }
  }
  vector<float> vec2( _dim, 0 );
  for ( auto const& w : words2 ) {
    auto const& it = vocab.find( w );
    if ( it == vocab.end() ){
      throw "unknown word '" + w + "'";
    }
    for ( size_t a = 0; a < _dim; ++a ){
      vec2[a] += it->second[a];
    }
  }

  // normalize the created vectors
  float len1 = 0;
  float len2 = 0;
  for ( size_t a = 0; a < _dim; ++a ) {
    len1 += vec1[a] * vec1[a];
    len2 += vec2[a] * vec2[a];
  }
  len1 = sqrt(len1);
  len2 = sqrt(len2);
  for ( size_t a = 0; a < _dim; ++a ) {
    vec1[a] /= len1;
    vec2[a] /= len2;
  }

  // now inproduct the two vectors
  double result = 0.0;
  for ( size_t a = 0; a < _dim; ++a ){
    result += vec1[a] * vec2[a];
  }
  return result;
}
