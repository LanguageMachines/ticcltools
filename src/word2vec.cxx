#include <iostream>
#include <cassert>
#include <cmath>
#include "ticcutils/StringOps.h"
#include "ticcl/word2vec.h"

using namespace std;

bool wordvec_tester::fill( const string& name ){
  FILE *f = fopen( name.c_str(), "rb");
  if (f == NULL) {
    return false;
  }
  size_t words = 0;
  size_t size = 0;
  if ( fscanf(f, "%lud", &words) != 1 ){
    cerr << "reading #words failed" << endl;
    return false;
  }
  if ( fscanf(f, "%lud", &size) != 1 ){
    cerr << "reading #dimension failed" << endl;
    return false;
  }
  _dim = size;
  for ( unsigned b = 0; b < words; b++) {
    string word;
    while (1) {
      char kar = fgetc(f);
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
    // ansd insert in the vocabulary
    vocab.insert( make_pair( word, vec ) );
  }
  fclose(f);
  return true;
}

bool wordvec_tester::lookup( const string& sentence, size_t num_vec,
			     vector<word_dist>& result ) const {
  result.clear();
  //  cerr << "looking up: '" << word << "'" << endl;
  vector<string> words;
  size_t num_words = TiCC::split( sentence, words );
  if ( num_words < 1 ){
    cerr << "empty searchterm" << endl;
    return false;
  }

  // create an aggregated vector of all the words
  vector<float> vec( _dim, 0 );
  for ( size_t b = 0; b < num_words; ++b ) {
    unordered_map<string,vector<float>>::const_iterator it = vocab.find( words[b] );
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
  for ( const auto& it: vocab ) {
    bool hit = false;
    for ( size_t b = 0; b < num_words; ++b ){
      if ( words[b] == it.first )
	hit = true;
    }
    if ( hit ) continue;
    float dist = 0;
    for ( size_t a = 0; a < _dim; ++a ){
      dist += vec[a] * it.second[a];
    }
    for ( size_t a = 0; a < num_vec; ++a ) {
      if (dist > result[a].d ) {
	// larger so shift the rest to the back
	// and insert
	for ( size_t d = num_vec - 1; d > a; d--) {
	  result[d] = result[d-1];
	}
	result[a].w = it.first;
	result[a].d = dist;
	break;
      }
    }
  }
  return true;
}

double wordvec_tester::distance( const string& s1, const string& s2 ) const {
  //  cerr << "looking up: '" << word << "'" << endl;
  vector<string> words1;
  vector<string> words2;
  size_t num_words1 = TiCC::split( s1, words1 );
  if ( num_words1 < 1 ){
    cerr << "empty searchterm" << endl;
    return false;
  }
  size_t num_words2 = TiCC::split( s2, words2 );
  if ( num_words2 < 1 ){
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
