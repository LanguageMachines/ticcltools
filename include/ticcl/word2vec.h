#ifndef WORD2VEC_H
#define WORD2VEC_H

#include <unordered_map>
#include <vector>
#include <string>

struct word_dist {
  std::string w;
  float d;
};

class wordvec_tester {
 public:
  bool fill( const std::string& );
  bool lookup( const std::string&,
	       size_t,
	       std::vector<word_dist>& ) const;
  double distance( const std::string&, const std::string& ) const;
  size_t size() const { return vocab.size(); };
  size_t dimension() const { return _dim; };
 private:
  std::unordered_map<std::string,std::vector<float>> vocab;
  size_t _dim;
};

#endif
