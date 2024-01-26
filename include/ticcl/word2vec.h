/*
  Copyright (c) 2019 - 2024
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
  wordvec_tester(): _dim(0){};
  bool fill( const std::string& );
  bool lookup( const std::string&,
	       size_t,
	       std::vector<word_dist>& ) const;
  double distance( const std::string&, const std::string& ) const;
  bool analogy( const std::vector<std::string>&,
		size_t,
		std::vector<word_dist>& );
  size_t size() const { return vocab.size(); };
  size_t dimension() const { return _dim; };
 private:
  std::unordered_map<std::string,std::vector<float>> vocab;
  size_t _dim;
};

#endif
