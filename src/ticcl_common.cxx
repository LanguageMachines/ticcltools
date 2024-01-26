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

#include "ticcutils/StringOps.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"

#include <cstdlib>
#include <string>
#include <vector>

using namespace icu;
using namespace std;

namespace ticcl {
  const bitType HonderdHash = high_five( 100 );
  const bitType HonderdEenHash = high_five( 101 );

  bitType high_five( int val ){
    bitType result = val;
    result *= val;
    result *= val;
    result *= val;
    result *= val;
    return result;
  }

  bitType hash( const UnicodeString& s,
		const map<UChar,bitType>& alphabet,
		bool debug ){
    UnicodeString us = s;
    us.toLower();
    bitType result = 0;
    bool multPunct = false;
    for( int i=0; i < us.length(); ++i ){
      auto it = alphabet.find( us[i] );
      if ( it != alphabet.end() ){
	result += it->second;
	if ( debug ){
	  cerr << "  CHAR, add " << UnicodeString( us[i] ) << " "
	       << it->second << " ==> " << result << endl;
	}
      }
      else {
	int8_t charT = u_charType( us[i] );
	if ( u_isspace( us[i] ) ){
	  continue;
	}
	else if ( ispunct( charT ) ){
	  if ( !multPunct ){
	    result += HonderdHash;
	    if ( debug ){
	      cerr << "PUNCT, add " << UnicodeString( us[i] ) << " "
		   << HonderdHash	 << " ==> " << result << endl;
	    }
	    multPunct = true;
	  }
	}
	else {
	  result += HonderdEenHash;
	  if ( debug ){
	    cerr << "   UNK, add " << UnicodeString( us[i] ) << " "
		 << HonderdHash << " ==> " << result << endl;
	  }
	}
      }
    }
    return result;
  }

  unsigned int ldCompare( const UnicodeString& s1, const UnicodeString& s2 ){
    const size_t len1 = s1.length(), len2 = s2.length();
    vector<unsigned int> col(len2+1), prevCol(len2+1);
    for ( unsigned int i = 0; i < prevCol.size(); ++i ){
      prevCol[i] = i;
    }
    for ( unsigned int i = 0; i < len1; ++i ) {
      col[0] = i+1;
      for ( unsigned int j = 0; j < len2; ++j )
	col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]),
			prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
      col.swap(prevCol);
    }
    unsigned int result = prevCol[len2];
    return result;
  }

  bool fillAlphabet( istream& is,
		     map<UChar,bitType>& alphabet,
		     int clip ){
    UnicodeString line;
    while ( TiCC::getline( is, line ) ){
      if ( line.length() == 0 || line[0] == '#' ){
	continue;
      }
      vector<UnicodeString> v = TiCC::split_at( line, "\t" );
      if ( v.size() != 3 ){
	throw runtime_error( "unsupported format for alphabet file" );
      }
      int freq = TiCC::stringTo<int>( v[1] );
      if ( freq > clip || freq == 0 ){
	// freq = 0 is special, for separator
	UnicodeString v0 = v[0];
	bitType hash_val = TiCC::stringTo<bitType>( v[2] );
	alphabet[v0[0]] = hash_val;
      }
    }
    return true;
  }

  set<bitType> read_bit_set( istream& is ){
    set<bitType> result;
    while ( is ){
      bitType bit;
      is >> bit;
      is.ignore( INT_MAX, '\n' );
      result.insert( bit );
    }
    return result;
  }

  set<bitType> read_anahash( istream& is,
			     const int& low,
			     const int& high,
			     size_t& skipped,
			     bool verbose ){
    set<bitType> result;
    UnicodeString line;
    while ( TiCC::getline( is, line ) ){
      vector<UnicodeString> parts = TiCC::split_at( line, "~" );
      if ( parts.size() > 1 ){
	bitType bit = TiCC::stringTo<bitType>( parts[0] );
	vector<UnicodeString> parts2 = TiCC::split_at( parts[1], "#" );
	if ( parts2.size() > 0 ){
	  UnicodeString firstItem = parts2[0];
	  if ( firstItem.length() >= low &&
	       firstItem.length() <= high ){
	    result.insert( bit );
	  }
	  else {
	    if ( verbose ){
	      cerr << "skip " << parts2[0] << endl;
	    }
	    ++skipped;
	  }
	}
      }
    }
    return result;
  }

  set<bitType> read_confusions( istream& is ){
    set<bitType> result;
    size_t count = 0;
    UnicodeString line;
    while ( TiCC::getline( is, line ) ){
      if ( ++count % 1000 == 0 ){
	cout << ".";
	cout.flush();
	if ( count % 50000 == 0 ){
	  cout << endl << count << endl;;
	}
      }
      vector<UnicodeString> parts = TiCC::split_at( line, "#" );
      if ( parts.size() > 0 ){
	bitType bit = TiCC::stringTo<bitType>( parts[0] );
	result.insert( bit );
      }
      else {
	cerr << "problems with line " << line << endl;
	cerr << "bail out " << endl;
	exit(1);
      }
    }
    return result;
  }

} // namespace ticcl
