/*
  Copyright (c) 2006 - 2022
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
#include <cassert>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/PrettyPrint.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/Unicode.h"
#include "ticcl/ticcl_common.h"

#include "config.h"

using namespace	std;
using namespace icu;

bool verbose = false;

enum S_Class { UNDEF, UNK, PUNCT, IGNORE, CLEAN };
ostream& operator<<( ostream& os, const S_Class& cl ){
  switch ( cl ){
  case UNDEF:
    os << "Undefined";
    break;
  case CLEAN:
    os << "Clean";
    break;
  case IGNORE:
    os << "Ignore";
    break;
  case UNK:
    os << "Unknown";
    break;
  case PUNCT:
    os << "Punctuated";
    break;
  default:
    os << "WTF";
  }
  return os;
}

bool fillAlpha( istream& is, set<UChar>& alphabet ){
  int l_cnt = 0;
  int u_cnt = 0;
  int s_cnt = 0;
  string line;
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> v = TiCC::split( line );
    if ( v.size() != 3 ){
      cerr << "unsupported format for alphabet file" << endl;
      exit(EXIT_FAILURE);
    }
    UnicodeString us = TiCC::UnicodeFromUTF8( v[0] );
    us.toLower();
    if ( alphabet.find( us[0] ) == alphabet.end() ){
      alphabet.insert( us[0] );
      ++l_cnt;
    }
    us.toUpper();
    if ( alphabet.find( us[0] ) == alphabet.end() ){
      alphabet.insert( us[0] );
      ++u_cnt;
    }
    else {
      ++s_cnt;
    }
    // for now, we don't use the other fields
  }
  cout << "read an alphabet with " << u_cnt << " uppercase characters, "
       << l_cnt-s_cnt << " lowercase characters and " << s_cnt
       << " other symbols." << endl;
  return true;
}

bool fillHemp( istream& is, set<UnicodeString>& hemps ){
  string line;
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    hemps.insert( TiCC::UnicodeFromUTF8(line) );
  }
  return true;
}

bool is_ticcl_punct( UChar uc ){
  switch ( uc ){
  case '_':
    return true;
  case '^':
    return true;
  case '<':
    return true;
  case '>':
    return true;
  case '-':
    return false;
  default:
    if ( ticc_ispunct( uc ) ){
      return true;
    }
    else if ( ticc_isother( uc ) ){
      return true;
    }
    else {
      UnicodeString us(uc);
      if ( us == "\u00a0" ){
	return true;
      }
    }
  }
  return false;
}

bool depunct( const UnicodeString& us, UnicodeString& result ){
  result.remove();
  int i = 0;
  for ( ; i < us.length(); ++i ){
    // skip leading punctuation and spaces
    if ( !( is_ticcl_punct( us[i] ) || u_isspace( us[i] ) ) ){
      if ( i < us.length()-1 && us[i] == '-' && !ticc_isletter(us[i+1]) ){
	++i;
	continue;
      }
      break;
    }
  }
  int j = us.length()-1;
  for ( ; j >= 0; j-- ){
    // skip trailing punctuation and spaces
    if ( !( is_ticcl_punct( us[j] ) || u_isspace( us[j] ) ) ){
      if ( j > 0 && us[j] == '-' && !ticc_isletter(us[j-1]) ){
	--j;
	continue;
      }
      break;
    }
  }
  if ( i == 0 && j == us.length()-1 ){
    return false; // no leading/trailing puncts
  }
  else {
    for ( int k = i; k <= j; ++k ){
      result += us[k];
    }
    return true;
  }
}

static TiCC::UniFilter filter;

bool normalize_weird( const UnicodeString& in, UnicodeString& result ){
  result = filter.filter( in );
  return result != in;
}

bool is_roman( const UnicodeString& word ){
  static UnicodeString pattern = "^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$";
  static TiCC::UnicodeRegexMatcher roman_detect( pattern, "roman" );
  UnicodeString pre, post;
  bool debug = false; //(word == "IX");
  if ( debug ){
    {
      roman_detect.set_debug(debug);
      cerr << "IS Roman: test pattern = " << roman_detect.Pattern() << endl;
      cerr << "op " << word << endl;
    }
  }
  bool test;
  {
    test = roman_detect.match_all( word, pre, post );
  }
  if ( test ){
    if ( debug ){
      {
	cerr << "FOUND roman number: " << word << endl;
      }
    }
    return true;
  }
  return false;
}

S_Class classify( const UnicodeString& word,
		  const set<UChar>& alphabet ){
  int is_digit = 0;
  int is_punct = 0;
  int is_letter = 0;
  int is_space = 0;
  int is_out = 0;
  int word_len = word.length();
  if ( word_len < 2 ){
    if ( verbose ){
      cerr << "UITGANG 0: word too short: Ignore" << endl;
    }
    return IGNORE;
  }
  if ( is_roman( word ) ){
    if ( verbose ){
      cerr << "UITGANG 1: Roman number: Ignore" << endl;
    }
    return IGNORE;
  }
  for ( int i=0; i < word_len; ++i ){
    UChar uchar = word[i];
    if ( u_isspace( uchar ) ){
      ++is_space; // ignored atm
    }
    else {
      int8_t charT = u_charType( uchar );
      if ( alphabet.empty() ){
	if ( verbose ){
	  cerr << "bekijk karakter " << UnicodeString(uchar) << " van type " << toString(charT) << endl;
	}
	if ( ticc_isletter( charT ) ){
	  ++is_letter;
	}
	else if ( ticc_isdigit( charT ) ){
	  ++is_digit;
	}
	else if ( ticc_ispunct( charT ) ){
	  ++is_punct;
	}
	else if ( ticc_isother( charT ) ){
	  ++is_out;
	  // OUT
	}
	else {
	  cerr << "Warning: karakter '" << UnicodeString(uchar) << "' ("
	       << TiCC::format_nonascii( TiCC::UnicodeToUTF8(UnicodeString(uchar)) )
	       << ") is van onbekend type " << toString(charT) << endl;
	  ++is_out;
	}
      }
      else {
	if ( verbose ){
	  cerr << "bekijk karakter " << UnicodeString(uchar) << " van type " << toString(charT) << endl;
	}
	if ( alphabet.find( uchar ) != alphabet.end() ) {
	  if ( verbose ){
	    cerr << "'" << UnicodeString(uchar) << "' is IN het alfabet" << endl;
	  }
	  ++is_letter;
	}
	else if ( charT == U_DECIMAL_DIGIT_NUMBER ){
	  if ( verbose ){
	    cerr << "'" << UnicodeString(uchar) << "' is DIGIT" << endl;
	  }
	  ++is_digit;
	}
	else if ( uchar == '.' ){
	  if ( verbose ){
	    cerr << "'" << UnicodeString(uchar) << "' is PUNCT" << endl;
	  }
	  ++is_punct;
	}
	else {
	  if ( verbose ){
	    cerr << "'" << UnicodeString(uchar) << "' is OUT het alfabet" << endl;
	  }
	  ++is_out;
	}
      }
    }
  }
  if ( verbose ){
    cerr << "Classify: " << word << " IN=" << is_letter << " OUT= " << is_out << " DIG=" << is_digit << " PUNCT=" << is_punct << endl;
  }
  if ( is_digit == word_len ){
    // Filter A: gewone getallen. Worden ongemoeid gelaten, worden dus niet
    // ge-unkt of ge-ticcled, worden ook niet geteld of opgenomen in de
    // frequentielijst
    if ( verbose ){
      cerr << "UITGANG 2: Ignore" << endl;
    }
    return IGNORE;
  }
  else if ( is_digit >= is_letter + is_out ){
    // Filter B: dingen als datums, floats, of combinatie getal + een of andere
    // geldaanduiding : zelfde als getallen
    // <martin> Komt erop neer dat indien meer cijfers dan iets anders.
    if ( verbose ){
      cerr << "UITGANG 3: Ignore" << endl;
    }
    return IGNORE;
  }
  else if ( word_len >= 4 && double(is_letter + is_digit)/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 4: Clean" << endl;
    }
    return CLEAN;
  }
  else if (word_len == 3 && double(is_letter + is_digit)/word_len >= 0.66 ){
    if ( verbose ){
      cerr << "UITGANG 5: Clean" << endl;
    }
    return CLEAN;
  }
  else if (word_len < 3 && (is_out+is_digit) <1 ){
    if ( verbose ){
      cerr << "UITGANG 6: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( word_len >= 4 && double( is_letter )/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 7: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( word_len == 3 && double( is_letter )/word_len >= 0.66 ){
    if ( verbose ){
      cerr << "UITGANG 8: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( ( word_len >= 8
	      && (is_letter + is_punct) > (word_len - 2) )
	    || ( word_len >= 4
		 && (is_letter + is_punct) > (word_len - 1) ) ){
    if ( verbose ){
      cerr << "UITGANG 9: Clean" << endl;
    }
    return CLEAN;
  }
  else {
    if ( verbose ){
      cerr << "UITGANG 10: UNK" << endl;
    }
    return UNK;
  }
}

S_Class classify( const UnicodeString& us,
		  const set<UChar>& alphabet,
		  UnicodeString& punct ){
  S_Class result = CLEAN;
  punct.remove();
  UnicodeString ps;
  if ( verbose ){
    cerr << "classify:" << us << endl;
  }
  if ( depunct( us, ps  ) ){
    if ( verbose ){
      cerr << "depuncted '" << us << "' to '" << ps << "'" << endl;
    }
    if ( ps.length() == 0 ){
      // Filter C: strings met alleen maar punctuatie > UNK
      if ( us.length() < 3 )
	result = IGNORE;
      else
	result = UNK;
    }
    else {
      result = classify( ps, alphabet );
      if ( verbose ){
	cerr << "classify(" << ps << ") result:" << result << endl;
      }
      if ( result != IGNORE ){
	if ( result == CLEAN ){
	  if ( verbose ){
	    cerr << "transpose CLEAN word (" <<  us << ") to PUNCT ("
		 << ps << ")" << endl;
	  }
	  punct = ps;
	  result = PUNCT;
	}
	else {
	  result = UNK;
	}
      }
    }
  }
  else {
    if ( verbose ){
      cerr << "NOT depuncted '" << us << "'" << endl;
    }
    result = classify( us, alphabet );
  }
  if ( verbose ){
    cerr << "classify(" << us << ") ==> " << result << endl;
  }
  return result;
}

bool isAcro( const vector<UnicodeString>& parts,
	     set<UnicodeString>& result ){
  static UnicodeString pattern =  "(?:de|het|een)" + US_SEPARATOR
    + "(\\p{Lu}+)-{0,1}(?:\\p{L}*)";
  static TiCC::UnicodeRegexMatcher acro_detect( pattern, "acro_detector" );
  result.clear();
  for ( size_t i = 0; i < parts.size() -1; ++i ){
    UnicodeString us = parts[i] + US_SEPARATOR + parts[i+1];
    UnicodeString pre, post;
    // acro_detect.set_debug(1);
    // cerr << "IS ACRO: test pattern = " << acro_detect.Pattern() << endl;
    // cerr << "op " << us << endl;
    if ( acro_detect.match_all( us, pre, post ) ){
      // cerr << "IT Mached!" << endl;
      result.insert( acro_detect.get_match( 0 ) );
      using TiCC::operator<<;
      // cerr << "FOUND regexp acronym: " << result << endl;
    }
  }
  return !result.empty();
}

bool isAcro( const UnicodeString& word ){
  static UnicodeString pattern2 = "^(\\p{Lu}{1,2}\\.{1,2}(\\p{Lu}{1,2}\\.{1,2})*)(\\p{Lu}{0,2})$";
  static TiCC::UnicodeRegexMatcher acro_detect2( pattern2, "dot_alter" );
  UnicodeString pre, post;
  // acro_detect2.set_debug(1);
  // cerr << "IS ACRO: test pattern = " << acro_detect2.Pattern() << endl;
  // cerr << "op " << word << endl;
  return acro_detect2.match_all( word, pre, post );
}

UnicodeString filter_punct( const UnicodeString& us ){
  UnicodeString result;
  for ( int i=0; i < us.length(); ++i ){
    if ( !ticc_ispunct(us[i]) ){
      result += us[i];
    }
  }
  return result;
}

S_Class classify_n_gram( const vector<UnicodeString>& parts,
			 UnicodeString& end_pun,
			 unsigned int& lexclean,
			 const map<UnicodeString,unsigned int>& decap_clean_words,
			 const set<UChar>& alphabet ){
  if ( verbose ){
    cerr << "classify a " << parts.size() << "-gram" << endl;
  }
  S_Class end_cl = UNDEF;
  for ( auto const& wrd : parts ){
    if ( verbose ){
      cerr << "       classify next part: " << wrd << endl;
    }
    if ( parts.size() > 1
	 && wrd.length() == 1
	 && (!isalnum(wrd[0]) && wrd[0] != '-' ) ){
      end_cl = IGNORE;
      if ( verbose ){
	cerr << "1 letter part: " << wrd << " ==> " << end_cl << endl;
      }
    }
    if ( end_cl == IGNORE ){
      // no use to continue
      break;
    }
    S_Class cl;
    UnicodeString pun;
    UnicodeString us = wrd;
    us.toLower();
    if ( decap_clean_words.find( us ) != decap_clean_words.end() ){
      // no need to do a lot of work for already clean words
      ++lexclean;
      if ( verbose ){
	cerr << us << " in background ==> CLEAN" << endl;
      }
      cl = CLEAN;
    }
    else {
      cl = classify( wrd, alphabet, pun );
      UnicodeString us = pun;
      us.toLower();
      if ( decap_clean_words.find( us ) != decap_clean_words.end() ){
	// so the depunct word is lexically clean
	++lexclean;
	if ( verbose ){
	  cerr << pun << " in background ==> CLEAN" << endl;
	}
      }
    }
    if ( verbose ){
      cerr << "end_cl=" << end_cl << " ADD " << cl << endl;
    }
    switch( cl ){
    case IGNORE:
      if ( end_cl == CLEAN
	   && parts.size() > 2
	   && &wrd == &parts[1]
	   && wrd.length() == 1
	   && wrd[0] == '-' ){
      }
      else {
	end_cl = IGNORE;
      }
      break;
    case CLEAN:
      switch ( end_cl ){
      case UNDEF:
	end_cl = CLEAN;
	break;
      case CLEAN:
	break;
      case UNK:
	end_cl = IGNORE;
	break;
      case PUNCT:
	break;
      case IGNORE:
	break;
      }
      break;
    case PUNCT:
      switch ( end_cl ){
      case UNDEF:
	end_cl = PUNCT;
	break;
      case CLEAN:
	end_cl = PUNCT;
	break;
      case UNK:
	end_cl = IGNORE;
	break;
      case PUNCT:
	break;
      case IGNORE:
	break;
      }
      break;
    case UNK:
      switch ( end_cl ){
      case UNDEF:
	end_cl = UNK;
	break;
      case CLEAN:
	end_cl = IGNORE;
	break;
      case UNK:
	// Multiple UNK words in a sequence
	// maybe wo should add all those parts to the unk file?
	// NO: every part of an n-gram is also available as a unigram
	//     unless someone messes with the freqency lists
	end_cl = IGNORE;
	break;
      case PUNCT:
	break;
      case IGNORE:
	break;
      }
      break;
    case UNDEF:
      throw logic_error( "undef value returned by classify()" );
    }
    if ( verbose ){
      cerr << "resulting end_cl=" << end_cl << endl;
    }
    if ( pun.isEmpty() ){
      pun = wrd;
    }
    end_pun += pun + US_SEPARATOR;
  }
  if ( end_pun.length() > 0 ){
    end_pun.remove( end_pun.length()-1 );
    // remove last SEPARATOR
  }
  if ( verbose ){
    cerr << "final end_pun = " << end_pun << endl;
  }
  return end_cl;
}

void classify_one_entry( const UnicodeString& orig_word, unsigned int freq,
			 map<UnicodeString,unsigned int>& clean_words,
			 const map<UnicodeString,unsigned int>& decap_clean_words,
			 map<UnicodeString,unsigned int>& unk_words,
			 map<UnicodeString,UnicodeString>& punct_words,
			 map<UnicodeString,unsigned int>& punct_acro_words,
			 map<UnicodeString,unsigned int>& compound_acro_words,
			 bool doAcro,
			 const set<UChar>& alphabet,
			 size_t artifreq ){
  UnicodeString word;
  bool normalized = normalize_weird( orig_word, word );
  if ( verbose ){
    cerr << endl << "Run UNK on : " << orig_word;
    if ( normalized ){
      cerr << " normalized to: : " << word;
    }
    cerr << endl << endl;
  }
  vector<UnicodeString> parts = TiCC::split_at( word, US_SEPARATOR );
  if ( parts.size() == 0 ){
    return;
  }
  if ( word[0] == US_SEPARATOR[0] ){
    parts[0] = US_SEPARATOR + parts[0];
  }
  if ( word[word.length()-1] == US_SEPARATOR[0] ){
    parts.back() += US_SEPARATOR;
  }
  if ( parts.size() == 2
       && word.length() < 6 ){
    if ( verbose ){
      cerr << "too short bigram: " << word << endl;
    }
    return;
  }
  else if ( parts.size() == 3
	    && word.length() < 8 ){
    if ( verbose ){
      cerr << "too short trigram: " << word << endl;
    }
    return;
  }
  unsigned int lexclean = 0;
  UnicodeString end_pun;
  S_Class end_cl = classify_n_gram( parts, end_pun,
				    lexclean, decap_clean_words, alphabet );

  switch ( end_cl ){
  case IGNORE:
    break;
  case CLEAN:
    {
      clean_words[word] += freq;
      if ( clean_words[word] < artifreq
	   && lexclean == parts.size() ){
	clean_words[word] += artifreq;
      }
      if ( normalized ){
	punct_words[orig_word] = word;
      }
      set<UnicodeString> acros;
      if ( doAcro && isAcro( word ) ){
	if ( verbose ){
	  cerr << "CLEAN ACRO: " << word << endl;
	}
	punct_acro_words[word] += freq;
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    cerr << "CLEAN ACRO: (regex)" << word << "/" << acro << endl;
	  }
	  compound_acro_words[acro] += freq;
	}
      }
      else if ( verbose ){
	cerr << "CLEAN word: " << word << endl;
      }
    }
    break;
  case UNK:
    {
      set<UnicodeString> acros;
      if ( doAcro && isAcro( word ) ){
	if ( verbose ){
	  cerr << "UNK ACRO: " << word << endl;
	}
	clean_words[word] += freq;
	punct_acro_words[word] += freq;
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    cerr << "UNK ACRO: " << word << "/" << acro << endl;
	  }
	  compound_acro_words[acro] += freq;
	}
      }
      else {
	if ( verbose ){
	  cerr << "UNK word: " << orig_word << endl;
	}
	unk_words[orig_word] += freq;
      }
    }
    break;
  case PUNCT:
    {
      set<UnicodeString> acros;
      if ( doAcro && isAcro( end_pun ) ){
	if ( verbose ){
	  cerr << "PUNCT ACRO: " << end_pun << endl;
	}
	punct_acro_words[end_pun] += freq;
	punct_words[end_pun] = word;
	clean_words[word] += freq;
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    cerr << "PUNCT ACRO: (regex) " << word << "/" << acro << endl;
	  }
	  compound_acro_words[acro] += freq;
	}
      }
      else {
	if ( verbose ){
	  cerr << "PUNCT word: " << word << " depunct to: " << end_pun << endl;
	}
	clean_words[end_pun] += freq;
	if ( clean_words[end_pun] < artifreq
	     && lexclean == parts.size() ){
	  clean_words[end_pun] += artifreq;
	}
	punct_words[orig_word] = end_pun;
      }
    }
    break;
  case UNDEF:
    throw logic_error( "this is realy odd" );
  }
}

void format( const UnicodeString& line ){
  cerr << "\t\t\t";
  for ( int i=0; i < line.length(); ++i ){
    cerr << UnicodeString(line[i]);
    if ( line[i] == ';' ){
      cerr << endl << "\t\t\t";
    }
  }
}

UnicodeString default_filter = "æ >ae;"
  "Æ } [:Uppercase Letter:]* > AE;"
  "Æ > Ae;"
  "œ > oe;"
  "Œ } [:Uppercase Letter:]+ > OE;"
  "Œ > Oe;"
  "ĳ > ij;"
  "Ĳ > IJ;"
  "ﬂ > fl;"
  "ﬁ > fi;"
  "ﬀ > ff;"
  "ﬃ > ffi;"
  "ﬄ > ffl;"
  "ﬅ > st;"
  "ß > ss;"
  "'~' > '*';"
  "'#' > '*';"
  "[\u00a0] > '_';"
  "[[:Hyphen:][:Dash:]]+ > '-';"
  "[•·]  > '.';";

void usage( const string& name ){
  cerr << name << " [options] frequencyfile" << endl;
  cerr << "\t" << name << " will filter a wordfrequency list (in FoLiA-stats format) " << endl;
  cerr << "\t\t The output will be a cleaned wordfrequency file, an unknown "
       << endl;
  cerr << "\t\t words list and a list of words having leading/trailing "
       << endl;
  cerr << "\t\t punctuation paired with their clean variants." << endl;
  cerr << "\t-o 'name'\t create outputfile(s) with prefix 'name'" << endl;
  cerr << "\t--alph='file'\t name of the alphabet file" << endl;
  cerr << "\t--background='file'\t validated lexicon file" << endl;
  cerr << "\t--hemp='file'\t name of a Historical Emphases file" << endl;
  cerr << "\t--corpus='file'\t DEPRECATED, use --background" << endl;
  cerr << "\t--artifrq='value'\t the default value for missing frequencies"
       << endl;
  cerr << "\t\t in the validated lexicon. (default = 0)" << endl;
  cerr << "\t--acro\t also create an acronyms file. (experimental)" << endl;
  cerr << "\t--filter='file'\t use rules from 'file' to transliterate  (experimental)" << endl;
  cerr << "\t\t see http://userguide.icu-project.org/transforms/general/rules for information about rules." << endl;
  cerr << "\t\t default the following filter is used: " << endl;
  format(default_filter);
  cerr << endl;
  cerr << "\t-v\t be verbose " << endl;
  cerr << "\t-h or --help\t this message " << endl;
  cerr << "\t-V or --version\t show version " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "acro,alph:,corpus:,background:,artifrq:,filter:,help,version,hemp:" );
    opts.parse_args( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('V' ) || opts.extract("version") ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('h' ) || opts.extract("help") ){
    usage(progname);
    exit(EXIT_SUCCESS);
  }
  if ( argc < 2	){
    usage( progname );
    exit(EXIT_FAILURE);
  }
  string alphafile;
  string background_file;
  unsigned int artifreq = 0;
  bool doAcro = opts.extract("acro");
  verbose = opts.extract('v');
  opts.extract("background", background_file);
  if ( background_file.empty() ){
    opts.extract("corpus", background_file);
    if ( !background_file.empty() ){
      cerr << "WARNING!: you used the deprecated --corpus option. Please change to --background" << endl;
    }
  }
  string hemp_file;
  opts.extract("hemp", hemp_file);
  opts.extract("alph", alphafile);
  string value;
  if ( opts.extract( "artifrq", value ) ){
    if ( !TiCC::stringTo(value,artifreq) ) {
      cerr << "illegal value for --artifrq (" << value << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
  string output_name;
  opts.extract( 'o', output_name );
  string filter_file_name;
  opts.extract( "filter", filter_file_name );
  if ( !filter_file_name.empty() ){
    filter.fill( filter_file_name, "user_defined_filter" );
    UnicodeString r = filter.get_rules();
  }
  else {
    filter.init( default_filter, "default_filter" );
  }

  if ( !opts.empty() ){
    cerr << "unsupported options : " << opts.toString() << endl;
    usage(progname);
    exit(EXIT_FAILURE);
  }

  vector<string> fileNames = opts.getMassOpts();
  if ( fileNames.size() == 0 ){
    cerr << "missing frequency inputfile" << endl;
    exit(EXIT_FAILURE);
  }
  if ( fileNames.size() > 1 ){
    cerr << "only one frequency inputfile may be specified" << endl;
    exit(EXIT_FAILURE);
  }
  string file_name = fileNames[0];
  ifstream is( file_name );
  if ( !is ){
    cerr << "unable to find or open frequency file: " << file_name << endl;
    exit(EXIT_FAILURE);
  }
  if ( output_name.empty() ){
    output_name = file_name;
  }
  string unk_file_name = output_name + ".unk";
  string fore_clean_file_name = output_name + ".fore.clean";
  string all_clean_file_name = output_name + ".clean";
  string punct_file_name = output_name + ".punct";
  string acro_file_name = output_name + ".acro";

  if ( !TiCC::createPath( all_clean_file_name ) ){
    cerr << "unable to open output file: " << all_clean_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream acs( all_clean_file_name );
  if ( !background_file.empty() ){
    if ( !TiCC::createPath( fore_clean_file_name ) ){
      cerr << "unable to open output file: " << fore_clean_file_name << endl;
      exit(EXIT_FAILURE);
    }
  }
  if ( !TiCC::createPath( unk_file_name ) ){
    cerr << "unable to open output file: " << unk_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream us( unk_file_name );
  if ( !TiCC::createPath( punct_file_name ) ){
    cerr << "unable to open output file: " << punct_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream ps( punct_file_name );
  if ( doAcro ){
    if ( !TiCC::createPath( acro_file_name ) ){
      cerr << "unable to open output file: " << acro_file_name << endl;
      exit(EXIT_FAILURE);
    }
  }

  set<UChar> alphabet;

  if ( !alphafile.empty() ){
    ifstream as( alphafile );
    if ( !as ){
      cerr << "unable to open alphabet file: " << alphafile << endl;
      exit(EXIT_FAILURE);
    }
    cout << "read an alphabet from: " << alphafile << endl;
    if ( !fillAlpha( as, alphabet ) ){
      cerr << "serious problems reading alphabet file: " << alphafile << endl;
      exit(EXIT_FAILURE);
    }
  }

  set<UnicodeString> hemps;
  if ( !hemp_file.empty() ){
    ifstream hs( hemp_file );
    if ( !hs ){
      cerr << "unable to read historical emphasis file: " << hemp_file << endl;
      exit(EXIT_FAILURE);
    }
    cout << "reading Historical Emphases: " << hemp_file << endl;
    if ( !fillHemp( hs, hemps ) ){
      cerr << "serious problem reading hemps file: " << hemp_file << endl;
      exit(EXIT_FAILURE);
    }
  }
  map<UnicodeString,unsigned int> all_clean_words;
  map<UnicodeString,unsigned int> fore_clean_words;
  map<UnicodeString,unsigned int> decap_clean_words;
  map<UnicodeString,unsigned int> unk_words;
  map<UnicodeString,unsigned int> punct_acro_words;
  map<UnicodeString,unsigned int> compound_acro_words;
  map<UnicodeString,UnicodeString> punct_words;
  map<UnicodeString,unsigned int> back_lexicon;
  //  hemps.insert("F_1_o_r_e_n_t_ij_n_e_r.");
  if ( !hemps.empty() ){
    cout << "start classifying the Historical Emphases with "
	 << hemps.size() << " entries"<< endl;
    for ( const auto& hemp : hemps ){
      vector<UnicodeString> uparts = TiCC::split_at( hemp, US_SEPARATOR );
      UnicodeString clean;
      for ( const auto& u : uparts ){
	clean += u;
      }
      map<UnicodeString,UnicodeString> dummy_puncts;
      classify_one_entry( clean, 1,
			  fore_clean_words, decap_clean_words,
			  unk_words, dummy_puncts,
			  punct_acro_words, compound_acro_words,
			  doAcro, alphabet, artifreq );
      UnicodeString punct = dummy_puncts[clean];
      if ( !punct.isEmpty() ){
	punct_words[hemp] = punct;
      }
      else {
	punct_words[hemp] = clean;
      }
    }
  }

  if ( !background_file.empty() ){
    if ( artifreq == 0 ){
      cerr << "a background file is specified (--background option), but artifreq is NOT set "
	   << "(--artifrq option)" << endl;
      exit(EXIT_FAILURE);
    }
    ifstream extra( background_file );
    if ( !extra ){
      cerr << "unable to open background file: " << background_file << endl;
      exit(EXIT_FAILURE);
    }
    string line;
    while ( getline( extra, line ) ){
      vector<string> v = TiCC::split_at( line, "\t" );
      if ( v.empty() ){
	// empty line, just ignore
	continue;
      }
      if ( v.size() > 2 ){
	cerr << "background file in strange format!" << endl;
	cerr << "offending line: " << line << endl;
	exit(EXIT_FAILURE);
      }
      unsigned int freq;
      if ( v.size() == 2 ){
	if ( !TiCC::stringTo(v[1],freq) ){
	  cerr << "value of " << v[1] << " is too big to fit in an unsigned int"
	       << endl;
	  cerr << "offending line: " << line << endl;
	  exit(EXIT_FAILURE);
	}
      }
      else {
	freq = artifreq;
      }
      UnicodeString word = TiCC::UnicodeFromUTF8(v[0]);
      back_lexicon[word] = freq;
    }
    cout << "read a background lexicon with " << back_lexicon.size()
	 << " entries." << endl;

    for ( const auto& it : back_lexicon ){
      UnicodeString w = it.first;
      all_clean_words[w] += it.second;
      w.toLower();
      decap_clean_words[w] += it.second;
    }
  }
  string line;
  size_t line_cnt = 0 ;
  size_t err_cnt = 0;
  map<UnicodeString,unsigned> fore_lexicon;
  while ( getline( is, line ) ){
    ++line_cnt;
    line = TiCC::trim( line );
    if ( line.empty() ){
      continue;
    }
    vector<string> v = TiCC::split_at( line, "\t" );
    if ( v.empty() ){
      // empty line. just ignore
      continue;
    }
    if ( v.size() < 2 ){
      cerr << "error in line #" << line_cnt
	   << " content='" << line << "'" << endl;
      if ( ++err_cnt > 10 ){
	cerr << "frequency file seems to be in wrong format!" << endl;
	cerr << "too many errors, bailing out" << endl;
	exit(EXIT_FAILURE);
      }
      continue;
    }

    UnicodeString word = TiCC::UnicodeFromUTF8(v[0]);
    unsigned int freq = TiCC::stringTo<unsigned int>(v[1]);
    fore_lexicon[word] = freq;
  }
  cout << "start classifying the foreground lexicon with "
       << fore_lexicon.size() << " entries"<< endl;
  for ( const auto& wf : fore_lexicon ){
    classify_one_entry( wf.first, wf.second,
			fore_clean_words, decap_clean_words,
			unk_words, punct_words,
			punct_acro_words, compound_acro_words,
			doAcro, alphabet, artifreq );
  }
  cout << "generating output files" << endl;
  cout << "using artifrq=" << artifreq << endl;
  if ( !background_file.empty() ){
    ofstream fcs( fore_clean_file_name );
    map<unsigned int, set<UnicodeString> > wf;
    for ( const auto& it : fore_clean_words ){
      unsigned int freq = it.second;
      auto back_it = back_lexicon.find( it.first );
      if ( back_it != back_lexicon.end() ){
	// add background frequency to the foreground
	freq += back_it->second;
      }
      if ( freq > artifreq && (freq -  artifreq) > artifreq ){
      	freq -= artifreq;
      }
      wf[freq].insert( it.first );
    }
    auto wit = wf.rbegin();
    while ( wit != wf.rend() ){
      for ( const auto& sit : wit->second ){
	fcs << sit << "\t" << wit->first << endl;
      }
      ++wit;
    }
    cout << "created separate " << fore_clean_file_name << endl;
    for ( const auto& it : fore_clean_words ){
      unsigned int f1 = all_clean_words[it.first];
      unsigned int freq = it.second;
      if ( freq > artifreq && f1 >= artifreq ){
	freq -= artifreq;
      }
      all_clean_words[it.first] += freq;
    }
    wf.clear();
    for ( const auto& it : all_clean_words ){
      wf[it.second].insert( it.first );
    }
    wit = wf.rbegin();
    while ( wit != wf.rend() ){
      for ( const auto& sit : wit->second ){
	acs << sit << "\t" << wit->first << endl;
      }
      ++wit;
    }
    cout << "created " << all_clean_file_name << endl;
  }
  else {
    map<unsigned int, set<UnicodeString> > wf;
    for ( const auto& it : fore_clean_words ){
      wf[it.second].insert( it.first );
    }
    map<unsigned int, set<UnicodeString> >::const_reverse_iterator wit = wf.rbegin();
    while ( wit != wf.rend() ){
      for ( const auto& sit : wit->second ){
	acs << sit << "\t" << wit->first << endl;
      }
      ++wit;
    }
    cout << "created " << all_clean_file_name << endl;
  }
  map<unsigned int, set<UnicodeString> > wf;
  for ( const auto& uit : unk_words ){
    wf[uit.second].insert( uit.first );
  }
  auto wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for ( const auto& sit : wit->second ){
      us << sit << "\t" << wit->first << endl;
    }
    ++wit;
  }
  cout << "created " << unk_file_name << endl;

  if ( doAcro ){
    for ( const auto& ait : punct_acro_words ){
      UnicodeString ps = ait.first;
      UnicodeString us = filter_punct( ps );
      if ( compound_acro_words.find( us ) != compound_acro_words.end() ){
	// the 'dotted' word is a true acronym
	// add to the list
	compound_acro_words[ps] += ait.second;
      }
      else {
	// mishit: add to the punct file??
	//	punct_words[ait.first] += ait.first;
	//	cerr << "refuse: " << ait.first << endl;
      }
    }
    ofstream as( acro_file_name );
    for ( const auto& ait : compound_acro_words ){
      as << ait.first << "\t" << ait.second << endl;
    }
    cout << "created " << acro_file_name << endl;
  }
  for ( const auto& pit : punct_words ){
    ps << pit.first << "\t" << pit.second << endl;
  }
  cout << "created " << punct_file_name << endl;
  cout << "done!" << endl;
}
