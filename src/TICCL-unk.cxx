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
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcutils/FileUtils.h"
#include "ticcutils/Unicode.h"
#include "ticcl/unicode.h"

#include "config.h"
#undef HAVE_OPENMP
#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace	std;

const string SEPARATOR = "_";

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
  string line;
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> v;
    int n = TiCC::split( line, v );
    if ( n != 3 ){
      cerr << "unsupported format for alphabet file" << endl;
      exit(EXIT_FAILURE);
    }
    UnicodeString us = TiCC::UnicodeFromUTF8( v[0] );
    us.toLower();
    alphabet.insert( us[0] );
    us.toUpper();
    alphabet.insert( us[0] );
    // for now, we don't use the other fields
  }
  return true;
}

bool is_ticcl_punct( UChar uc ){
  switch (  uc ){
  case '^':
    return true;
  case '-':
    return false;
  default:
    return ticc_ispunct( uc );
  }
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
    if ( verbose ){
      cerr << "depunct '" << us << "' ==> '" << result << "'" << endl;
    }
    return true;
  }
}

static TiCC::UniFilter filter;

bool normalize_weird( const UnicodeString& in, UnicodeString& result ){
  UnicodeString us = filter.filter( in );
  result = us;
  return result != in;
}

bool is_roman( const UnicodeString& word ){
  static string pattern = "^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$";
  static TiCC::UnicodeRegexMatcher roman_detect( TiCC::UnicodeFromUTF8(pattern), "roman" );
  UnicodeString pre, post;
  bool debug = false; //(word == "IX");
  if ( debug ){
    //#pragma omp critical (debug)
    {
      roman_detect.set_debug(debug);
      cerr << "IS Roman: test pattern = " << roman_detect.Pattern() << endl;
      cerr << "op " << word << endl;
    }
  }
  bool test;
  //#pragma omp critical (roman)
  {
    test = roman_detect.match_all( word, pre, post );
  }
  if ( test ){
    if ( debug ){
      //#pragma omp critical (debug)
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
      cerr << "UITGANG 0: word to short: Ignore" << endl;
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
      ++is_space;
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

S_Class classify( const string& word, const set<UChar>& alphabet,
		  string& punct ){
  S_Class result = CLEAN;
  punct.clear();
  UnicodeString us = TiCC::UnicodeFromUTF8( word );
  UnicodeString ps;
  if ( depunct( us, ps  ) ){
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
	cerr << "classify: " << result << endl;
      }
      if ( result != IGNORE ){
	if ( result == CLEAN ){
	  if ( verbose ){
	    cerr << "transpose CLEAN word (" <<  us << ") to PUNCT ("
		 << ps << ")" << endl;
	  }
	  punct = TiCC::UnicodeToUTF8( ps );
	  result = PUNCT;
	}
	else {
	  result = UNK;
	}
      }
    }
  }
  else {
    result = classify( us, alphabet );
  }
  if ( verbose ){
    cerr << "classify(" << word << ") ==> " << result << endl;
  }
  return result;
}

bool isAcro( const vector<string>& parts,
	     set<string>& result ){
  static string pattern =  "(?:de|het|een)" + SEPARATOR + "(\\p{Lu}+)-(?:\\p{L}*)";
  static TiCC::UnicodeRegexMatcher acro_detect( TiCC::UnicodeFromUTF8(pattern), "acro_detector" );
  result.clear();
  for ( size_t i = 0; i < parts.size() -1; ++i ){
    UnicodeString us = TiCC::UnicodeFromUTF8( parts[i] + SEPARATOR + parts[i+1] );
    UnicodeString pre, post;
    //#pragma omp critical (debug)
    //    {
    //       acro_detect.set_debug(1);
    //    }
    //  cerr << "IS ACRO: test pattern = " << acro_detect.Pattern() << endl;
    //  cerr << "op " << us << endl;
    //#pragma omp critical (acro)
    {
      if ( acro_detect.match_all( us, pre, post ) ){
	//    cerr << "IT Mached!" << endl;
	result.insert( TiCC::UnicodeToUTF8(acro_detect.get_match( 0 )) );
	//    cerr << "FOUND regexp acronym: " << result << endl;
      }
    }
  }
  return !result.empty();
}

bool isAcro( const string& word ){
  static string pattern2 = "^(\\p{Lu}{1,2}\\.{1,2}(\\p{Lu}{1,2}\\.{1,2})*)(\\p{Lu}{0,2})$";
  static TiCC::UnicodeRegexMatcher acro_detect2( TiCC::UnicodeFromUTF8(pattern2), "dot_alter" );
  UnicodeString pre, post;
  //  acro_detect2.set_debug(1);
  UnicodeString us = TiCC::UnicodeFromUTF8( word );
  //  cerr << "IS ACRO: test pattern = " << acro_detect2.Pattern() << endl;
  //  cerr << "op " << us << endl;
  bool test;
  //#pragma omp critical (acro2)
  {
    test = acro_detect2.match_all( us, pre, post );
  }
  return test;
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

S_Class classify_n_gram( const vector<string>& parts,
			 string& end_pun,
			 unsigned int& lexclean,
			 const map<UnicodeString,unsigned int>& decap_clean_words,
			 const set<UChar>& alphabet ){
  S_Class end_cl = UNDEF;
  for ( auto const& wrd : parts ){
    if ( parts.size() > 1
	 && wrd.size() == 1
	 && (!isalnum(wrd[0]) && wrd[0] != '-' ) ){
      end_cl = IGNORE;
      if ( verbose ){
	cerr << "1 letter part: " << wrd << " ==> " << end_cl << endl;
      }
    }
    if ( end_cl == IGNORE ){
      break;
    }
    S_Class cl;
    string pun;
    UnicodeString us = TiCC::UnicodeFromUTF8( wrd );
    us.toLower();
    if ( decap_clean_words.find( us ) != decap_clean_words.end() ){
      // no need to do a lot of work for already clean words
      ++lexclean;
      cl = CLEAN;
    }
    else {
      cl = classify( wrd, alphabet, pun );
    }
    if ( verbose ){
      cerr << "end_cl=" << end_cl << " ADD " << cl << endl;
    }
    switch( cl ){
    case IGNORE:
      if ( end_cl == CLEAN
	   && parts.size() > 2
	   && &wrd == &parts[1]
	   && wrd.size() == 1
	   && wrd[0] == '-' ){
      }
      else {
	end_cl = IGNORE;
      }
      break;
    case CLEAN:
      if ( end_cl == UNDEF ){
	end_cl = CLEAN;
      }
      else if ( end_cl == CLEAN ){
      }
      else if ( end_cl == UNK ){
	end_cl = IGNORE;
      }
      else if ( end_cl == PUNCT ){
      }
      break;
    case PUNCT:
      if ( end_cl == UNDEF ){
	end_cl = PUNCT;
      }
      else if ( end_cl == CLEAN ){
	end_cl = PUNCT;
      }
      else if ( end_cl == UNK ){
	end_cl = IGNORE;
      }
      else if ( end_cl == PUNCT ){
      }
      break;
    case UNK:
      if ( end_cl == UNDEF ){
	end_cl = UNK;
      }
      else if ( end_cl == UNK ){
	end_cl = IGNORE;
      }
      else if ( end_cl == CLEAN ){
	end_cl = IGNORE;
      }
      break;
    case UNDEF:
      throw logic_error( "undef value returned by classify()" );
    }
    if ( verbose ){
      cerr << "resulting end_cl=" << end_cl << endl;
    }
    if ( pun.empty() ){
      pun = wrd;
    }
    end_pun += pun + SEPARATOR;
  }
  end_pun = TiCC::trim_back( end_pun, SEPARATOR );
  return end_cl;
}

void classify_one_entry( const string& orig_word, unsigned int freq,
			 map<string,unsigned int>& clean_words,
			 const map<UnicodeString,unsigned int>& decap_clean_words,
			 map<string,unsigned int>& unk_words,
			 map<string,string>& punct_words,
			 map<string,unsigned int>& punct_acro_words,
			 map<string,unsigned int>& compound_acro_words,
			 bool doAcro,
			 const set<UChar>& alphabet,
			 size_t artifreq ){
  UnicodeString us = TiCC::UnicodeFromUTF8( orig_word );
  UnicodeString nus;
  bool normalized = normalize_weird( us, nus );
  string word = TiCC::UnicodeToUTF8( nus );
  if ( verbose ){
    //#pragma omp critical (verbose)
    {
      cerr << endl << "Run UNK on : " << orig_word;
      if ( normalized ){
	cerr << " normalized to: : " << word;
      }
      cerr << endl << endl;
    }
  }
  vector<string> parts;
  TiCC::split_at( word, parts, SEPARATOR );
  if ( parts.size() == 0 ){
    return;
  }
  else if ( parts.size() == 2
	    && word.size() < 6 ){
    if ( verbose ){
      //#pragma omp critical (verbose)
      {
	cerr << "to short bigram: " << word << endl;
      }
    }
    return;
  }
  else if ( parts.size() == 3
	    && word.size() < 8 ){
    if ( verbose ){
      //#pragma omp critical (verbose)
      {
	cerr << "to short trigram: " << word << endl;
      }
    }
    return;
  }
  unsigned int lexclean = 0;
  string end_pun;
  S_Class end_cl = classify_n_gram( parts, end_pun,
				    lexclean, decap_clean_words, alphabet );

  switch ( end_cl ){
  case IGNORE:
    break;
  case CLEAN:
    {
      //#pragma omp critical (classify)
      {
	clean_words[word] += freq;
	if ( clean_words[word] < artifreq
	     && lexclean == parts.size() ){
	  clean_words[word] += artifreq;
	}
	if ( normalized ){
	  punct_words[orig_word] = word;
	}
      }
      set<string> acros;
      if ( doAcro && isAcro( word ) ){
	if ( verbose ){
	  //#pragma omp critical (verbose)
	  {
	    cerr << "CLEAN ACRO: " << word << endl;
	  }
	}
	//#pragma omp critical (classify)
	{
	  punct_acro_words[word] += freq;
	}
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    //#pragma omp critical (verbose)
	    {
	      cerr << "CLEAN ACRO: (regex)" << word << "/" << acro << endl;
	    }
	  }
	  //#pragma omp critical (classify)
	  {
	    compound_acro_words[acro] += freq;
	  }
	}
      }
      else if ( verbose ){
	//#pragma omp critical (verbose)
	{
	  cerr << "CLEAN word: " << word << endl;
	}
      }
    }
    break;
  case UNK:
    {
      set<string> acros;
      if ( doAcro && isAcro( word ) ){
	if ( verbose ){
	  //#pragma omp critical (verbose)
	  {
	    cerr << "UNK ACRO: " << word << endl;
	  }
	}
	//#pragma omp critical (classify)
	{
	  clean_words[word] += freq;
	  punct_acro_words[word] += freq;
	}
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    //#pragma omp critical (verbose)
	    {
	      cerr << "UNK ACRO: " << word << "/" << acro << endl;
	    }
	  }
	  //#pragma omp critical (classify)
	  {
	    compound_acro_words[acro] += freq;
	  }
	}
      }
      else {
	if ( verbose ){
	  //#pragma omp critical (verbose)
	  {
	    cerr << "UNK word: " << orig_word << endl;
	  }
	}
	//#pragma omp critical (classify)
	{
	  unk_words[orig_word] += freq;
	}
      }
    }
    break;
  case PUNCT:
    {
      set<string> acros;
      if ( doAcro && isAcro( end_pun ) ){
	if ( verbose ){
	  //#pragma omp critical (verbose)
	  {
	    cerr << "PUNCT ACRO: " << end_pun << endl;
	  }
	}
	//#pragma omp critical (classify)
	{
	  punct_acro_words[end_pun] += freq;
	  punct_words[end_pun] = word;
	  clean_words[word] += freq;
	}
      }
      else if ( doAcro && isAcro( parts, acros ) ){
	for ( const auto& acro : acros ){
	  if ( verbose ){
	    //#pragma omp critical (verbose)
	    {
	      cerr << "PUNCT ACRO: (regex) " << word << "/" << acro << endl;
	    }
	  }
	  //#pragma omp critical (classify)
	  {
	    compound_acro_words[acro] += freq;
	  }
	}
      }
      else {
	//#pragma omp critical (classify)
	{
	  clean_words[end_pun] += freq;
	  punct_words[orig_word] = end_pun;
	  if ( verbose ){
	    //#pragma omp critical (verbose)
	    {
	      cerr << "PUNCT word: " << word << endl;
	    }
	  }
	}
      }
    }
    break;
  case UNDEF:
    throw logic_error( "this is realy odd" );
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
				    "ﬀ > ff;"
				    "ﬃ > ffi;"
				    "ﬄ > ffl;"
				    "ﬅ > st;"
				    "ß > ss;"
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
  cerr << "\t--corpus='file'\t validated lexicon file" << endl;
  cerr << "\t--artifrq='value'\t the default value for missing frequencies"
       << endl;
  cerr << "\t\t in the validated lexicon. (default = 0)" << endl;
  cerr << "\t--acro\t also create an acronyms file. (experimental)" << endl;
  cerr << "\t--filter='file'\t use rules from 'file' to transliterate  (experimental)" << endl;
  cerr << "\t\t see http://userguide.icu-project.org/transforms/general/rules for information about rules." << endl;
  cerr << "\t\t default the following filter is used: " << endl
       << default_filter << endl;
#ifdef HAVE_OPENMP
  cerr << "\t-t <threads>\n\t--threads <threads> Number of threads to run on." << endl;
  cerr << "\t\t\t If 'threads' has the value \"max\", the number of threads is set to a" << endl;
  cerr << "\t\t\t reasonable value. (which can be set with OMP_NUM_TREADS environment variable.)" << endl;
#endif
  cerr << "\t-h\t this message " << endl;
  cerr << "\t-v\t be verbose " << endl;
  cerr << "\t-V\t show version " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:t:" );
    opts.set_long_options( "acro,alph:,corpus:,artifrq:,filter:,help,version,threads:" );
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
  string corpusfile;
  unsigned int artifreq = 0;
  bool doAcro = opts.extract("acro");
  verbose = opts.extract('v');
  opts.extract("corpus", corpusfile);
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
  value = "1";
  if ( !opts.extract( 't', value ) ){
    opts.extract( "threads", value );
  }
#ifdef HAVE_OPENMP
  int numThreads = 1;
  if ( TiCC::lowercase(value) == "max" ){
    numThreads = omp_get_max_threads();
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
  string clean_file_name = output_name + ".clean";
  string punct_file_name = output_name + ".punct";
  string acro_file_name = output_name + ".acro";

  if ( !TiCC::createPath( clean_file_name ) ){
    cerr << "unable to open output file: " << clean_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream cs( clean_file_name );
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
    if ( !fillAlpha( as, alphabet ) ){
      cerr << "serious problems reading alphabet file: " << alphafile << endl;
      exit(EXIT_FAILURE);
    }
  }

  map<string,unsigned int> clean_words;
  map<UnicodeString,unsigned int> decap_clean_words;
  map<string,unsigned int> unk_words;
  map<string,unsigned int> punct_acro_words;
  map<string,unsigned int> compound_acro_words;
  map<string,string> punct_words;
  if ( !corpusfile.empty() ){
    if ( artifreq == 0 ){
      cerr << "a corpus file is specified (--corpus option), but artifreq is NOT set "
	   << "(--artifrq option)" << endl;
      exit(EXIT_FAILURE);
    }
    ifstream extra( corpusfile );
    if ( !extra ){
      cerr << "unable to open corpus file: " << corpusfile << endl;
      exit(EXIT_FAILURE);
    }
    string line;
    while ( getline( extra, line ) ){
      vector<string> v;
      int n = TiCC::split_at( line, v, "\t" );
      if ( n == 0 ){
	// empty line, just ignore
	continue;
      }
      if ( n > 2 ){
	cerr << "corpus file in strange format!" << endl;
	cerr << "offending line: " << line << endl;
	exit(EXIT_FAILURE);
      }
      unsigned int freq;
      if ( n == 2 ){
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
      clean_words[v[0]] += freq;
      UnicodeString us = TiCC::UnicodeFromUTF8( v[0] );
      us.toLower();
      decap_clean_words[us] += freq;
    }
  }
  string line;
  size_t line_cnt = 0 ;
  size_t err_cnt = 0;
  map<string,unsigned> my_lexicon;
  while ( getline( is, line ) ){
    ++line_cnt;
    line = TiCC::trim( line );
    if ( line.empty() ){
      continue;
    }
    vector<string> v;
    int n = TiCC::split_at( line, v, "\t" );
    if ( n == 0 ){
      // empty line. just ignore
      continue;
    }
    if ( n < 2 ){
      cerr << "error in line #" << line_cnt
	   << " content='" << line << "'" << endl;
      if ( ++err_cnt > 10 ){
	cerr << "frequency file seems to be in wrong format!" << endl;
	cerr << "too many errors, bailing out" << endl;
	exit(EXIT_FAILURE);
      }
      continue;
    }

    string orig_word = v[0];
    unsigned int freq = TiCC::stringTo<unsigned int>(v[1]);
    my_lexicon[orig_word] = freq;
  }
  cout << "read a lexion with " << line_cnt << " entries"<< endl;
  //#pragma omp parallel for shared(clean_words,unk_words,punct_words,punct_acro_words,compound_acro_words) schedule(dynamic,1)
  for ( size_t i=0; i < my_lexicon.size(); ++i ){
    auto wf = my_lexicon.begin();
    advance( wf, i );
    classify_one_entry( wf->first, wf->second,
			clean_words, decap_clean_words,
			unk_words, punct_words,
			punct_acro_words, compound_acro_words,
			doAcro, alphabet, artifreq );
  }
  cout << "generating output files" << endl;
  cout << "using artifrq=" << artifreq << endl;
  map<unsigned int, set<string> > wf;
  for ( const auto& it : clean_words ){
    wf[it.second].insert( it.first );
  }
  map<unsigned int, set<string> >::const_reverse_iterator wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for ( const auto& sit : wit->second ){
      cs << sit << "\t" << wit->first << endl;
    }
    ++wit;
  }
  cout << "created " << clean_file_name << endl;
  wf.clear();
  for ( const auto& uit : unk_words ){
    wf[uit.second].insert( uit.first );
  }
  wit = wf.rbegin();
  while ( wit != wf.rend() ){
    for ( const auto& sit : wit->second ){
      us << sit << "\t" << wit->first << endl;
    }
    ++wit;
  }
  cout << "created " << unk_file_name << endl;

  if ( doAcro ){
    for ( const auto& ait : punct_acro_words ){
      UnicodeString us = TiCC::UnicodeFromUTF8(ait.first);
      us = filter_punct( us );
      string acro = TiCC::UnicodeToUTF8( us );
      if ( compound_acro_words.find( acro ) != compound_acro_words.end() ){
	// the 'dotted' word is a true acronym
	// add to the list
	compound_acro_words[ait.first] += ait.second;
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
  for ( const auto pit : punct_words ){
    ps << pit.first << "\t" << pit.second << endl;
  }
  cout << "created " << punct_file_name << endl;
  cout << "done!" << endl;
}
