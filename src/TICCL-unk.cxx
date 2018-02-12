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
#include "ticcutils/Unicode.h"
#include "ticcl/unicode.h"

#include "config.h"

using namespace	std;
//#define DEBUG

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

bool fillSimpleAlpha( istream& is, set<UChar>& alphabet ){
  string line;
  while ( getline( is, line ) ){
    if ( line.size() == 0 || line[0] == '#' ){
      continue;
    }
    vector<string> vec;
    TiCC::split( line, vec );
    line = vec[0];
    UnicodeString us = TiCC::UnicodeFromUTF8( line );
    us.toLower();
    for( int i=0; i < us.length(); ++i )
      alphabet.insert( us[i] );
    us.toUpper();
    for( int i=0; i < us.length(); ++i )
      alphabet.insert( us[i] );
  }
  return true;
}

bool depunct( const UnicodeString& us, UnicodeString& result ){
  result.remove();
  int i = 0;
  for ( ; i < us.length(); ++i ){
    // skip leading punctuation and spaces
    if ( !( u_ispunct( us[i] ) || u_isspace( us[i] ) ) )
      break;
  }
  int j = us.length()-1;
  for ( ; j >= 0; j-- ){
    // skip trailing punctuation and spaces
    if ( !( u_ispunct( us[j] ) || u_isspace( us[j] ) ) )
      break;
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


S_Class classify( const UnicodeString& word,
		  set<UChar>& alphabet ){
  int is_upper = 0;
  int is_digit = 0;
  int is_punct = 0;
  int is_letter = 0;
  int is_out = 0;
  int is_space = 0;
  int word_len = word.length();
  if ( word_len < 2 ){
#ifdef DEBUG
    cerr << "UITGANG 0: Ignore" << endl;
#endif
    return IGNORE;
  }
  for ( int i=0; i < word_len; ++i ){
    UChar uchar = word[i];
    if ( u_isspace( uchar ) ){
      ++is_space;
    }
    else {
      int8_t charT = u_charType( uchar );
      if ( ticc_isupper( charT ) ){
	++is_upper;
      }
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
	  cerr << "Warning: karakter " << UnicodeString(uchar) << " is van onbekend type " << toString(charT) << endl;
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
	else {
	  if ( verbose ){
	    cerr << "'" << UnicodeString(uchar) << "' is OUT het alfabet" << endl;
	  }
	  ++is_out;
	}
      }
    }
  }
  //  word_len -= is_space;
  if ( verbose ){
    cerr << "Classify: " << word << " IN=" << is_letter << " OUT= " << is_out << " DIG=" << is_digit << " PUNCT=" << is_punct << endl;
  }
  if ( is_digit == word_len ){
    // Filter A: gewone getallen. Worden ongemoeid gelaten, worden dus niet
    // ge-unkt of ge-ticcled, worden ook niet geteld of opgenomen in de
    // frequentielijst
    if ( verbose ){
      cerr << "UITGANG 1: Ignore" << endl;
    }
    return IGNORE;
  }
  else if ( is_digit >= is_letter + is_out ){
    // Filter B: dingen als datums, floats, of combinatie getal + een of andere
    // geldaanduiding : zelfde als getallen
    // <martin> Komt erop neer dat indien meer cijfers dan iets anders.
    if ( verbose ){
      cerr << "UITGANG 2: Ignore" << endl;
    }
    return IGNORE;
  }
  else if ( word_len >= 4 && double(is_letter + is_digit)/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 3: Clean" << endl;
    }
    return CLEAN;
  }
  else if (word_len == 3 && double(is_letter + is_digit)/word_len >= 0.66 ){
    if ( verbose ){
      cerr << "UITGANG 4: Clean" << endl;
    }
    return CLEAN;
  }
  else if (word_len < 3 && (is_out+is_digit) <1 ){
    if ( verbose ){
      cerr << "UITGANG 5: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( word_len >= 4 && double( is_letter )/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 6: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( word_len == 3 && double( is_letter )/word_len >= 0.66 ){
    if ( verbose ){
      cerr << "UITGANG 7: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( ( word_len >= 8
	      && (is_letter + is_punct) > (word_len - 2) )
	    || ( word_len >= 4
		 && (is_letter + is_punct) > (word_len - 1) ) ){
    if ( verbose ){
      cerr << "UITGANG 8: Clean" << endl;
    }
    return CLEAN;
  }
  else {
    if ( verbose ){
      cerr << "UITGANG 9: UNK" << endl;
    }
    return UNK;
  }
}

S_Class classify( const string& word, set<UChar>& alphabet,
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
      if ( result != IGNORE ){
	if ( result == CLEAN ){
	  punct = TiCC::UnicodeToUTF8( ps );
	  result = PUNCT;
	}
	else
	  result = UNK;
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

bool isAcro( const vector<string>& parts, string& result ){
  static string pattern =  "(?:de|het|een)" + SEPARATOR + "([A-Z]+)-(?:[a-z]*)";
  static TiCC::UnicodeRegexMatcher acro_detect( TiCC::UnicodeFromUTF8(pattern), "newline_splitter" );
  if ( parts.size() != 2 ){
    return false;
  }
  UnicodeString us = TiCC::UnicodeFromUTF8( parts[0] + SEPARATOR + parts[1] );
  UnicodeString pre, post;
  //  acro_detect.set_debug(1);
  //  cerr << "IS ACRO: test pattern = " << acro_detect.Pattern() << endl;
  //  cerr << "op " << us << endl;
  if ( acro_detect.match_all( us, pre, post ) ){
    //    cerr << "IT Mached!" << endl;
    result = TiCC::UnicodeToUTF8(acro_detect.get_match( 0 ));
    cerr << "FOUND regexp acronym: " << result << endl;
    return true;
  }
  return false;
}

bool isAcro( const string& word, string& stripped ){
  stripped = "";
  UnicodeString us = TiCC::UnicodeFromUTF8( word );
  if ( u_ispunct( us[0] ) ){
    us = UnicodeString( us, 1 );
    stripped = TiCC::UnicodeToUTF8(us);
  }
  if ( us.length() < 6 ){
    UnicodeString Us = us;
    if ( Us.toUpper() == us ){
      return true;
    }
  }
  bool isOK = true;
  for ( int i=0; i < us.length(); ++i ){
    if ( u_ispunct( us[i] ) ){
      if ( isOK ){
	return false;
      }
      else {
	isOK = true;
      }
    }
    else if ( u_isupper( us[i] ) && isOK ){
      isOK = false;
    }
    else {
      return false;
    }
  };
  return true;
}

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
  cerr << "\t-h\t this message " << endl;
  cerr << "\t-v\t be verbose " << endl;
  cerr << "\t-V\t show version " << endl;
}

int main( int argc, char *argv[] ){
  TiCC::CL_Options opts;
  try {
    opts.set_short_options( "vVho:" );
    opts.set_long_options( "acro,alph:,corpus:,artifrq:" );
    opts.parse_args( argc, argv );
  }
  catch( TiCC::OptionError& e ){
    cerr << e.what() << endl;
    usage( argv[0] );
    exit( EXIT_FAILURE );
  }
  string progname = opts.prog_name();
  if ( opts.extract('V' ) ){
    cerr << PACKAGE_STRING << endl;
    exit(EXIT_SUCCESS);
  }
  if ( opts.extract('h' ) ){
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

  ofstream cs( clean_file_name );
  if ( !cs ){
    cerr << "unable to open output file: " << clean_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream us( unk_file_name );
  if ( !us ){
    cerr << "unable to open output file: " << unk_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream ps( punct_file_name );
  if ( !ps ){
    cerr << "unable to open output file: " << punct_file_name << endl;
    exit(EXIT_FAILURE);
  }
  if ( doAcro ){
    ofstream as( acro_file_name );
    if ( !as ){
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
  map<string,unsigned int> acro_words;
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
    unsigned int freq = TiCC::stringTo<unsigned int>(v[1]);

    S_Class end_cl = UNDEF;

    string end_pun;

    string word = v[0];
    vector<string> parts;
    TiCC::split_at( word, parts, SEPARATOR );
    if ( parts.size() == 0 ){
      continue;
    }
    unsigned int lexclean = 0;
    for ( auto const& wrd : parts ){
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
      switch( cl ){
      case IGNORE:
	if ( end_cl == UNDEF ){
	  end_cl = IGNORE;
	}
	else if ( end_cl == IGNORE ){
	}
	else if ( end_cl == CLEAN ){
	}
	else if ( end_cl == UNK ){
	}
	else if ( end_cl == PUNCT ){
	}
	break;
      case CLEAN:
	if ( end_cl == UNDEF ){
	  end_cl = CLEAN;
	}
	else if ( end_cl == IGNORE ){
	  end_cl = CLEAN;
	}
	else if ( end_cl == CLEAN ){
	}
	else if ( end_cl == UNK ){
	}
	else if ( end_cl == PUNCT ){
	}
	break;
      case PUNCT:
	if ( end_cl == UNDEF ){
	  end_cl = PUNCT;
	}
	else if ( end_cl == IGNORE ){
	  end_cl = PUNCT;
	}
	else if ( end_cl == CLEAN ){
	  end_cl = PUNCT;
	}
	else if ( end_cl == UNK ){
	}
	else if ( end_cl == PUNCT ){
	}
	break;
      case UNK:
	if ( end_cl == UNDEF ){
	  end_cl = UNK;
	}
	else if ( end_cl == IGNORE ){
	  end_cl = UNK;
	}
	else if ( end_cl == CLEAN ){
	  end_cl = UNK;
	}
	else if ( end_cl == UNK ){
	}
	else if ( end_cl == PUNCT ){
	  end_cl = UNK;
	}
	break;
      case UNDEF:
	throw logic_error( "undef value returned by classify()" );
      }
      if ( pun.empty() ){
	pun = wrd;
      }
      end_pun += pun + SEPARATOR;
    }
    end_pun = TiCC::trim_back( end_pun, SEPARATOR );

    switch ( end_cl ){
    case IGNORE:
      break;
    case CLEAN:
      {
	clean_words[word] += freq;
	if ( lexclean == parts.size() ){
	  clean_words[word] += artifreq;
	}
	string stripped;
	string acro;
	if ( doAcro && isAcro( word, stripped ) ){
	  if ( verbose ){
	    cerr << "CLEAN ACRO: " << word << "/" << stripped << endl;
	  }
	  acro_words[word] += 1;
	}
	else if ( doAcro && isAcro( parts, acro ) ){
	  if ( verbose ){
	    cerr << "CLEAN ACRO: (regex)" << word << "/" << acro << endl;
	  }
	  acro_words[acro] += 1;
	}
	else if ( verbose ){
	  cerr << "CLEAN word: " << word << endl;
	}
      }
      break;
    case UNK:
      {
	string stripped;
	string acro;
	if ( doAcro && isAcro( word, stripped ) ){
	  if ( verbose ){
	    cerr << "UNK ACRO: " << word << "/" << stripped << endl;
	  }
	  if ( !stripped.empty() ){
	    punct_words[word] = stripped;
	    clean_words[stripped] += freq;
	    acro_words[stripped] += 1;
	  }
	  else {
	    clean_words[word] += freq;
	    acro_words[word] += 1;
	  }
	}
	else if ( doAcro && isAcro( parts, acro ) ){
	  if ( verbose ){
	    cerr << "UNK ACRO: " << word << "/" << acro << endl;
	  }
	  acro_words[acro] += 1;
	}
	else {
	  if ( verbose ){
	    cerr << "UNK word: " << word << endl;
	  }
	  unk_words[word] += freq;
	}
      }
      break;
    case PUNCT:
      {
	punct_words[word] = end_pun;
	clean_words[end_pun] += freq;
	string punc = end_pun + ".";
	string stripped;
	string acro;
	if ( doAcro && isAcro( punc,stripped ) ){
	  if ( verbose ){
	    cerr << "PUNCT ACRO: " << word << "/" << end_pun << "/" << stripped << endl;
	  }
	  acro_words[end_pun] += 1;
	}
	else if ( doAcro && isAcro( parts, acro ) ){
	  if ( verbose ){
	    cerr << "PUNCT ACRO: (regex) " << word << "/" << acro << endl;
	  }
	  acro_words[acro] += 1;
	}
	else if ( verbose ){
	  cerr << "PUNCT word: " << word << endl;
	}
      }
      break;
    case UNDEF:
      throw logic_error( "this is realy odd" );
    }
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
  for ( const auto pit : punct_words ){
    ps << pit.first << "\t" << pit.second << endl;
  }
  cout << "created " << punct_file_name << endl;

  if ( doAcro ){
    ofstream as( acro_file_name );
    for ( const auto& ait : acro_words ){
      as << ait.first << "\t" << ait.second << endl;
    }
    cout << "created " << acro_file_name << endl;
  }

  cout << "done!" << endl;
}
