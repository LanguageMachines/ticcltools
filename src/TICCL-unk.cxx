#include <cstdlib>
#include <getopt.h>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "ticcutils/CommandLine.h"
#include "ticcutils/StringOps.h"
#include "ticcl/unicode.h"

#include "config.h"

using namespace	std;
//#define DEBUG

bool verbose = false;

enum S_Class { UNK, PUNCT, IGNORE, CLEAN };
ostream& operator<<( ostream& os, const S_Class& cl ){
  switch ( cl ){
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
    UnicodeString us = UTF8ToUnicode( v[0] );
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
    UnicodeString us = UTF8ToUnicode( line );
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
  int is_digit = 0;
  int is_punct = 0;
  int is_in = 0;
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
      if ( alphabet.empty() ){
	if ( verbose ){
	  cerr << "bekijk karakter " << UnicodeString(uchar) << " van type " << toString(charT) << endl;
	}
	if ( ticc_isletter( charT ) ){
	  ++is_in;
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
	  ++is_in;
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
    cerr << "Classify: " << word << " IN=" << is_in << " OUT= " << is_out << " DIG=" << is_digit << " PUNCT=" << is_punct << endl;
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
  else if ( is_digit >= is_in + is_out ){
    // Filter B: dingen als datums, floats, of combinatie getal + een of andere
    // geldaanduiding : zelfde als getallen
    // <martin> Komt erop neer dat indien meer cijfers dan iets anders.
    if ( verbose ){
      cerr << "UITGANG 2: Ignore" << endl;
    }
    return IGNORE;
  }
  else if ( word_len >= 4 && double(is_in + is_digit)/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 3: Clean" << endl;
    }
    return CLEAN;
  }
  else if (word_len == 3 && double(is_in + is_digit)/word_len >= 0.66 ){
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
  else if ( word_len >= 4 && double( is_in )/word_len >= 0.75 ){
    if ( verbose ){
      cerr << "UITGANG 6: Clean" << endl;
    }
    return CLEAN;
  }
  else if ( word_len == 3 && double( is_in )/word_len >= 0.66 ){
    if ( verbose ){
      cerr << "UITGANG 7: Clean" << endl;
    }
    return CLEAN;
  }
  else {
    if ( verbose ){
      cerr << "UITGANG 8: UNK" << endl;
    }
    return UNK;
  }
}

S_Class classify( const string& word, set<UChar>& alphabet,
		  string& punct ){
  S_Class result = CLEAN;
  punct.clear();
  UnicodeString us = UTF8ToUnicode( word );
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
	  punct = UnicodeToUTF8( ps );
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

bool isAcro( const string& word, string& stripped ){
  stripped = "";
  UnicodeString us = UTF8ToUnicode( word );
  if ( u_ispunct( us[0] ) ){
    us = UnicodeString( us, 1 );
    stripped = UnicodeToUTF8(us);
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
    opts.set_short_options( "vVh" );
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
  ifstream is( file_name.c_str() );
  if ( !is ){
    cerr << "unable to find or open frequency file: " << file_name << endl;
    exit(EXIT_FAILURE);
  }
  set<UChar> alphabet;

  if ( !alphafile.empty() ){
    ifstream as( alphafile.c_str() );
    if ( !as ){
      cerr << "unable to open alphabet file: " << alphafile << endl;
      exit(EXIT_FAILURE);
    }
    if ( !fillAlpha( as, alphabet ) ){
      cerr << "serious problems reading alphabet file: " << alphafile << endl;
      exit(EXIT_FAILURE);
    }
  }

  string unk_file_name = file_name + ".unk";
  string clean_file_name = file_name + ".clean";
  string punct_file_name = file_name + ".punct";
  string acro_file_name = file_name + ".acro";

  ofstream cs( clean_file_name.c_str() );
  if ( !cs ){
    cerr << "unable to open output file: " << clean_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream us( unk_file_name.c_str() );
  if ( !us ){
    cerr << "unable to open output file: " << unk_file_name << endl;
    exit(EXIT_FAILURE);
  }
  ofstream ps( punct_file_name.c_str() );
  if ( !ps ){
    cerr << "unable to open output file: " << punct_file_name << endl;
    exit(EXIT_FAILURE);
  }
  if ( doAcro ){
    ofstream as( acro_file_name.c_str() );
    if ( !as ){
      cerr << "unable to open output file: " << acro_file_name << endl;
      exit(EXIT_FAILURE);
    }
  }

  map<string,unsigned int> clean_words;
  map<string,unsigned int> unk_words;
  map<string,unsigned int> acro_words;
  map<string,string> punct_words;
  if ( !corpusfile.empty() ){
    if ( artifreq == 0 ){
      cerr << "a corpus file is specified (--corpus option), but artifreq is NOT set "
	   << "(--artifrq option)" << endl;
      exit(EXIT_FAILURE);
    }
    ifstream extra( corpusfile.c_str() );
    if ( !extra ){
      cerr << "unable to open corpus file: " << corpusfile << endl;
      exit(EXIT_FAILURE);
    }
    string line;
    while ( getline( extra, line ) ){
      vector<string> v;
      int n = TiCC::split_at( line, v, "\t" );
      if ( n < 1 || n > 2 ){
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
    }
  }
  string line;
  while ( getline( is, line ) ){
    line = TiCC::trim( line );
    if ( line.empty() ){
      continue;
    }
    vector<string> v;
    int n = TiCC::split_at( line, v, "\t" );
    if ( n < 2 ){
      cerr << "frequency file in wrong format!" << endl;
      cerr << "offending line: " << line << endl;
      exit(EXIT_FAILURE);
    }
    unsigned int freq = TiCC::stringTo<unsigned int>(v[1]);

    string pun;
    string word = v[0];
    S_Class cl;
    if ( clean_words.find( word ) != clean_words.end() ){
      // no need to do a lot of work for already clean words
      cl = CLEAN;
    }
    else {
      cl = classify( word, alphabet, pun );
    }
    switch ( cl ){
    case IGNORE:
      break;
    case CLEAN:
      {
	clean_words[word] += freq;
	string stripped;
	if ( doAcro && isAcro( word, stripped ) ){
	  if ( verbose ){
	    cerr << "CLEAN ACRO: " << word << "/" << stripped << endl;
	  }
	  acro_words[word] += 1;
	}
	else if ( verbose ){
	  cerr << "CLEAN word: " << word << endl;
	}
      }
      break;
    case UNK:
      {
	string stripped;
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
	punct_words[word] = pun;
	clean_words[pun] += freq;
	string punc = pun + ".";
	string stripped;
	if ( doAcro && isAcro( punc,stripped ) ){
	  if ( verbose ){
	    cerr << "PUNCT ACRO: " << word << "/" << pun << "/" << stripped << endl;
	  }
	  acro_words[pun] += 1;
	}
	else if ( verbose ){
	  cerr << "PUNCT word: " << word << endl;
	}
      }
      break;
    }
  }
  cout << "generating output files" << endl;
  cout << "using artifrq=" << artifreq << endl;
  map<unsigned int, set<string> > wf;
  map<string,unsigned int >::const_iterator it = clean_words.begin();
  while( it != clean_words.end()  ){
    wf[it->second].insert( it->first );
    ++it;
  }
  map<unsigned int, set<string> >::const_reverse_iterator wit = wf.rbegin();
  while ( wit != wf.rend() ){
    set<string>::const_iterator sit = wit->second.begin();
    while ( sit != wit->second.end() ){
      cs << *sit << "\t" << wit->first << endl;
      ++sit;
    }
    ++wit;
  }
  cout << "created " << clean_file_name << endl;
  wf.clear();
  it = unk_words.begin();
  while( it != unk_words.end()  ){
    wf[it->second].insert( it->first );
    ++it;
  }
  wit = wf.rbegin();
  while ( wit != wf.rend() ){
    set<string>::const_iterator sit = wit->second.begin();
    while ( sit != wit->second.end() ){
      us << *sit << "\t" << wit->first << endl;
      ++sit;
    }
    ++wit;
  }
  cout << "created " << unk_file_name << endl;
  map<string,string>::const_iterator it2 = punct_words.begin();
  while ( it2 != punct_words.end() ){
    ps << it2->first << "\t" << it2->second << endl;
    ++it2;
  }
  cout << "created " << punct_file_name << endl;

  if ( doAcro ){
    ofstream as( acro_file_name.c_str() );
    map<string,unsigned int>::const_iterator it3 = acro_words.begin();
    while ( it3 != acro_words.end() ){
      as << it3->first << "\t" << it3->second << endl;
      ++it3;
    }
    cout << "created " << acro_file_name << endl;
  }

  cout << "done!" << endl;
}
