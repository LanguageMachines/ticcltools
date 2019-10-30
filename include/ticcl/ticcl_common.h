#ifndef TICCL_COMMON_H
#define TICCL_COMMON_H

#include "unicode/unistr.h"
#include "unicode/ustream.h"
#include "unicode/uchar.h"

typedef uint64_t bitType;

const std::string S_SEPARATOR = "_";
const icu::UnicodeString US_SEPARATOR = "_";
const UChar U_SEPARATOR = '_';


bitType high_five( int val );

inline std::string toString( int8_t c ){
  switch ( c ){
  case 0:
    return "U_UNASSIGNED";
  case 1:
    return "U_UPPERCASE_LETTER";
  case 2:
    return "U_LOWERCASE_LETTER";
  case 3:
    return "U_TITLECASE_LETTER";
  case 4:
    return "U_MODIFIER_LETTER";
  case 5:
    return "U_OTHER_LETTER";
  case 6:
    return "U_NON_SPACING_MARK";
  case 7:
    return "U_ENCLOSING_MARK";
  case 8:
    return "U_COMBINING_SPACING_MARK";
  case 9:
    return "U_DECIMAL_DIGIT_NUMBER";
  case 10:
    return "U_LETTER_NUMBER";
  case 11:
    return "U_OTHER_NUMBER";
  case 12:
    return "U_SPACE_SEPARATOR";
  case 13:
    return "U_LINE_SEPARATOR";
  case 14:
    return "U_PARAGRAPH_SEPARATOR";
  case 15:
    return "U_CONTROL_CHAR";
  case 16:
    return "U_FORMAT_CHAR";
  case 17:
    return "U_PRIVATE_USE_CHAR";
  case 18:
    return "U_SURROGATE";
  case 19:
    return "U_DASH_PUNCTUATION";
  case 20:
    return "U_START_PUNCTUATION";
  case 21:
    return "U_END_PUNCTUATION";
  case 22:
    return "U_CONNECTOR_PUNCTUATION";
  case 23:
    return "U_OTHER_PUNCTUATION";
  case 24:
    return "U_MATH_SYMBOL";
  case 25:
    return "U_CURRENCY_SYMBOL";
  case 26:
    return "U_MODIFIER_SYMBOL";
  case 27:
    return "U_OTHER_SYMBOL";
  case 28:
    return "U_INITIAL_PUNCTUATION";
  case 29:
    return "U_FINAL_PUNCTUATION";
  default:
    return "OMG IK HEB GEEN IDEE!";
  }
}

inline bool ticc_ispunct( int8_t charT ){
  return ( charT == U_OTHER_PUNCTUATION ||
	   charT == U_INITIAL_PUNCTUATION ||
	   charT == U_FINAL_PUNCTUATION ||
	   charT == U_CONNECTOR_PUNCTUATION ||
	   charT == U_START_PUNCTUATION ||
	   charT == U_END_PUNCTUATION ||
	   charT == U_DASH_PUNCTUATION );
}

inline bool ticc_ispunct( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_ispunct( charT );
}

inline bool ticc_isdigit( int8_t charT ){
  return charT == U_DECIMAL_DIGIT_NUMBER;
}

inline bool ticc_isdigit( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_isdigit( charT );
}

inline bool ticc_isupper( int8_t charT ){
  return charT == U_UPPERCASE_LETTER;
}

inline bool ticc_isupper( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_isupper( charT );
}

inline bool ticc_islower( int8_t charT ){
  return charT == U_LOWERCASE_LETTER;
}

inline bool ticc_islower( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_islower( charT );
}

inline bool ticc_isother( int8_t charT ){
  return ( charT == U_MATH_SYMBOL ||
	   charT == U_OTHER_SYMBOL ||
	   charT == U_OTHER_NUMBER ||
	   charT == U_FORMAT_CHAR ||
	   charT == U_NON_SPACING_MARK ||
	   charT == U_CURRENCY_SYMBOL ||
	   charT == U_MODIFIER_SYMBOL ||
	   charT == U_CONTROL_CHAR );
}

inline bool ticc_isother( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_isother( charT );
}

inline bool ticc_isletter( int8_t charT ){
  return ( charT == U_LOWERCASE_LETTER ||
	   charT == U_UPPERCASE_LETTER );
}

inline bool ticc_isletter( UChar uc ){
  int8_t charT =  u_charType( uc );
  return ticc_isletter( charT );
}

#endif //  TICCL_COMMON_H
