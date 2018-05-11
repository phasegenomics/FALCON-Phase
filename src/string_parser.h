#ifndef STRING_PARSER_H
#define STRING_PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


 /**
  * [get_next_word description]
  *  This function requires you manage s & e
  * @param  str   the string you want to get the next word
  * @param  s     start index
  * @param  e     end index
  * @param  delim char
  * @return       the next word or null on failture, need to call free
  */
  char * get_next_word(char * str, int *s, int *e, const char delim);

  /**
   * [scan_s description]
   * @param  str      is a pointer to the beginning of your string. More specifically, it is a pointer to the position first character in some string that you want to parse.
   * @param  str_len  is the length of that string. This does not have to be the full string, just the bit you care about.
   * @param  s        points to the start of your token
   * @param  e        points to the end of your token
   * @param  delim    is the character that you want to delmit your string
   * @return          -1 end of string -2 something is wrong
   *
   *
   * Thanks Ryan Layer :
   *
   * http://layerlab.org/2018/01/30/how-i-tokenize-a-string-in-c.html
   */
   int scan_s(char *str, int str_len, int *s, int *e, const char delim);

#endif /* STRING_PARSER_H */
