#include "string_parser.h"

 int scan_s(char *str, int str_len, int *s, int *e, const char delim) {
     if (*e == str_len)
         return -1;
     for (*e = *s; *e <= str_len; *e+=1) {
         if ((str[*e] == delim) || (*e == str_len)) {
             return *e - *s;
         }
     }
     return -2;
 }

 char * get_next_word(char * str, int *s, int *e, const char delim){

   char * word = NULL;

   int len  = strlen(str);
   int flag = scan_s(str, len, s, e, delim);

   if(flag < 0) return word;

   word = (char *) malloc((*e-*s+1)*sizeof(char));


   char * endptr = stpncpy(word, str + *s, *e-*s) ;
   word[*e-*s] = '\0';
   *s = *e + 1;

   return word;
 }
