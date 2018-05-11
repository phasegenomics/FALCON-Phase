#include "dna.h"

char * reverse_comp(char * str){

      int i = strlen(str)-1,j=0;
      char ch;
      while (i>j) {
          ch = str[i];
          str[i]= str[j];
          str[j] = basemap[(int)ch];
          i--;
          j++;
      }
      return str;

}
