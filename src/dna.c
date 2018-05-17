/*
Author: Zev Kronenberg
Contact :zev@phasegenomics.com
Date: May 17th 2018

The Clear BSD + Attribution License

Copyright (c) 2018, Pacific Biosciences of California, Inc. and Phase Genomics, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are 
permitted (subject to the limitations in the disclaimer below) provided that the 
following conditions are met:

1.Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2.Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3.All advertising materials mentioning features or use of this software 
must display the following acknowledgement:
  This <product/service> includes <software/the use of software> developed 
  by Pacific Biosciences of California, Inc. and Phase Genomics, Inc.

  4.Distributions of data generated through the use of this software as 
part of a service provided by a for-profit organization must be accompanied 
by the above copyright notice, this list of conditions, the following 
acknowledgement and the disclaimer below:
  This data was generated using software developed by Pacific Biosciences of 
  California, Inc. and Phase Genomics, Inc.

  5.Neither the names of the copyright holders nor the names of their 
contributors may be used to endorse or promote products derived from this 
software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY’S PATENT RIGHTS ARE GRANTED BY 
THIS LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS 
OR THEIR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/





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
