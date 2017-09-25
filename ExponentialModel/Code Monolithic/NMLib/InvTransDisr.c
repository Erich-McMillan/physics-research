//
// InvTransDisr.c
// Erich Mcmillan
// 15 June 2016
//

#include <stdio.h>
#include <math.h>

#include "ran2.h"

/* global definitions */
float GAMMA, LOGSUM;
int SEQSIZE, I;

void InvTransDisr(long *_seed, int *_seq) {
  /* to use this function the following functions must be set in this order :setSEQSIZE();
  setGAMMA(); setLOGSUM(); After this is done this sequence may be used. InvTransDisr() will
  generate a power-law sequence with a gamma value and size as specified. The largest possible degree will be
  SEQSIZE. The area under the sequence histogram will be 1.
  ****** _seq must be allocated for enough space to fit SEQSIZE ints
  */

  /* declare variables */
  double logprob, logcumul, logcumultmp, rantemp;
  int curNode, n;
  /* log (random number + 1) */
  for(n = 0; n < SEQSIZE; n++) {
    rantemp = ran2(_seed);
    logprob = log1p(rantemp);

    /* set preconditions */
    logcumul = 0.0;
    curNode = 0;

    /* find degree of nth node */
    do {
      logcumultmp = GAMMA*log(++curNode)-LOGSUM;
      logcumul = fmax(logcumul,logcumultmp) + log1p(exp(-fabs(logcumul-logcumultmp)));
    } while (logprob>logcumul);

    /* store value of nth node */
    _seq[n] = curNode;

  }
  return;
}
void setLOGSUM() {
  /* calculates logsum used in InvTransDisr needs only be set before InvTransDisr
     is used and whenever gamma value or sequence size changes */
  int n;
  float logsum = 0.0, logsumtmp;

  for (n=2; n<SEQSIZE; n++) {
    logsumtmp = GAMMA*log(n);
    logsum = fmax(logsum,logsumtmp) + log1p(exp(-fabs(logsum-logsumtmp)));
  }

  LOGSUM = logsum;
  return;
}
void setSEQSIZE(int _seqSize) {
  /* sets SEQSIZE. Needs to be set before any other function is used besides setGAMMA()
  */
  SEQSIZE = _seqSize;
  return;
}
void setGAMMA(float _gamma) {
  /* sets GAMMA. Needs to be set before any other function is used besides setSEQSIZE()
  */
  GAMMA = _gamma;
  return;
}
int ErdosGallai(int *_seq, int _size) {
    /*
        Description: Test the graphicality of the input sequence according to the Erdos Gallai theorem. For this the sequence must be reordered in decending order. According to the Erdos-Gallai theroem the total number of edges must be even and must conform to the test of the theorem which can be found here: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem
        Inputs:
            _seq: the sequence to be checked for graphicality
            _size: the number of elements in the sequence
        Output:
          returns 1 if graphical
          returns 0 if non-graphical
    */

    int a, i, j, k, s, m, sum1, sum2, sumdeg, output;
    sumdeg = 0;
    output = 1;

    /* sort array in decending order using bubble sort */
    for( i = 0; i < _size; i++ ) {
        for( j = i + 1; j < _size; j++ ) {
            if(_seq[i] < _seq[j]) {
                a = _seq[i];
                _seq[i] = _seq[j];
                _seq[j] = a;
            }
        }
    }

    /* test that total degree is even */
    for( k = 0; k < _size; k++ ) {
        sumdeg += _seq[k];
    }
    if(sumdeg % 2 != 0) {
        return 0;
    }

    /* test graphicality of array using erdos gallai theorem */
    for( k = 1; k < _size; k++ ) {
        if(output == 1) {
            sum1 = 0;
            sum2 = 0;
            for( s = 0; s < k; s++ ) {
                sum1 += _seq[s];
            }
            for( m = k; m < _size; m++ ) {
                //printf("m is: %d\n", m );
                if(k > _seq[m]) {
                    sum2 += _seq[m];
                } else {
                    sum2 += k;
                }
            }

            //printf("sum1 is: %d	sum2 is: %d\n", sum1, sum2);
            if(sum1 > ((k * (k - 1)) + sum2)) {
                output = 0;
            }
        }
    }

    //printf("is Graphical: %d\n", output);
    return output;
}
int numEdgesEven(int _size, int *_seq) {
  // returns 0 if the value is odd returns 1 if even
  int i, a, j, sum=0;

  /* sort array in decending order using bubble sort */
  for( i = 0; i < _size; i++ ) {
      for( j = i + 1; j < _size; j++ ) {
          if(_seq[i] < _seq[j]) {
              a = _seq[i];
              _seq[i] = _seq[j];
              _seq[j] = a;
          }
      }
  }

  for(i = 0; i < _size; i++) {
    sum += _seq[i];
  }

  return (0 && (sum%2));
}
