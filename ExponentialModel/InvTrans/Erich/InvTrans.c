//
//  InvTrans.c
//  ExponentialGraphModel
//
//  Created by erich on 10/14/15.
//  Copyright (c) 2015 erich. All rights reserved.
//

// MARK: Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "InvTrans.h"
#include <math.h>
#include "ran2.h"

// MARK: Functions
float InvTrans(int _n, double _g, long *_p, int *_seq) {
    /*
        Description: Generates n random numbers m times with a power law distribution where p(x) = cx^gamma using a normally distributed random number generator and returns the degree sequence of the network. The inputs of the function are: n:the number of elements in the sequence, g:the exponent of the distribution, and p: a pointer of an negative integer used to seed the ran2 random number generator. *seq must be a pointer to an array 1*n integers long. The probability/area under the output curve will be ~1. The sequence will also be graphical according to the Erdos_Gallai theorem.
        Inputs:
            _n: the size and max value of the desired sequence
            _g: the exponent gamma of the desired sequence
            _p: a pointer used to seed the random number generator
            _seq: where the generated sequence is stored for use in the program and for output
        Outputs:
            _seq: this is where the generated sequence will reside after the function has completed.
            C: this value is directly returned from the function which indicates the c in the model necessary to make the area under the curve 1 for this sequence size and gamma value.
    */

    /* declare variables */
    const double g1 = _g+1;
    double C, temp;

    C = 1.0/((pow((double)_n-.5, g1)/g1)-(pow(.5, g1)/g1));

    /* generate sequences until one which is graphical is found */
    do {
        for(int i = 0; i < _n; i++) {
            temp = pow((ran2(_p)/C)*(g1)+pow(.5, g1),1.0/g1)+.5;
            _seq[i] = temp;
        }
    } while(!ErdosGallai(_seq, _n));

    return C;
}
float InvTransMaxIts(int _n, double _g, long *_p, int *_seq, unsigned int _maxits) {
    /*
     Description: Generates n random numbers m times with a power law distribution where p(x) = cx^gamma using a normally distributed random number generator and returns the degree sequence of the network. The inputs of the function are: n:the number of elements in the sequence, g:the exponent of the distribution, and p: a pointer of an negative integer used to seed the ran2 random number generator. *seq must be a pointer to an array 1*n integers long. The probability/area under the output curve will be ~1. The sequence will also be graphical according to the Erdos_Gallai theorem.
     Inputs:
     _n: the size and max value of the desired sequence
     _g: the exponent gamma of the desired sequence
     _p: a pointer used to seed the random number generator
     _seq: where the generated sequence is stored for use in the program and for output
     Outputs:
     _seq: this is where the generated sequence will reside after the function has completed.
     C: this value is directly returned from the function which indicates the c in the model necessary to make the area under the curve 1 for this sequence size and gamma value.
     returns 0 if the maxits is exceeded this means no suitable sequence was found after checking _maxits sequences
     */

    /* declare variables */
    const double g1 = _g+1;
    double C, temp;
    unsigned int its = 0;

    C = 1.0 / ( (pow( (double)_n-.5, g1 ) / g1) - (pow(.5, g1)/g1));
    /* generate sequences until one which is graphical is found */
    do {
        for(int i = 0; i < _n; i++) {
            temp = pow((ran2(_p)/C)*(g1)+pow(.5, g1), 1.0/g1)+.5;
            _seq[i] = temp;
        }
        its++;
    } while(!ErdosGallai(_seq, _n) && its < _maxits);

    if(its >= _maxits) return 0;
    else return C;
}
float InvTransNoCheck(int _n, double _g, long *_p, int *_seq) {
    /*
        Description: performs the same aso InvTrans() but does not check if the sequence is graphical according to erdos gallai theorem
    */

    const double g1 = _g+1;
    double C, temp;

    C = 1.0/((pow((double)_n-.5, g1)/g1)-(pow(.5, g1)/g1));
    // generate a sequence
    for(int i = 0; i < _n; i++) {
        temp = pow((ran2(_p)/C)*(g1)+pow(.5, g1),1.0/g1) +.5;
        _seq[i] = temp;
    }

    return C;
}
int ErdosGallai(int *_seq, int _size) {
    /*
        Description: Test the graphicality of the input sequence according to the Erdos Gallai theorem. For this the sequence must be reordered in decending order. According to the Erdos-Gallai theroem the total number of edges must be even and must conform to the test of the theorem which can be found here: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem
        Inputs:
            _seq: the sequence to be checked for graphicality
            _size: the number of elements in the sequence
    */

    int a, sum1, sum2, sumdeg, output;
    sumdeg = 0;
    output = 1;

    /* sort array in decending order using bubble sort */
    for( int i = 0; i < _size; i++ ) {
        for( int j = i + 1; j < _size; j++ ) {
            if(_seq[i] < _seq[j]) {
                a = _seq[i];
                _seq[i] = _seq[j];
                _seq[j] = a;
            }
        }
    }

    /* test that total degree is even */
    for( int k = 0; k < _size; k++ ) {
        sumdeg += _seq[k];
    }
    if(sumdeg % 2 != 0) {
        return 0;
    }

    /* test graphicality of array using erdos gallai theorem */
    for( int k = 1; k < _size; k++ ) {
        if(output == 1) {
            sum1 = 0;
            sum2 = 0;
            for( int s = 0; s < k; s++ ) {
                sum1 += _seq[s];
            }
            for( int m = k; m < _size; m++ ) {
                //printf("m is: %d\n", m );
                if(k > _seq[m]) {
                    sum2 += _seq[m];
                } else {
                    sum2 += k;
                }
            }

            //printf("sum1 is: %d\t sum2 is: %d\n", sum1, sum2);
            if(sum1 > ((k * (k - 1)) + sum2)) {
                output = 0;
            }
        }
    }

    //printf("is Graphical: %d\n", output);
    return output;
}
