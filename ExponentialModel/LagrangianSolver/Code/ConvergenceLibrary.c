//
//  ConvergenceLibrary.c
//  ExponentialGraphModel
//
//  Created by erich on 4/8/16.
//  Copyright © 2016 erich. All rights reserved.
//

#include "ConvergenceLibrary.h"
#include "nr.h"
#include "nrutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


// MARK: Internal Function Declarations
void resetSEQSEQQUANT();
void resetFVEC();
void resetXVEC();
void setFVEC(float* fvec, int fvecsize);
void setCONVPROB(int errorcode);
int arrayContainsInternal(int, int*, int, int*);
int arrayContainsInternalF(float _num, float *_seq, int _size, int *_index);


// MARK: External Definitions
int *CONVSEQ, *SEQQUANT, MAXSIZE, CURRSIZE, SIMPEXP, NEWTON, CONVPROB, ITS;
float *FVEC, *XVEC;
double TIMEELAPSED;

// MARK: Library Initialization Function
void initializeLibrary(int _maxSequenceSize) {
    /*
        Description: Allocates enough space to fit inputs up to _maxSequenceSize
        Input:
            _maxSequenceSize: the maximum size to allocate space
    */

    MAXSIZE = _maxSequenceSize;
    CURRSIZE = 0;
    CONVSEQ = malloc(sizeof(int)*_maxSequenceSize);
    SEQQUANT = malloc(sizeof(int)*_maxSequenceSize);
    FVEC = vector(1, _maxSequenceSize);
    XVEC = vector(1, _maxSequenceSize);
    SIMPEXP = 0;
    NEWTON = 0;
    CONVPROB = 0;
    return;
}
void resetLibrary(int _newMaxSeqSize) {
    freeLibrary();
    initializeLibrary(_newMaxSeqSize);
    return;
}
void reallocLibrary(int _newMaxSeqSize) {
    /*
        Description: Saves the current settings of the library and then resets the library and sets the settings back
        Inputs:
            _newMaxSeqSize: the new size for which the library is allocated
    */

    // Grab settings of library
    int tempSimp, tempNewt, tempConv;
    tempSimp = SIMPEXP;
    tempNewt = NEWTON;
    tempConv = CONVPROB;

    // reset library
    resetLibrary(_newMaxSeqSize);

    // reset settings of library
    SIMPEXP = tempSimp;
    NEWTON = tempNewt;
    CONVPROB = tempConv;

    return;
}
void freeLibrary() {
    /*
        Description: Frees the memory from all arrays used in this library should be used when one is finished using the library
    */

    free_vector(FVEC, 1, MAXSIZE);
    free_vector(XVEC, 1, MAXSIZE);
    free(CONVSEQ);
    free(SEQQUANT);
    return;
}

// MARK: Systems of Equations Setter
void funcv(int _n, float _x[], float _f[]) {
    /*
         Description: Populates the x[] and f[] arrays with the functions describing the system of equations to be solved by either the Broyden method or Newton method. To do this f[] becomes the function such as y=x1+2(x2) where x1 and x2 are independent variables in the equation set which would be represented as f[equation number] = x[1]+2x[2]; For the simplifed exponential model the form of the equations is simplifed to reduce computation requirements as such repeated roots are represented by changing the coefficents in front of the part of the equation which represents this root value. See example below for simplification of 4 equations:
             For the sequence {1, 2, 1, 1} one is repeated 3 times therefore the equations can be simplified from 4 to 2. This may not seem like much however as sequence sizes become larger and due to the properites of power-law distributed networks where a great percentage of the roots are repeated the advantages become clearer enabling a 1000 equation system to reduced to less than 50 equations. For the sequence above the simplifed exponential model is used to represent how the reduction of equation is done.

             this sequence {1,2,1,1} would be represented by 4 equations where y- is the degree of the node:

             1 = exp(x1+x2)+exp(x1+x3)+exp(x1+x4)
             2 = exp(x2+x1)+exp(x2+x3)+exp(x2+x4)
             1 = exp(x3+x1)+exp(x3+x2)+exp(x3+x4)
             1 = exp(x4+x1)+exp(x4+x2)+exp(x4+x3)

             however since we know that same degree nodes should have the same lagrange solutions let x1 = x3 = x4 and simplify the equations to the following:

             1 = exp(x1+x2)+2exp(2(x1))
             2 = 3exp(x2+x1)

         Inputs:
             _n: the number of equations to be solves
             _x: the variables in the equation which are initialized outside of the function by x=vector(1, _n);
             _f: the equations which are initialized similar to _x
             Globals:
             these globals are used because the funcv is not allowed to be changed due to the way the numerical recipes library interacts with it (it is defined inside the library as this form)
             SEQ: the degree of the nodes
             QU: the quantity of each node
             SEXP: a boolean determines if the simplifed exponential form is used or the full exponential model is used
         Outputs:
             populates f[] with the representations of the equations to be solved
     */

    int i, j;

    for( i = 1; i <= _n; i++ ) {
        _f[i] = 0;
        for( j = 1; j <= _n; j++ ) {
            if( j == i ) {
                if(SIMPEXP) {
                    _f[i] += ((SEQQUANT[j-1])-1)*exp(2*_x[i]);
										//_f[i] += ((SEQQUANT[j-1]))*exp(2*_x[i]);
                } else {
                    _f[i] += ((SEQQUANT[j-1])-1)/(1+exp(-2*_x[i]));
										//_f[i] += ((SEQQUANT[j-1]))/(1+exp(-2*_x[i]));
                }
            } else {
                if(SIMPEXP) {
                    _f[i] += SEQQUANT[j-1]*exp(_x[i]+_x[j]);
                } else {
                    _f[i] += SEQQUANT[j-1]/(1+exp(-1*(_x[i]+_x[j])));
                }
            }
        }
        _f[i] -= CONVSEQ[i-1];
    }
}

// MARK: Main Function
void attemptConvergence(int _maxIts) {
    /*
         Description: takes an input sequence stored in CONVSEQ whose degree distribution conforms to the equation y = cx^gamma and solves either for the solutions to the non-linear sets of the simplified exponential model whose equations are of the form y=exp(-x1+x2) or the Exponential model whose equations y=1/(1+exp(x1+x2)). To solve these non-linear systems the broyden or the Newton method is used.
         Inputs:
            TO USE THIS FUNCTION:
                First initializeLibrary() passing it a integer which will represent the maximum size of any array used
                Second call setNewSize() passing it the size of the array which you will pass in first. This should be called everytime you send in an array which is different sized than the previous one and always after initializeLibrary() if the size you pass to setNewSize() is larger than that passed to initializedLibrary() it will reset the library allocating enough space to accomodate the new specified size.
                Third set the sequence to be solved for. To set the sequence pass it to setSequence() along with its corresponding quantity array
                Fourth set the guesses for the solutions. To set the guesses for the solution set pass the guesses to setGuesses() this can be set to zero an it should converge to a solution however it will take sometime and may exceed the maximum iterations allowed in the convergence algorithms.
                Fifth set whether you wish to use the FULL or SIMPLIFIED exponetial models using setForSimpExp() and setForFullExp().
                Sixth set for newton or broydn method. Use setForNewton() and setForBroydn() broydn will be faster but may or may not converge as nicely undercertain conditions.
         Outputs:
            XVEC will contain the solutions to obtain these use getSolutions()
            FVEC will contain values which indicate how close the solutions are see getConvValues()
            CONVPROB will contain an indication of whether the convergence has encountered issues see getConditionCode()
     */

    // variable declaration
    int i, check=0;
    ITS = 0;
    struct timeval t1, t2;
    CONVPROB = 0;

    gettimeofday(&t1, NULL);
    // attempt convergence until max its is reached
    for(i = 0; i < _maxIts; i++) {
        // run broydn or newton
        if(NEWTON) {
            newt(XVEC, CURRSIZE, &check, funcv);
            ITS += nrITSN;
        } else {
            broydn(XVEC, CURRSIZE, &check, funcv);
            ITS += nrITSB;
        }
        funcv(CURRSIZE, XVEC, FVEC);
        // set error code
        setCONVPROB(check);
    }
    gettimeofday(&t2, NULL);

    // determine number of miliseconds elapsed
    TIMEELAPSED = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    TIMEELAPSED += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
}

// MARK: Public Setter Functions
int setSequence(int *_seq, int *_seqQa, int _seqSize) {
    /*
        Description: sets CONVSEQ and SEQQUANT with the values in _seq and _seqQa respectively. If the seqSize is greater than the MAXSIZE or is not equal to CURRSIZE then returns 1 and does not set the values
        Inputs:
            _seq: the sequence to be converged upon
            _seqQa: the quantity of each value in _seq. This allows compressed sequences to be entered
            _seqSize: the size of the arrays both should be equal
        Ouputs:
            returns 1 if seqSize is greater than MAXSIZE or not equal to CURRSIZE
            returns 0 if set is sucessful
    */
    int i;
    // check if the library has enough allocated space for the input guesses
    if(_seqSize > MAXSIZE) {
        printf("setSequence: Size of input guesses is too large for current allocation please reset library for more space!\n");
        return 1;
    } else if(_seqSize != CURRSIZE) {
        printf("setSequence: Current size does not match with size of guess array enter a different array\n");
        return 1;
    } else {
        // copy the input array into XVEC
        for( i = 0; i < _seqSize; i++) {
            CONVSEQ[i] = _seq[i];
            SEQQUANT[i] = _seqQa[i];
        }
        return 0;
    }
}
int setGuesses(float *_guesses, int _numGuesses) {
    /*
        Description: Sets the XVEC with the guesses provided to the function. Better guesses will result in quicker and more accurate convergence to the actual solutions.
        Inputs:
            _guesses: an array of floats with the respective guesses for each value in the CONVSEQ array
            _numGuesses: the size of the input array
        Outputs:
            returns 0 if no problem encountered
            returns 1 if the size of the input array does not correspond to the current size of the library or if the size of input array is too large for the current allocation of the library
    */
    int i;
    // check if the library has enough allocated space for the input guesses
    if(_numGuesses > MAXSIZE) {
        printf("setGuesses: Size of input guesses is too large for current allocation please reset library for more space!\n");
        return 1;
    } else if(_numGuesses != CURRSIZE) {
        printf("setGuesses: Current size does not match with size of guess array enter a different array\n");
        return 1;
    } else {
        // copy the input array into XVEC
        for( i = 0; i < _numGuesses; i++) {
            XVEC[i+1] = _guesses[i];
            //printf("\t\ti%f\n", _guesses[i]);

        }
        return 0;
    }
}
void setNewSize(int _size) {
    /*
        Description: Sets the new size of the library to the inputs size. If the input size is greater than max size then reallocates the library to fit the size. Should be set if any new sequence of different size to the previous sequence in the library.
        Input:
            _size: the new size of the library:
    */

    // check if the size is less than MAXSIZE
    if(_size > MAXSIZE) {
        reallocLibrary(_size);
    } else {
        CURRSIZE = _size;
        resetFVEC();
        resetXVEC();
        resetSEQSEQQUANT();
    }
    return;
}
void setForNewton() {
    NEWTON = 1;
}
void setForBroydn() {
    NEWTON = 0;
}
void setForSimpExp() {
    SIMPEXP = 1;
}
void setForFullExp() {
    SIMPEXP = 0;
}

// MARK Internal Setter Functions
void resetSEQSEQQUANT() {
    /*
        Description: resets all CONVSEQ and SEQQUANT values to zero
    */
    int i;
    for( i = 0; i < MAXSIZE; i++) {
        CONVSEQ[i] = 0;
        SEQQUANT[i] = 0;
    }
    return;
}
void resetFVEC() {
    /*
     Description: resets all FVEC values to zero
     */
    int i;
    for( i = 1; i < MAXSIZE+1; i++) {
        FVEC[i] = 0;
    }
    return;
}
void resetXVEC() {
    /*
     Description: resets all XVEC values to zero
     */
    int i;
    for( i = 1; i < MAXSIZE+1; i++) {
        XVEC[i] = 0;
    }
    return;
}
void setFVEC(float *_f, int _fsize) {
    /*
        Description: Sets FVEC to the value sin _f. _f must have been defined using vector function.
    */
    int i;
    for( i = 1; i < _fsize+1; i++) {
        FVEC[i] = _f[i];
    }
    return;
}
void setCONVPROB(int _errorcode) {
  /*
    if(_errorcode == 0) {
        CONVPROB = 0;
    } else if(NEWTON && _errorcode) {
        CONVPROB = 1;
    } else if(!NEWTON && _errorcode) {
        CONVPROB = 2;
    } else {
        CONVPROB = 3;
    }
  */
  //CONVPROB *= 10;
  if(_errorcode == 0) {
    CONVPROB = 0;
  } else if(NEWTON && _errorcode) {
    CONVPROB = 1;
  } else if(!NEWTON && _errorcode == 1) {
    CONVPROB = 2;
  } else if(!NEWTON && _errorcode == 3) {
    CONVPROB = 1;
  } else if (_errorcode == 5) {
    CONVPROB = 4;
  } else {
    CONVPROB = 3;
  }
}

// MARK: Getter Functions
int getConvergenceType() {
    /*
        Description: Returns the value of NEWTON
        Outputs:
            returns 1 if Newtons method is to be used
            returns 0 if broydns method is to be used
    */

    return NEWTON;
}
int getModelType() {
    /*
        Description: Returns the value of SIMEXP
        Output:
            returns 1 if simplified exponential model is being used
            returns 0 if full exponential model is being used
    */

    return SIMPEXP;
}
int getConditionCode() {
    /*
        Description: returns the value of CONVPROB
        Outputs:
            each 10s place in CONVPROB represents the convergence code that occurred at
            that specific iteration. If there are 4 iterations then the condition code should look
            similar to: 1004. This indicates that on the first iteration an issue with
            lnsrch occurred or TOLMIN was exceeded in newton, the second and third iterations
            had no issues and the 4th iteration had a ludcmp error.
            returns 0 if no issues encountered
            returns 1 if an issue with lnsrch occurred or if TOLMIN in newt was exceeded
            returns 2 if a singluar jacobian was encountered in broydn
            returns 3 if MAXITS was exceeded in newt or broydn
            returns 4 if ludcmp error
    */

    return CONVPROB;
}
int getSetSize() {
    /*
        Description: returns CURRSIZE
    */
    return CURRSIZE;
}
int getMaxSize() {
    /*
        Description: returns MAXSIZE
    */
    return MAXSIZE;
}
int getSequence(int *_seq, int *_seqQa) {
    /*
        Description: places the sequence and sequence quant into the arrays given. These arrays should be allocated for enough int elements as CURRSIZE defines
        Inputs:
            _seq: where CONVSEQ will be copied
            _seqQa: where SEQQUANT will be copied
        Outputs:
            returns CURRSIZE
            copies CONVSEQ and SEQQUANT into the respective arrays
    */
    int i;
    for( i = 0; i < CURRSIZE; i++) {
        _seq[i] = CONVSEQ[i];
        _seqQa[i] = SEQQUANT[i];
    }

    return CURRSIZE;
}
int getSolutions(float *_solu) {
    /*
        Description: places XVEC into the array given. This array should be allocated for enough int elements as CURRSIZE defines
        Inputs:
            _solu: where XVEC will be copied
        Outputs:
            returns CURRSIZE
            copies XVEC into the array
     */
    int i;
    for( i = 0; i < CURRSIZE; i++) {
        _solu[i] = XVEC[i+1];
    }

    return CURRSIZE;
}
int getConvValues(float *_convVal) {
    /*
        Description: places FVEC into the array given. This array should be allocated for enough int elements as CURRSIZE defines. FVEC tells how close to the actual solution you have arrived. If all FVEC values are zero then this is the solution, if the values differ greatly from zero then you have not converged upon the solutions.
        Inputs:
            _convVal: where FVEC will be copied
        Outputs:
            returns CURRSIZE
            copies FVEC into the array
     */
    int i;
    for( i = 0; i < CURRSIZE; i++) {
        _convVal[i] = FVEC[i+1];
    }

    return CURRSIZE;
}
void getAccOfSol(float *_minInacc, float *_maxInacc, float *_avgInacc, float *_avgInaccAbs) {
    /*
         Description: Will find the minimum and maximum divergence from the solutions as well as the average divergence from the solutions
             the input array should be the lagrange solutions from the broydn or newton methods.
         Inputs:
             _minInacc: where the lowest value of _sol array is stored
             _maxInacc: where the largest value of the _sol array is stored
             _avgInacc: where the average of all the values in the array is stored
             _avgInaccAbs: where the average value of the absolute value of each element in the _sol array is stored
         Outputs:
             the minimum value in the array is stored in minInacc, the max value of the array is stored in maxInacc, the avg value of the array is stored in avgInacc and the avg of the abs values of the array is stored in avgInaccAbs
     */

    // declare values
    float val, min=-1, max=-1, sum=0, sumAbs=0;
    int i, mi=1, ma=1;

    // loop through array
    for(i = 0; i < CURRSIZE; i++) {
        /* get val */
        val = FVEC[i+1];
        /* add element to sum */
        sum += val;
        /* add abs val of element to sum */
        sumAbs += fabs(val);
        /* determine if the value is min or max */
        if(mi || val < min) {
            min = val;
            mi = 0;
        } else if(ma || val > max) {
            max = val;
            ma = 0;
        } else {

        }
    }

    // get avg
    *_avgInacc = sum/(float)CURRSIZE;
    *_avgInaccAbs = sumAbs/(float)CURRSIZE;

    // assign min and max
    *_minInacc= min;
    *_maxInacc = max;

    return;
}
int getITS() {
  return ITS;
}
double getTimeElapsed() {
  // returns number of miliseconds elapsed during convergence
  return TIMEELAPSED;
}

// MARK: Printing Functions
void printFVECXVEC() {
  printf("FUNCV\n");
  printf("XVEC\tFVEC\n");
  int j;
  for(j = 1; j <= getSetSize();j++){
    printf("%f\t%f\n", XVEC[j], FVEC[j]);
  }
}

// MARK: Sequence Compression/Decompression
int compressSequence(int *_seq, int _size, int *_seq_comp, int *_seq_quant) {
    /*
     Description: compressesses the input _seq into only its unique elements which are placed in _seq_comp and then placing the corresponding quantity of each in a seperate _seq_quant

     example: the following sequence would be compressed as such
     input: _seq {1,1,1,2,3,5,5}
     output: _seq_comp {5,3,2,1}, _seq_quant {2,1,1,3}
     Inputs:
     _seq: the sequence to be compressed
     _size: the size of the input sequence
     _seq_comp: a pointer to the compressed sequence, must be allocated for a least as much space as _seq needs
     _seq_quant: a pointer to the quantity of each element in the compressed sequence must be allocated for as much space as _seq needs
     Outputs:
     the compressed sequence and corresponding quantities of each are placed in _seq_comp and _seq_quant
     the function also returns the size of the _seq_comp and _seq_quant; basically the new size of the compressed sequence
     */

    /* declare variables */
    int num, temp, index, i;

    /* define variables */
    num = 0;
    index = 0;

    for(i = 0; i < _size; i++) {
      _seq_quant[i] = 0;
    }
    /* compress sequence */
    for( i = 0; i < _size; i++) {
        temp = _seq[i];
        if(arrayContainsInternal(temp, _seq_comp, num, &index)) {
            _seq_quant[index] += 1;
        } else {
            _seq_comp[num] = temp;
            _seq_quant[num] += 1;
            num++;
        }
    }

    /* return the size of the compressed sequence */
    return num;
}
int decompressSequence(int *_seq_comp, int *_seq_quant, int _compLen, int *_decomp) {
    /*
     Description: Will decompress a sequence provided with its corresponding quantity array. To be used in conjunction with compressSequence function will undo the effects of this function.

     Example:
     Input:
     _seq_comp = {1, 3, 4}
     _seq_quant = {3, 1, 2}

     Output:
     _decomp = {1,1,1,3,4,4}
     will return 6 to indicate the length of the decompressed sequence
     Inputs:
     _seq_comp: The compressed sequence containing only the unique numbers of the original sequence
     _seq_quant: The quanitity of the corresponding unique number in _seq_comp refers to the quantity of this number which were found in the original uncompressed sequence and therefore the quantity of this number which will be found in the decompressed version of the sequence returned by this function.
     _compLen: The length of the compressed sequence.
     _decomp: The array to which the decompressed sequence will be stored. The memory necessary to store this will be calculated based on the inputs to _seq_quant and then allocated to _decomp
     Outputs:
     The decompressed sequence will be stored in the _decomp array. The resulting size of this decompressed sequence will be returned by the function as an integer.
     */

    /* declare variables */
    int i, j, num, currind, quant;

    /* define variables */
    quant = 0;
    currind = 0;
    num = 0;

    /* determine the size decomp will need to be */
    for(i = 0; i < _compLen; i++) {
        num += _seq_quant[i];
    }

    /* allocate memory to _decomp necessary to hold the decompressed sequence */
    _decomp = malloc(sizeof(int)*num);

    /* expand the compressed sequence into decomp */
    for( i = 0; i < num; i++) {
        quant += _seq_quant[i];
        for(j = currind; j < quant; j++) {
            _decomp[j] = _seq_comp[i];
        }
        currind = quant;
    }

    /* return the size of the decompressed sequence */
    return num;
}
int arrayContainsInternal(int _num, int *_seq, int _size, int *_index) {
    /*
     Description: Determines if the element _num is contained in the array _seq and returns the index of the first occurance of this element
     Inputs:
     _num: the element for which to search
     _seq: the array in which to search
     _size: the size of the array in which to search
     _index: is populated with the index at which the element is found
     Outputs:
     returns 0 if the element is not found and 1 if it is found.
     index becomes the index at which the element is found
     */

    /* variable declaration */
    int output, i;

    /* variable initialization */
    output = 0;
    i = 0;

    /* find index of number if exists in sequence */
    while(!output && i < _size) {
        if(_seq[i] == _num) {
            *_index = i;
            output = 1;
        }
        i++;
    }

    return output;
}
int compressSequenceF(float *_seq, int _size, float *_seq_comp, int *_seq_quant) {
    /*
     Description: compressesses the input _seq into only its unique elements which are placed in _seq_comp and then placing the corresponding quantity of each in a seperate _seq_quant

     example: the following sequence would be compressed as such
     input: _seq {1,1,1,2,3,5,5}
     output: _seq_comp {5,3,2,1}, _seq_quant {2,1,1,3}
     Inputs:
     _seq: the sequence to be compressed
     _size: the size of the input sequence
     _seq_comp: a pointer to the compressed sequence, must be allocated for a least as much space as _seq needs
     _seq_quant: a pointer to the quantity of each element in the compressed sequence must be allocated for as much space as _seq needs
     Outputs:
     the compressed sequence and corresponding quantities of each are placed in _seq_comp and _seq_quant
     the function also returns the size of the _seq_comp and _seq_quant; basically the new size of the compressed sequence
     */

    /* declare variables */
    int index, i, num;
    float temp;
    /* define variables */
    num = 0;
    index = 0;

    /* compress sequence */
    for( i = 0; i < _size; i++) {
        temp = _seq[i];
        if(arrayContainsInternalF(temp, _seq_comp, num, &index)) {
            _seq_quant[index] += 1;
        } else {
            _seq_comp[num] = temp;
            _seq_quant[num] += 1;
            num++;
        }
    }

    /* return the size of the compressed sequence */
    return num;
}
int arrayContainsInternalF(float _num, float *_seq, int _size, int *_index) {
    /*
     Description: Determines if the element _num is contained in the array _seq and returns the index of the first occurance of this element
     Inputs:
     _num: the element for which to search
     _seq: the array in which to search
     _size: the size of the array in which to search
     _index: is populated with the index at which the element is found
     Outputs:
     returns 0 if the element is not found and 1 if it is found.
     index becomes the index at which the element is found
     */

    /* variable declaration */
    int output, i;

    /* variable initialization */
    output = 0;
    i = 0;

    /* find index of number if exists in sequence */
    while(!output && i < _size) {
        if(_seq[i] == _num) {
            *_index = i;
            output = 1;
        }
        i++;
    }

    return output;
}

// MARK: Lagrange Estimation/Reverse Lagrange to Degree sequence
void CL_LagrangeEstimation(int _seqSize, int *_seq, float *_expLG, int *_quantLG, float *_Cvals) {
  /* uses Dr. Bassler's proofs for the ChungLu model to estimate the Lagrange
     multipliers for a degree sequence */
     int j, i, *temp;
     float C=0, *allLG, sumDegree=0.0;

     allLG = (float*) malloc(sizeof(float) * _seqSize);
     temp = (int*) malloc(sizeof(int) * _seqSize);

     for(j = 0; j < _seqSize; j++) {
       sumDegree += _seq[j]*_quantLG[j];
     }

     //C = pow(threeFact*(1.0/3.0)*sumDegree, 2.0/3.0);
     sumDegree = pow(sumDegree, 1.0/2.0);
     //printf("%f\n", sumDegree);
     for(i = 0; i < _seqSize; i++) {
       //printf("%f\t%d\t%d\n", log((float)_seq[i]/C)/-2.0, _seq[i], _quantLG[i]);
       //_expLG[i] = log((float)_seq[i]/C)/-2.0;
       _expLG[i] = log((float)_seq[i]/sumDegree);
     }

     /* compress sequence */
     for(i = 0; i < _seqSize; i++) {
       _Cvals[i] = (float)_seq[i]/C;
     }

     return;
}
void LagrangeToDegreeSeq(int _seqSize, int *_expseq, float *_LG, int *_quant, int *_resultSeq) {
  /* converts the lagrange solutions to a degree sequence, if the lagrange values
     are correct then they should return the same sequence as was used to find them */

     int j, i, sumDegree=0, threeFact = 6;
     for(j = 0; j < _seqSize; j++) {
       sumDegree += _expseq[j]*_quant[j];
       //printf("%d\t%d\n", _expseq[j], _quant[j]);
     }

     //float C = pow(threeFact*(1.0/3.0)*sumDegree, 2.0/3.0);
     printf("\tLagrangeToDegreeSeq\t%d\t%d\n", sumDegree, _seqSize);

     sumDegree = pow(sumDegree, .5);
     for(i = 0; i < _seqSize; i++) {
       //printf("%f\n", exp(-2.0*_LG[i])*C);
       //_resultSeq[i] = exp(-2.0*_LG[i]) * C + .5;
       _resultSeq[i] = (sumDegree * exp(_LG[i])) + .5;
       printf("%d\t%d\t%f\n", _expseq[i], _resultSeq[i], _LG[i]);
     }

}
int CheckLagrangeAccuracy(int _seqSize, int *_seq, int *_quant, float *_LG, float _tol, float *_acc) {
  /* uses funcv to determine how accurate the lagrange solutions provided are */

  int *tempSEQQ, *tempCONVSEQ, error = 0;
  float *tempacc, *tmpLG;
  int i, j, tempsize;
  tempSEQQ = (int*) malloc(getSetSize() * sizeof(int));
  tempCONVSEQ = (int*) malloc(getSetSize() *sizeof(int));
  tempacc = (float*) malloc((_seqSize+1) * sizeof(float));
  tmpLG = (float*) malloc((_seqSize+1) * sizeof(float));

	*_acc = 0;

  for(i = 0; i < getSetSize(); i++) {
    tempSEQQ[i] = SEQQUANT[i];
    tempCONVSEQ[i] = CONVSEQ[i];
  }

  for(i = 0; i < _seqSize; i++ ) {
    tmpLG[i+1] = _LG[i];
  }

  tempsize = getSetSize();

  setNewSize(_seqSize);
  setSequence(_seq, _quant, _seqSize);

  funcv(_seqSize, tmpLG, tempacc);

  //printf("\tFUNCVACC\n");
  for(i = 1; i <= _seqSize; i++) {
    if(fabsf(tempacc[i]) > _tol) error = 1;
    *_acc += fabsf(tempacc[i]);

    /*
    if(_acc[i] != _acc[i]) {
      // print the sequence and the corresponding quant values
      printf("Seqvals\tseqquant\tseqsize=%d\n", _seqSize);
      for(j = 0; j < _seqSize; j++) {
        printf("%f\t%f\t%f\n", _seq[j], _quant[j], _LG[i]);
      }
      // break from loop
      error = 1;
      break;
    }
    */
    //printf("\t%d\t%f\t%f\n", _seq[i-1], _LG[i-1], _acc[i]);
  }
	*_acc /= _seqSize;
  setNewSize(tempsize);
  setSequence(tempCONVSEQ, tempSEQQ, tempsize);
  return error;
}
void estimateBetaValues(int _seqSize, int *_seqComp, int *_seqQnt, float *_guess, float _gamma, long *_seed) {
	/*
	*/

	unsigned int i, seqSum=0;
	double denomCL;

	if(_gamma < -3 || 1 == 1) {
		/* provide ChungLu estimated Beta Values */
		/* find the sequence summation */
			for(i = 0; i < _seqSize; i++) {
				seqSum += _seqComp[i] * _seqQnt[i];
			}

		/* calcuate the denominator for the calculation */
			denomCL = pow(seqSum, 1.0/2.0);

		/* determine the beta values for the unique sequence values */
			for(i = 0; i < _seqSize; i++) {
				_guess[i] = log((double)_seqComp[i]/denomCL);
			}

	} else {
		/* randomly generate coefficient values within the range specified for _gamma and _seqSize */
	}
}
