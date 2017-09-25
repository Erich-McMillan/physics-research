//
//  ConvergenceLibrary.h
//  ExponentialGraphModel
//
//  Created by erich on 4/8/16.
//  Copyright © 2016 erich. All rights reserved.
//

#ifndef ConvergenceLibrary_h
#define ConvergenceLibrary_h

// MARK: Includes
#include <stdio.h>

// MARK: Globals
extern int *CONVSEQ, *SEQQUANT, MAXSIZE, CURRSIZE, SIMPEXP, NEWTON, CONVPROB;
extern float *FVEC, *XVEC;

// MARK: Initialize Functions
void initializeLibrary(int maxseqsize);
void freeLibrary();
void reallocLibrary(int newmaxseqsize);
void resetLibrary(int newmaxseqsize);

// MARK: Main Function
void attemptConvergence(int maximumiterations);

// MARK: Public Setter Functions
int setSequence(int* sequence, int* sequencequant, int seqsize);
int setGuesses(float* guesses, int numguesses);
void setNewSize(int newsize);
void setForNewton();
void setForBroydn();
void setForSimpExp();
void setForFullExp();

// MARK: Public Getter Functions
int getConvergenceType();
int getModelType();
int getConditionCode();
int getSetSize();
int getMaxSize();
int getSequence(int* seqstorage, int* seqquantstorage);
int getSolutions(float* solstorage);
int getConvValues(float* convstorage);
void getAccOfSol(float* mininacc, float* maxinacc, float* avginacc, float* avginaccabs);

// MARK: Sequence Compression/Decompression Functions
int compressSequence(int* sequence, int size, int* compressedseq, int* seqquant);
int decompressSequence(int* compressedseq, int* seqquant, int compressedlength, int* decompressedseq);

// MARK: End
#endif /* ConvergenceLibrary_h */
