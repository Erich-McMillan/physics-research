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
extern int *CONVSEQ, *SEQQUANT, MAXSIZE, CURRSIZE, SIMPEXP, NEWTON, CONVPROB, ITS;
extern float *FVEC, *XVEC;
extern double TIMEELAPSED;

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
int getITS();
double getTimeElapsed();

// MARK: Sequence Compression/Decompression Functions
int compressSequence(int* sequence, int size, int* compressedseq, int* seqquant);
int decompressSequence(int* compressedseq, int* seqquant, int compressedlength, int* decompressedseq);
int compressSequenceF(float* _seq, int _size, float* _seq_comp, int* _seq_quant);

// MARK: Printing Functions
void printFVECXVEC();

void CL_LagrangeEstimation(int _seqSize, int *_seq, float *_expLG, int *_quantLG, float *_Cvals);
void LagrangeToDegreeSeq(int _seqSize, int *_expseq, float *_LG, int *_quant, int *_resultSeq);
int CheckLagrangeAccuracy(int _seqSize, int *_seq, int *_quant, float *_LG, float _tol, float *_acc);
void estimateBetaValues(int _seqSize, int *_seqComp, int *_seqQant, float *_guess, float _gamma, long *_seed);
// MARK: End
#endif /* ConvergenceLibrary_h */
