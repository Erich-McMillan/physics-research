//
//  InvTransDisr.h
//  ExponentialGraphModel
//
//  Created by erich on 6/15/16.
//  Copyright (c) 2016 erich. All rights reserved.
//

// MARK: Definitions
#ifndef __ExponentialGraphModel__InvTransDisr__
#define __ExponentialGraphModel__InvTransDisr__

// MARK: Includes
#include <stdio.h>

/* MARK: External Defintions */
extern float GAMMA, LOGSUM;
extern int SEQSIZE, I;

// MARK: Function Defintions
void InvTransDisr(long* seedpointer, int* sequence);
void setLOGSUM();
void setSEQSIZE(int _seqSize);
void setGAMMA(float _gamma);
int ErdosGallai(int *seq, int size);
int numEdgesEven(int size, int *seq);

// MARK: End
#endif /* defined(__ExponentialGraphModel__InvTrans__) */
