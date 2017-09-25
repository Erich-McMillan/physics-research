//
//  InvTrans.h
//  ExponentialGraphModel
//
//  Created by erich on 10/14/15.
//  Copyright (c) 2015 erich. All rights reserved.
//

// MARK: Definitions
#ifndef __ExponentialGraphModel__InvTrans__
#define __ExponentialGraphModel__InvTrans__

// MARK: Includes
#include <stdio.h>

// MARK: Function Defintions
float InvTrans(int seqsize, double gamma, long* seedpointer, int* sequence);
float InvTransMaxIts(int seqsize, double gamma, long* seedpointer, int* sequence, unsigned int maxIts);
float InvTransNoCheck(int seqsize, double gamma, long* seedpointer, int* sequence);
int ErdosGallai(int *seq, int size);

// MARK: End
#endif /* defined(__ExponentialGraphModel__InvTrans__) */
