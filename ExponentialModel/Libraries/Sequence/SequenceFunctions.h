//
//  SequenceFunctions.h
//  ExponentialGraphModel
//
//  Created by erich on 4/8/16.
//  Copyright © 2016 erich. All rights reserved.
//

#ifndef SequenceFunctions_h
#define SequenceFunctions_h

// MARK: Includes
#include <stdio.h>

// MARK: Array Sort/Find Functions
int arrayContains(int element, int* sequence, int seqsize, int* indexofelement);
int quantityInArray(float element, float* sequence , int seqsize);
void minimumValue(int size, float* seq, int* index, float* value);
// MARK: End
#endif /* SequenceFunctions_h */
