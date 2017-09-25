//
//  SequenceFunctions.c
//  ExponentialGraphModel
//
//  Created by erich on 4/8/16.
//  Copyright © 2016 erich. All rights reserved.
//

#include "SequenceFunctions.h"

#include <stdlib.h>
#include <stdio.h>

// MARK: Array Sort/Find Functions
int arrayContains(int _num, int *_seq, int _size, int *_index) {
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
int quantityInArray(float _num, float *_arr, int _size) {
    /*
     Description: the function will find the quantity of the specified element contained in the array.
     Inputs:
     _num: the element for which to search
     _arr: the array in which to search
     _size: the size of the array to search
     Outputs:
     returns the quanitity of the element _num found in the array. If it doesn't exist 0 will be returned
     */
    
    int numArray, i;
    
    numArray = 0;
    
    for(i = 0; i < _size; i++) {
        if(_arr[i] == _num) {
            numArray++;
        }
    }
    
    return numArray;
}
void minimumValue(int _size, float *_seq, int *_index, float *_val) {
    /*
        Description: finds the minimum value in the array and returns the min value in _val and the index at which it occurs in _index
    */
    
    float currMin = _seq[0];
    int index = 0;
    
    for(int i = 1; i < _size; i++) {
        if(currMin > _seq[i]) {
            currMin = _seq[i];
            index = i;
        }
    }
    
    *_val = currMin;
    *_index = index;
    
    return;
}

