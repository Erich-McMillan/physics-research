//
//  HelperFunctions.c
//  ExponentialGraphModel
//
//  Created by erich on 3/28/16.
//  Copyright © 2016 erich. All rights reserved.
//

// MARK: Includes
#include "FileFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// MARK: Internal Function Declarations
void getFilepath(char* filepath);

// MARK: Functions
FILE* openFile(char *_filePath, char *_openFormat) {
    /*
        Description: Opens the file in the format specified by openFormat at the location filepath and points fp to the file.
        Inputs:
            _filepath: the location of the file
            _fp: the pointer which is populated with the file
            _openFormat: the format in which the file is opened; such as "rb"
        Outputs:
            returns the file pointer if the file is opened successfully NULL if a problem occurred
    */
    
    // declare variables
    FILE *fp;
    
    // Open the file. If opening fails return null
    if((fp = fopen(_filePath, _openFormat)) == NULL) {
        return NULL;
    } else {
        return fp;
    }
}
FILE* openFilePrompt(char *_fileName, char *_openFormat) {
    /*
        Description: will prompt the user for a destination to store the file named _fileName and then attempt to open the file there. If this fails it will notify the user and reprompt for a different destination.
        Inputs:
            _fileName: the name of the file
            _openFormat: the way the file is opened. "w", "wb", "r" ect.
        Outputs:
            returns a filepointer to the file once an appropriate file name is given
    */
    
    // declare variables
    char filepath[500];
    FILE *fp;
    
    do {
        for(int i = 0; i < 500; i++) {
            filepath[i] = '\0';
        }
        getFilepath(filepath);
        strcat(filepath, _fileName);
        fp = openFile(filepath, _openFormat);
        if(fp == NULL) printf("Error opening file try a different filepath:\n");
    } while(fp == NULL);
    
    return fp;
}

// MARK: Internal Functions
void getFilepath(char *_filepath) {
    /*
        Description: will prompt user for a filepath and place it in the input char array
        Inputs:
            _filepath: this is where the input filepath will be stored should be allocated for at least 100 characters
        Outputs:
            will place the filename into _filepath
     
    */
    
    // declare variables
    char t;
    int index = 0;
    
    // prompt user
    printf("Enter the filepath excluding the filename at which the program results will be stored:\n");
    
    // grab filepath value
    scanf("%c", &t);
    while(t != '\n') {
        _filepath[index] = t;
        index++;
        scanf("%c", &t);
    }
    
    // make sure the string is terminated with '/'
    if(_filepath[strlen(_filepath)-1]) {
        _filepath[strlen(_filepath)] = '/';
    }
        
    return;
}
