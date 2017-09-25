//
//  solveLagr.c
//  ExponentialGraphModel
//
//  Created by erich on 6/28/2017.
//  Copyright © 2017 erich. All rights reserved.
//

/* MARK: Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>

#include "ConvergenceLibrary.h"

/* MARK: Function Declarations */
int ImprovedConvergence(int *_seqCompr, int *_seqQnty, int _seqcomplgth, float _gamma, float _tol, long *_seed, float *_convSol, float *convFVEC, float *_convacc, float *_estimated);

/* MARK: Main */
int main( int argc, char **argv ) {
	/*
			Description:
				This program will use broyden and newton methods to find the lagrangian multipliers for the ERGM
				model based on an input sequence which conforms to the degree distribution y=Cx^-gamma

			Inputs:
				argv[1]: (int) _sizeNtw; the size of the network being evaluated
				argv[2]: (int) _seqNum; the sequence identifying number
				argv[3]: (float) _tol; the min accuracy tolerance of the convergence program
				argv[4]: (float) _gamma; the exponent of the network being evaluated
				argv[5]: (string) _pathSeq; the full path to the network being evaluated
				argv[6]: (string) _pathLagr; the full path to where the results shall be stored

			Outputs:
				The program will write the degree, resulting lagrangian multipliers and the convergence accuracy
				to a file in the specified location, each set of values to its own column and separated by \t
	*/

	/* Display help statement if no inputs and exit */
	if(argc < 2) {
		printf("\n\nDescription:\nThis program will use broyden and newton methods to find the lagrangian multipliers for the ERGM\nmodel based on an input sequence which conforms to the degree distribution y=Cx^-gamma\n\nInputs:\nargv[1]: (int) _sizeNtw; the size of the network being evaluated\nargv[2]: (int) _seqNum; the sequence identifying number\nargv[3]: (float) _tol; the min accuracy tolerance of the convergence program\nargv[4]: (float) _gamma; the exponent of the network being evaluated\nargv[5]: (string) _pathSeq; the full path to the network being evaluated\nargv[6]: (string) _pathLagr; the full path to where the results shall be stored\nOutputs:\nThe program will write the degree, resulting lagrangian multipliers and the convergence accuracy\nto a file in the specified location, each set of values to its own column and separated by \\t\n\n\n");

		return 1;
	}

	/* Declare argument variables */
	int _sizeNtw, _seqNum;
	float _tol, _gamma;
	long _seed = -123;
	char _pathSeq[200], _pathLagr[200];

	/* get argument variables */
	_sizeNtw = atoi(argv[1]);
	_seqNum = atoi(argv[2]);
	_tol = atof(argv[3]);
	_gamma = atof(argv[4]);
	strcpy(_pathSeq, argv[5]);
	strcpy(_pathLagr, argv[6]);

	/* declare program variables */
	int i, j, k, s, bytadd, tmpval, ensedgcnt, ensGCC, ensTRI, ensSQR, seqcomplgth;
	int *seqIn, *ensGenSq, *ensLL3, *seqQnty, *seqCompr;
	int convcode; // the returned value from ImprovedConvergence() if 0 then no convergence error
	char nameOutLagr[200]={};
	float convacc;
	float *convSols, *convFVEC, *estimated;
	FILE *inSeq, *outLagr;

	/* allocate memory */
	seqIn      = (int*)     malloc(  sizeof(int) * _sizeNtw);
	seqQnty    = (int*)     malloc(  sizeof(int) * _sizeNtw);
	seqCompr   = (int*)     malloc(  sizeof(int) * _sizeNtw);
	convSols   = (float*)   malloc(  sizeof(int) * _sizeNtw);
	convFVEC   = (float*)   malloc(  sizeof(int) * _sizeNtw);
	estimated  = (float*) malloc(sizeof(float) * _sizeNtw);


	/* Open files */
	bytadd = sprintf(nameOutLagr, _pathLagr);
	sprintf(nameOutLagr+bytadd, "lagrMultiN%dr%.1fSq%d.txt",_sizeNtw, fabsf(_gamma), _seqNum);
	if((inSeq = fopen(_pathSeq, "r")) == NULL) {
		printf("Error opening seqFile exiting to system:\n");
		return 2;
	}
	if((outLagr = fopen(nameOutLagr, "w")) == NULL) {
		printf("Error opening file %s\n", nameOutLagr);
		return 2;
	}

	/* set up ConvergenceLibrary */
	initializeLibrary(_sizeNtw);
	setForFullExp();

	/* Generate/Retrieve Sequence */
	rewind(inSeq);
	/* load sequence from file */
	for(s = 0; s < _sizeNtw; s++) {
		fscanf (inSeq, "%d", &tmpval);
		seqIn[s] = tmpval;
	}

	/* Compress Sequence */
	seqcomplgth = compressSequence(seqIn, _sizeNtw, seqCompr, seqQnty);
	setNewSize(seqcomplgth);
	setSequence(seqCompr, seqQnty, seqcomplgth);

	/* Attempt Convergence */
	convcode = ImprovedConvergence(seqCompr, seqQnty, seqcomplgth, _gamma, _tol, &_seed, convSols, convFVEC, &convacc, estimated);

	printf("\n\n_______________\n\n");
	printf("\tConvVal\t\tAccConv\n");
	for(i = 0; i < seqcomplgth; i++) {
		printf("\t%f\t%f\n", convSols[i], convFVEC[i]);
	}

	/* write results to file */
	for(i = 0; i < seqcomplgth; i++) {
		fprintf(outLagr, "%d\t%.15f\t%d\t%.15f\t%.15f\n", seqCompr[i], convSols[i], seqQnty[i], convFVEC[i], estimated[i] );
	}

	/* Close files */
	fclose(inSeq);
	fclose(outLagr);

	/* Clear arrays */
	free(seqIn);
	free(seqQnty);
	free(seqCompr);
	free(convSols);
	free(convFVEC);

	/* Exit */
	return 0;
}

/* MARK: ImprovedConvergence */
int ImprovedConvergence(int *_seqCompr, int *_seqQnty, int _seqcomplgth, float _gamma, float _tol, long *_seed, float *_convSol, float *_convFVEC, float *_convacc, float *_estimated) {
	/*
			Description:
				This function will attempt to converge upon the solution to the non-linear system of equations described by the exponential random graph model. To do so it will use the guesses for the model provided by the function 'estimateBetaValues()' using these guesses convergence attempts will be made using both the newton's an broyden's methods for solving systems of equations. Due to the difficulty of converging to solutions for power-law sequences between gamma = -3 and -2 the method will attempt to converge multiple times (this is due to the fact that 'estimateBetaValues()' will, when gamma is between -3 and -2, provide randomly changing guesses to the solutions based around the average known curve fit of the beta values for the given sequence size and gamma value. )

			Inputs:
				_seqCompr: an array containing the unique values in the power-law sequence
				_seqQnty: the corresponding quantity of each unique value in _seqCompr
				_seqcomplgth: the length of _seqCompr and _seqQnty
				   _gamma: the exponent value of the power-law distributed sequences
					   _tol: the maximum permitted error for any individual equation in the system. It is expected that for a complete convergence every equation will evaluate to 0, however due to singular jacobian matricies it is more convenient to converge to a tolerance. The closer the tolerance to 0 the more accurate the solutions however a convergence to small _tol values is not guarenteed.
						_seed: the seed necessary for the random number generator
				 _convSol: the array into which the solutions to the beta values will be stored
				_convFVEC: the solutions to the system of equations the closer each value to zero the better the convergence
				 _convacc: the average of _convFVEC if small and the integer returned by this function has a ones digit of 0 then the convergence worked with the _tol specified

			Outputs:
				returns a two digit integer. The 10's digit indicates the error code encountered by the systems-solver the 1's digit indicates if _tol was exceeded by the found solution. If the value is 0 then the convergence worked with the _tol specified to note the error codes see the 'getConditionCode()' function in ConvergenceLibrary.c
	*/

	//printf(" INITIALIZATION ImprovedConvergence\n");
	/* MARK: INITIALIZATION */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */
			/* MARK: Variable declarations */
				float *bestSol, *bestFVEC, *currSolN, *currFVECN, *currSolB, *currFVECB, *currGuess;
				float accBst=DBL_MAX, accB, accN;
				int itsMax = _seqcomplgth*2, itsCur = 0, convIts=100, i, j, prev;
				/* Declare flag variables */
				int errBst, tolBst, errB, tolB, errN, tolN, errOut;
				/* Allocate Memory */
				bestSol   = (float*) malloc(sizeof(float) * _seqcomplgth);
				bestFVEC  = (float*) malloc(sizeof(float) * _seqcomplgth);
				currSolN  = (float*) malloc(sizeof(float) * _seqcomplgth);
				currFVECN = (float*) malloc(sizeof(float) * _seqcomplgth);
				currSolB  = (float*) malloc(sizeof(float) * _seqcomplgth);
				currFVECB = (float*) malloc(sizeof(float) * _seqcomplgth);
				currGuess = (float*) malloc(sizeof(float) * _seqcomplgth);

	/* MARK: MAIN BODY */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */
	//printf(" MAINBODY ImprovedConvergence\n");
	estimateBetaValues(_seqcomplgth, _seqCompr, _seqQnty, _estimated, _gamma, _seed);
			/* MARK: check gamma value */
				if(_gamma < -3.0) {
					/* MARK: Get ChungLu guesses and feed solutions back as guesses */
						estimateBetaValues(_seqcomplgth, _seqCompr, _seqQnty, currGuess, _gamma, _seed);
						setGuesses(currGuess, _seqcomplgth);
						setForNewton();


						for(i = 0; i < 100; i++) {
							//printf(" %d\n", i);
							attemptConvergence(convIts);
							getSolutions(currSolN);
							setGuesses(currSolN, _seqcomplgth);
						}

						getConvValues(currFVECN);
						errN = getConditionCode();
						tolN = CheckLagrangeAccuracy(_seqcomplgth, _seqCompr, _seqQnty, currSolN, _tol, &accN);

						accBst = accN;
						errBst = errN;
						tolBst = tolN;
						for(i = 0; i < _seqcomplgth; i++) {
							bestSol[i] = currSolN[i];
						  bestFVEC[i] = currFVECN[i];
						}
				} else {
					/* MARK: Loop for convergence */
						do {
							/* Estimate Lagragne Multipliers */
								estimateBetaValues(_seqcomplgth, _seqCompr, _seqQnty, currGuess, _gamma, _seed);
								setGuesses(currGuess, _seqcomplgth);

							/* Attempt Convergence for Broyden using currGuess */
								setForBroydn();
								attemptConvergence(convIts);
								getSolutions(currSolB);
								getConvValues(currFVECB);
								errB = getConditionCode();
								tolB = CheckLagrangeAccuracy(_seqcomplgth, _seqCompr, _seqQnty, currSolB, _tol, &accB);

							/* Attempt Convergence for Newton using CurrGuess */
								setForNewton();
								attemptConvergence(convIts);
								getSolutions(currSolN);
								getConvValues(currFVECN);
								errN = getConditionCode();
								tolN = CheckLagrangeAccuracy(_seqcomplgth, _seqCompr, _seqQnty, currSolN, _tol, &accN);

							/* Test for best convergence */
								if(accN > accBst && accB > accBst) {
									/* do nothing */
								} else if(accN < accB) {
									/* newton's converged to better answer than Best & broyden */
									accBst = accN;
									errBst = errN;
									tolBst = tolN;
									for(i = 0; i < _seqcomplgth; i++) {
										bestSol[i] = currSolN[i];
										bestFVEC[i] = currFVECN[i];
									}
								} else {
									/* broyden's converged to better answer than Best & newton */
									accBst = accB;
									errBst = errB;
									tolBst = tolB;
									for(i = 0; i < _seqcomplgth; i++) {
										bestSol[i] = currSolB[i];
										bestFVEC[i] = currFVECB[i];
									}
								}

							/* Increment its counter */
								itsCur ++;

						} while(itsCur < itsMax);
					}

				prev = 0;
				for(i = 0; i < _seqcomplgth; i++) {
						_convSol[i] = bestSol[i];
						_convFVEC[i] = bestFVEC[i];
				}
				*_convacc = accBst;
				errOut = errBst*10;
				errOut += tolBst;
	/* MARK: END */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */
			/* MARK: Free memory */
				free(bestSol);
				free(bestFVEC);
				free(currSolN);
				free(currFVECN);
				free(currSolB);
				free(currFVECB);
				free(currGuess);

			/* MARK: Return */
				return errOut;
}
