/* MARK: Includes */
	#include <stdio.h>
	#include <stdlib.h>
	#include <unistd.h>
	#include <ctype.h>
	#include <float.h>
	#include <string.h>
	#include <math.h>
	#include <sys/types.h>
  #include <sys/wait.h>

	#include "ConvergenceLibrary.h"
	#include "InvTransDisr.h"
	#include "ran2.h"


/* MARK: Function Declarations */
	int ImprovedConvergence(int *_seqCompr, int *_seqQnty, int _seqcomplgth, float _gamma, float _tol, long *_seed, float *_convSol, float *convFVEC, float *_convacc);
	int GenerateGraph(int _N, long *_sd, float **_eP, int **_adt, int *_eSq, int**_aMx, int** _eLt);
	void** malloc2d(size_t numRows, size_t rowSize);
	void free2d(void **a);
	void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_L3_local);
	void L4(int *_zdegseq, int **_zat, int **_zam, int **_zel, int _N, int _edgenum, int *_L4_sum);

/* MARK: main */
	int main (int argc, char **argv) {
		/*
		*/

		/* MARK: INITIALIZATION */
		/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */
				/* MARK: Declare argument variables */
					int argVal, _grphCheck=0, _sizeNtw, _sizeEns, _numEns=1, _seqNum;
					float _tol, _gamma;
					long _seed;
					char _pathSeq[200], _pathStrg[200];

				/* MARK: Parse arguments from command line */
					opterr = 0;
					while ((argVal = getopt (argc, argv, "N:e:g:t:r:p:s:n:")) != -1) {
						switch (argVal) {
							case 'N':
								_sizeNtw = atoi(optarg);
								break;
							case 'e':
								_sizeEns = atoi(optarg);
								break;
							case 'g':
								_gamma = atof(optarg);
								break;
							case 't':
								_tol = atof(optarg);
								break;
							case 'r':
								_seed = atol(optarg);
								break;
							case 'p':
								strcpy(_pathStrg, optarg);
								break;
							case 's':
							 	strcpy(_pathSeq,optarg);
								break;
							case 'n':
								_seqNum = atoi(optarg);
								break;
							case '?':
								if (optopt == 'c')
								fprintf(stderr, "Option -%c requires an argument.\n", optopt);
								else if (isprint(optopt))
								fprintf(stderr, "Unknown option `-%c'.\n", optopt);
								else
								fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
								return 1;
							default:
								abort ();
						}
					}
					printf(" _sizeNtw = %d\n _sizeEns = %d\n _gamma = %f\n _tol = %f\n _seed = %ld\n",
					_sizeNtw, _sizeEns, _gamma, _tol, _seed);
					printf(" _pathSeq = %s\n _pathStrg = %s\n",
					_pathSeq, _pathStrg);

				/* MARK: Set-up libraries and runtime variables */
					/* set up InvTransDisr */
					setSEQSIZE(_sizeNtw);
					setGAMMA(_gamma);
					setLOGSUM();
					/* set up ConvergenceLibrary */
					initializeLibrary(_sizeNtw);
					setForFullExp();

				/* MARK: Declare function variables */
					/* declare numeric variables */
					int i, j, k, s, bytadd, tmpval, ensedgcnt, ensGCC, ensTRI, ensSQR, seqcomplgth;
					int *seqIn, *ensGenSq, *ensLL3, *seqQnty, *seqCompr;
					int **ensAdjTb, **ensAdjMtx, **ensEdgeLst;
					char nameOutDL[200]={}, nameOutTSG[200]={}, nameOutAdj[200]={};
					float convacc;
					float *convSols, *convFVEC;
					float **edgeProb;
					/* Declare flag variables */
					int convcode; // the returned value from ImprovedConvergence() if 0 then no convergence error
					/* Declare filevariables */
					FILE *inSeq, *outDL, *outTSG/*, *outAdj*/;

					/* allocate Memory */
					seqIn      = (int*)     malloc(  sizeof(int) * _sizeNtw);
					ensGenSq   = (int*)     malloc(  sizeof(int) * _sizeNtw);
					ensLL3     = (int*)     malloc(  sizeof(int) * _sizeNtw);
					seqQnty    = (int*)     malloc(  sizeof(int) * _sizeNtw);
					seqCompr   = (int*)     malloc(  sizeof(int) * _sizeNtw);
					ensEdgeLst = (int**)    malloc2d(_sizeNtw * _sizeNtw , sizeof(int) * 2);
					ensAdjTb   = (int**)    malloc2d(_sizeNtw            , sizeof(int) * _sizeNtw);
					ensAdjMtx  = (int**)    malloc2d(_sizeNtw            , sizeof(int) * _sizeNtw);
					convSols   = (float*)   malloc(  sizeof(int) * _sizeNtw);
					convFVEC   = (float*)   malloc(  sizeof(int) * _sizeNtw);
					edgeProb   = (float**)  malloc2d(_sizeNtw            , sizeof(float) * _sizeNtw);

				/* MARK: Open files for writing */
					bytadd = sprintf(nameOutDL, _pathStrg);
					sprintf(nameOutDL+bytadd, "/DegreeLCC_%d_%.4d_%.3f_%d.txt",_sizeNtw, _seqNum, _gamma, _sizeEns);
					bytadd = sprintf(nameOutTSG, _pathStrg);
					sprintf(nameOutTSG+bytadd, "/TriSqrGcc_%d_%.4d_%.3f_%d.txt",_sizeNtw, _seqNum, _gamma, _sizeEns);
					bytadd = sprintf(nameOutAdj, _pathStrg);
					sprintf(nameOutAdj+bytadd, "/AdjMat_%d_%d_%.4f_%d.txt",_sizeNtw, _seqNum, _gamma, _sizeEns);
					if((inSeq = fopen(_pathSeq, "r")) == NULL) {
						printf("Error opening seqFile exiting to system:\n");
						return 2;
					}
					if((outDL = fopen(nameOutDL, "w")) == NULL) {
						printf("Error opening file %s\n", nameOutDL);
						return 2;
					}
					if((outTSG = fopen(nameOutTSG, "w")) == NULL) {
						printf("Error opening file %s\n", nameOutTSG);
						return 2;
					}
					/*
					if((outAdj = fopen(nameOutAdj, "w")) == NULL) {
						printf("Error opening file %s\n", nameOutAdj);
						return 2;
					}
					*/
		printf("\n BEGINNING MAIN\n");
		/* MARK: BEGINNING MAIN */
		/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */

				/* MARK: GET/GENERATE SEQUENCE AND CONVERGENCE UPON SOLUTIONS */
						printf(" RETRIEVING SEQUENCE\n");
						/* Generate/Retrieve Sequence */
						rewind(inSeq);
						/* load sequence from file */
						for(s = 0; s < _sizeNtw; s++) {
							fscanf (inSeq, "%d", &tmpval);
							seqIn[s] = tmpval;
  					}

						printf(" COMPRESSING SEQUENCE\n");
						/* Compress Sequence */
						seqcomplgth = compressSequence(seqIn, _sizeNtw, seqCompr, seqQnty);
						setNewSize(seqcomplgth);
						setSequence(seqCompr, seqQnty, seqcomplgth);

						printf(" ATTEMPTING CONVERGENCE\n");
						/* Attempt Convergence */
						convcode = ImprovedConvergence(seqCompr, seqQnty, seqcomplgth, _gamma, _tol, &_seed, convSols, convFVEC, &convacc);
						/*
						printf("\n\n_______________\n\n");
						printf("\tConvVal\t\tAccConv\n");
						for(i = 0; i < _sizeNtw; i++) {
							printf("\t%f\t%f\n", convSols[i], convFVEC[i]);
						}
						printf("\n");
						*/
				/* MARK: GENERATE ENSEMBLE */
						printf(" CALCULATING EDGE PROBABILITIES\n");
						/* Calculate probablities of edges */
						for(j = 0; j < _sizeNtw; j++) {
							for(k = 0; k < _sizeNtw; k++) {
								edgeProb[j][k] = 1.0/(1.0+exp((double)(-1.0*(convSols[j]+convSols[k]))));
							}
						}

						printf(" ENSEMBLE STARTED\n");
						/* Generate Ensemble of Sequences */
						/* loop through number of networks in the ensemble */
						for(i = 0; i < _sizeEns; i++) {
							//printf("	en: %d\n", i);
							// generate graph
							ensedgcnt = GenerateGraph(_sizeNtw, &_seed, edgeProb, ensAdjTb, ensGenSq, ensAdjMtx, ensEdgeLst);
							//printf("\t%d\n", ensedgcnt);
							// get triangles and squares from
							L3(ensGenSq, ensAdjTb, ensAdjMtx, _sizeNtw, &ensTRI, ensLL3);
							L4(ensGenSq, ensAdjTb, ensAdjMtx, ensEdgeLst, _sizeNtw, ensedgcnt, &ensSQR);
							// write info to files
							fprintf(outTSG, "%d\t%d\n", ensTRI, ensSQR);
							for(j = 0; j < _sizeNtw; j++) {
								fprintf(outDL, "%d\t%d\n", ensGenSq[j], ensLL3[j]);
							}
							/*
							for(j = 0; j < _sizeNtw; j++) {
								for(k = 0; k < _sizeNtw; k++) {
									fprintf(outAdj, "%d\t", ensAdjMtx[j][k]);
								}
								fprintf(outAdj, "\n");
							}
							*/
						}

		/* MARK: CLEANUP AND END */
		/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */
				/* Close files */
				fclose(inSeq);
				fclose(outDL);
				fclose(outTSG);
				//fclose(outAdj);

				/* Clear arrays */
				free(seqIn);
				free(ensGenSq);
				free(ensLL3);
				free(seqQnty);
				free(seqCompr);
				free2d((void**)ensAdjTb);
				free2d((void**)ensAdjMtx);
				free2d((void**)ensEdgeLst);
				free(convSols);
				free(convFVEC);
				free2d((void**)edgeProb);

				/* End return */
				return 0;
	}

/* MARK: malloc2d */
	void** malloc2d(size_t numRows, size_t rowSize) {
		void **a;
		size_t i;

		/* a is an array of void * pointers that point to the rows */
		/* The last element is 0, so free2d can detect the last row */
		a = malloc(sizeof(void *) * (numRows + 1));        /* one extra for sentinel */
		if(a == 0) {
				/* malloc failed */
				return 0;
		}

		/* now allocate the actual rows */
		for(i = 0; i < numRows; i++) {
				a[i] = malloc(rowSize);
				if(a[i] == 0) {
						/* note that 0 in a[i] will stop freed2d after it frees previous rows */
						free2d(a);
						return 0;
				}
		}

		/* initialize the sentinel value */
		a[numRows] = 0;

		return a;
	}
/* MARK: free2d */
	void free2d(void **a) {
	    void **row;

	    /* first free rows */
	    for(row = a; *row != 0; row++) {
	        free(*row);
	    }

	    /* then free array of rows */
	    free(a);
	}
/* MARK: L3 */
	void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_L3_local) {
	    int i,j,k,a,b;
			*_L3_sum = 0;
	    for(k=0;k<_N;k++)
	    {

					_L3_local[k]=0;
					for(i=0;i<_zdegseq[k];i++)
	        {
	            a=_zat[k][i];
	            for(j=i+1;j<_zdegseq[k];j++)
	            {
	                b=_zat[k][j];
									*_L3_sum += _zam[a][b];
									_L3_local[k] += _zam[a][b];
	            }
	        }
	    }
	    *_L3_sum=*_L3_sum/3;
	    return;
	}
/* MARK: L4 */
	void L4(int *_zdegseq, int **_zat, int **_zam, int **_zel, int _N, int _edgenum, int *_L4_sum) {
	    int i,j,k,p,q,a,b;
			*_L4_sum = 0;
	    for(k=0;k<_edgenum;k++)
	    {
	        p=_zel[k][0];
	        q=_zel[k][1];

	        for(i=0;i<_zdegseq[p];i++)
	        {
	            a=_zat[p][i];
	            if(a==q)continue;
	            else
	            {
	                for(j=0;j<_zdegseq[q];j++)
	                {
	                    b=_zat[q][j];
	                    if(b==p)continue;
	                    else {
												//printf("\t%d\t%d\n", *_L4_sum, _zam[a][b]);
												*_L4_sum +=_zam[a][b];
												//printf("\t%d\t%d\n\n", *_L4_sum, _zam[a][b]);
											}
	                }
	            }
	        }
	    }
	    *_L4_sum=*_L4_sum/4;
	    //printf("L4=%d\n",L4_sum);
	    return;
	}
/* MARK: GenerateGraph */
	int GenerateGraph(int _N, long *_sd, float **_eP, int **_adt, int *_eSq, int**_aMx, int** _eLt) {
		/* Generates a single graph and places the results into an adjacency table, adjmatrix, generated sequence and edge list
				returns the number of edges generated
		*/
		int i, j, edgecount=0;
		/* MARK: Reset ensemble */
				for(i = 0; i < _N; i++) {
					_eSq[i] = 0;
					for(j = 0; j < _N; j++) {
						_adt[i][j] = -1;
						_aMx[i][j] = 0;
						//_eLt[i][j] = 0;
					}
				}
		/* MARK: Generate Graph */
			/* loop through all pairs of verticies in matrix to determine if edges exist */
			for(i = 0; i < _N; i++) {
				/* loop through all pairs of vertices in matrix to determine if edges exist */
				for(j = i+1; j < _N; j++) {
					/* generate random number, if less than probablity of edge then add edge to network */
					if(ran2(_sd) < _eP[i][j]) {
						// increment edge counter
						edgecount++;
						// populate adj table
						_adt[i][_eSq[i]] = j;
						_adt[j][_eSq[j]] = i;
						// populate generated sequence
						_eSq[i] += 1;
						_eSq[j] += 1;
						// populate adjmatrix
						_aMx[i][j] = 1;
						_aMx[j][i] = 1;
						// populate edge list
						_eLt[edgecount][0] = i;
						_eLt[edgecount][1] = j;
						//printf("\t%d\t%d\n", _eLt[edgecount][0], _eLt[edgecount][1]);
					}
				}
			}
			return edgecount;
	}
/* MARK: ImprovedConvergence */
	int ImprovedConvergence(int *_seqCompr, int *_seqQnty, int _seqcomplgth, float _gamma, float _tol, long *_seed, float *_convSol, float *_convFVEC, float *_convacc) {
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
					int itsMax = _seqcomplgth*10, itsCur = 0, convIts=100, i, j, prev;
					/* Declare flag variables */
					int errBst, tolBst, errB, tolB, errN, tolN, errOut;
					/* Allocate Memory */
					printf(" ALLOCATE ImprovedConvergence\n");
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

						for(j = prev; j < (prev+_seqQnty[i]); j++) {
							_convSol[j] = bestSol[i];
							_convFVEC[j] = bestFVEC[i];
						}
						prev += _seqQnty[i];
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
