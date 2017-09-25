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

#include "ran2.h"


/* MARK: Function Declarations */
int GenerateGraph(int _N, long *_sd, float **_eP, int **_adt, int *_eSq, int**_aMx, int** _eLt);
void** malloc2d(size_t numRows, size_t rowSize);
void free2d(void **a);
void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_L3_local, int *_triples);
void L4(int *_zdegseq, int **_zat, int **_zam, int **_zel, int _N, int _edgenum, int *_L4_sum);

/* MARK: main */
int main (int argc, char **argv) {
	/*
			Description:

			Inputs:

			Outputs:
	*/

	/* Initialization */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */

		/* declare argument variables */
		int argVal, _sizeNtw, _sizeEns, _seqNum;
		float _gamma;
		long _seed;
		char _pathSeq[200]={}, _pathStrg[200]={};

		/* Parse arguments from command line */
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

		printf(" _sizeNtw = %d\n _sizeEns = %d\n _gamma = %f\n _seed = %ld\n",
		_sizeNtw, _sizeEns, _gamma, _seed);
		printf(" _pathSeq = %s\n _pathStrg = %s\n", _pathSeq, _pathStrg);

		/* declare numeric variables */
		int i, j, k, s, bytadd, tmpdeg, tmpqnty, ensedgcnt, ensGCC, ensTRI, ensSQR, seqcomplgth, expandcnt=0;
		int *seqIn, *ensGenSq, *ensGenAvgSq, *ensLL3, numTriples=0;
		int **ensAdjTb, **ensAdjMtx, **ensEdgeLst;
		char nameOutDL[200]={}, nameoutTS[200]={}, nameOutAdj[200]={};
		float convacc, tmpsol, tmpfvec;
		float *convSols;
		float **edgeProb;
		/* Declare flag variables */
		int convcode; // the returned value from ImprovedConvergence() if 0 then no convergence error
		/* Declare filevariables */
		FILE *inSeq, *outDL, *outTS/*, *outAdj*/;

		/* allocate Memory */
		seqIn      = (int*)     malloc(  sizeof(int) * _sizeNtw);
		ensGenSq   = (int*)     malloc(  sizeof(int) * _sizeNtw);
		ensGenAvgSq= (int*)     calloc(  _sizeNtw, sizeof(int) );
		ensLL3     = (int*)     malloc(  sizeof(int) * _sizeNtw);
		ensEdgeLst = (int**)    malloc2d(_sizeNtw * _sizeNtw , sizeof(int) * 2);
		ensAdjTb   = (int**)    malloc2d(_sizeNtw            , sizeof(int) * _sizeNtw);
		ensAdjMtx  = (int**)    malloc2d(_sizeNtw            , sizeof(int) * _sizeNtw);
		convSols   = (float*)   malloc(  sizeof(int) * _sizeNtw);
		edgeProb   = (float**)  malloc2d(_sizeNtw            , sizeof(float) * _sizeNtw);

		/* Open files for writing */
		bytadd = sprintf(nameOutDL, _pathStrg);
		sprintf(nameOutDL+bytadd, "/DegreeLCC_%d_%d_%.1f_%d.txt",_sizeNtw, _seqNum, _gamma, _sizeEns);
		bytadd = sprintf(nameoutTS, _pathStrg);
		sprintf(nameoutTS+bytadd, "/TriSqrGcc_%d_%d_%.1f_%d.txt",_sizeNtw, _seqNum, _gamma, _sizeEns);
		if((inSeq = fopen(_pathSeq, "r")) == NULL) {
			printf("Error opening seqFile exiting to system:\n");
			return 2;
		}
		if((outDL = fopen(nameOutDL, "w")) == NULL) {
			printf("Error opening file %s\n", nameOutDL);
			return 2;
		}
		if((outTS = fopen(nameoutTS, "w")) == NULL) {
			printf("Error opening file %s\n", nameoutTS);
			return 2;
		}



	/* Main code */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */

		/* Generate/Retrieve Sequence */
		rewind(inSeq);
		/* load sequence from file */
		while(fscanf(inSeq, "%d\t%f\t%d\t%f", &tmpdeg, &tmpsol, &tmpqnty, &tmpfvec) != -1) {

			for(j = expandcnt ; j < expandcnt + tmpqnty; j++ ) {
				seqIn[j] = tmpdeg;
				convSols[j] = tmpsol;
			}

			expandcnt = j;
		}

		/* Calculate probablities of edges */
		for(j = 0; j < _sizeNtw; j++) {
			for(k = 0; k < _sizeNtw; k++) {
				edgeProb[j][k] = 1.0/(1.0+exp((double)(-1.0*(convSols[j]+convSols[k]))));
			}
		}

		/* loop through number of networks in the ensemble */
		for(i = 0; i < _sizeEns; i++) {
			// generate graph
			ensedgcnt = GenerateGraph(_sizeNtw, &_seed, edgeProb, ensAdjTb, ensGenSq, ensAdjMtx, ensEdgeLst);

			for(j = 0; j < _sizeNtw*_sizeNtw; j++) {
				//printf("%d\t%d\n", ensEdgeLst[j][0], ensEdgeLst[j][1]);
			}

			for(j = 0; j < _sizeNtw; j++) {
				ensGenAvgSq[j] += ensGenSq[j];
				//printf("%d\n", ensGenSq[j]);
			}

			// get triangles and squares from
			L3(ensGenSq, ensAdjTb, ensAdjMtx, _sizeNtw, &ensTRI, ensLL3, &numTriples);
			//L4(ensGenSq, ensAdjTb, ensAdjMtx, ensEdgeLst, _sizeNtw, ensedgcnt, &ensSQR);
			ensSQR = -1;

			// write info to files
			fprintf(outTS, "%d\t%d\t%d\t%f\n", ensedgcnt, ensTRI, ensSQR, 3*ensTRI/(double)numTriples);
			for(j = 0; j < _sizeNtw; j++) {
				fprintf(outDL, "%d\t%d\n", ensGenSq[j], ensLL3[j]);
			}

		}

		// for(i = 0; i < _sizeNtw; i++) {
		// 	//printf("%d\t%f\n", seqIn[i], ensGenAvgSq[i]/(float)_sizeEns);
		// }



	/* Cleanup and end */
	/* ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ */

			/* Close files */
			fclose(inSeq);
			fclose(outDL);
			fclose(outTS);

			/* Clear arrays */
			free(seqIn);
			free(ensGenSq);
			free(ensLL3);
			free2d((void**)ensAdjTb);
			free2d((void**)ensAdjMtx);
			free2d((void**)ensEdgeLst);
			free(convSols);
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
void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_L3_local, int *_numTriples) {
    int i,j,k,a,b;
		*_L3_sum = 0;
		*_numTriples = 0;
    for(k=0;k<_N;k++)
    {

				_L3_local[k]=0;
				for(i=0;i<_zdegseq[k];i++)
        {
            a=_zat[k][i];
            for(j=i+1;j<_zdegseq[k];j++)
            {
                b=_zat[k][j];
								*_numTriples += 1;
								//printf("%d\n", *_numTriples);
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
		//printf("STARTED\n");
    for(k=0;k<_edgenum;k++)
    {
			//printf("\tK LEVEL\n");
      p=_zel[k][0];
      q=_zel[k][1];
			//printf("\tp=%d\tq=%d\n", p, q);

      for(i=0;i<_zdegseq[p];i++)
      {
				//printf("\tI LEVEL\n");
        a=_zat[p][i];
				//printf("\ta=%d\n", a);
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
				}
			}

			for(i = 0; i < _N*_N; i++) {
				_eLt[i][0] = 0;
				_eLt[i][1] = 0;
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
