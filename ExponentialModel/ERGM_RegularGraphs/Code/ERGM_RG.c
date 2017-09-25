#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ran2.h"

int GenerateGraph(int _N, long *_sd, float _eP, int **_adt, int *_eSq, int**_aMx, int** _eLt);
void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_triples);
void** malloc2d(size_t numRows, size_t rowSize);
void free2d(void **a);

int main(int argc, char** argv) {
	/* declare variables */
	int N, R, i, j, deg=0, numtri, bytadd, numTriples=0;
	long seed;
	int *ensGenSq;
	int **ensAdjTb, **ensAdjMtx, **ensEdgeLst;
	double beta, prob;
	FILE *fresult;
	char nameresult[200]={}, pathStrg[200]={};


	if(argc != 5) {
		printf("first argument: Network size (N)\n");
		printf("second argument: Number of networks to generate (R)\n");
		printf("third argument: Seed for random number generator-must be negative integer (seed)\n");
		printf("fourth argument: location for storage of results (pathStrg)\n");
		return 0;
	}

	/* command line inputs */
	N = atoi(argv[1]);
	R = atoi(argv[2]);
	seed = atol(argv[3]);
	strcpy(pathStrg, argv[4]);

	/* create filename */
	sprintf(nameresult, "ERGM_RG_N%dr%d.txt", N, R);

	/* open file for writing */
	if((fresult=fopen(nameresult, "w"))==NULL) {
		printf("error opening file: exiting\n");
		return 1;
	}

	/* allocate memory */
	ensEdgeLst = (int**) malloc2d(N * N, sizeof(int) * 2);
	ensGenSq   = (int*) malloc(sizeof(int) * N);
	ensAdjTb   = (int**) malloc2d(N, sizeof(int) * N);
	ensAdjMtx  = (int**) malloc2d(N, sizeof(int) * N);

	/* get the avg degree for every possible degree between 1 and N-1*/
	for(j = 0; j < N; j++) {
		/* create filename */
		bytadd = sprintf(nameresult, pathStrg);
		sprintf(bytadd+nameresult, "ERGM_RG_N%dr%dd%.3d.txt", N, R,j);

		printf("%s\n", nameresult);

		/* open file for writing */
		if((fresult=fopen(nameresult, "w"))==NULL) {
			printf("error opening file: exiting\n");
			return 1;
		}

		/* calculate edge probability and Beta values */
		beta = log((N-1)/(double)j - 1.0);
	  prob = 1.0/(1.0+exp(beta));

		for(i = 0; i < R; i++) {
			/* generate graph */
			deg = GenerateGraph(N, &seed, prob, ensAdjTb, ensGenSq, ensAdjMtx, ensEdgeLst);

			numtri = 0;
			/* get number of triangles */
			L3(ensGenSq, ensAdjTb, ensAdjMtx, N, &numtri, &numTriples);

			/* write result to file */
			fprintf(fresult, "%f\t%d\t%f\n", 2*deg/(double)(N), numtri, (3*numtri)/(double)numTriples);
		}

		/* close file */
		fclose(fresult);

	}

	/* free arrays */
	free2d((void**)ensEdgeLst);
	free2d((void**)ensAdjTb);
	free2d((void**)ensAdjMtx);
	free(ensGenSq);

	return 0;
}


/* MARK: GenerateGraph */
int GenerateGraph(int _N, long *_sd, float _eP, int **_adt, int *_eSq, int**_aMx, int** _eLt) {
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
				if(ran2(_sd) < _eP) {
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
/* MARK: L3 */
void L3(int *_zdegseq, int **_zat, int **_zam, int _N, int *_L3_sum, int *_numTriples) {
  int i,j,k,a,b;
	*_L3_sum = 0;
	*_numTriples = 0;

  for(k=0;k<_N;k++)
  {
		for(i=0;i<_zdegseq[k];i++)
    {
      a=_zat[k][i];
      for(j=i+1;j<_zdegseq[k];j++)
      {
				*_numTriples += 1;
        b=_zat[k][j];
				*_L3_sum += _zam[a][b];
      }
    }
  }
  *_L3_sum=*_L3_sum/3;
  return;
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
