#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void** malloc2d(size_t numRows, size_t rowSize);
void free2d(void **a);

int main(int argc, char **argv) {
	/* MARK: INITIALIZATION */
		/* Declare argv variables */
		char _pathToFile[200], _pathToStrg[200];
		int _numSeq, _sizeNtw, _numEns, _sizeEns, _seqNum;
		float _gamma;

		/* Assign argv variables */
		strcpy(_pathToFile, argv[1]);
		_sizeNtw = atoi(argv[2]);
		_numEns = atoi(argv[3]);
		_sizeEns = atoi(argv[4]);
		_gamma = atof(argv[5]);
		_seqNum = atoi(argv[6]);
		strcpy(_pathToStrg, argv[7]);

		/* Display argv variables */
		printf("*******************************************************************************************");
		printf("\n FILE.exe\n");
		printf("\t_pathToFile = %s\n\t_sizeNtw = %d\n\t_numEns = %d\n\t_sizeEns = %d\n\t_gamma = %f\n",
					_pathToFile, _sizeNtw, _numEns, _sizeEns, _gamma);

		/* Declare function variables */
		FILE *adjmatrices, *output;
		int h, i, j, k, edge;
		char outfilename[200];
		/* Allocate memory */

		/* Open adjmatrices file */
			if((adjmatrices = fopen(_pathToFile, "r")) == NULL) {
				printf("ERROR opening %s exiting to system:\n", _pathToFile);
			}

			sprintf(outfilename, "/NAME_N%d_G%.3f_S%d_E%d.txt", _sizeNtw, _gamma, _seqNum, _sizeEns);
			strcat(_pathToStrg, outfilename);
			if((output = fopen(_pathToStrg, "w")) == NULL) {
				printf("Error opening file, exiting to system:\n");
				return 1;
			}

			fprintf(output,"Val\n");


			/* MARK: MAIN CODE */
				/* Set adjmatrices to zero */
					for(i = 0; i < _sizeNtw; i++) {
						seq[i] = 0;
					}
				/* loop through all networks in ensemble */
					for(h = 0; h < _sizeEns; h++) {
						/* load next network */
						for(i = 0; i < _sizeNtw; i++) {
							for(j = 0; j < _sizeNtw; j++) {
								fscanf(adjmatrices, "%d", &edge);
							}
						}
					}

				fclose(adjmatrices);
				fclose(output);
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
