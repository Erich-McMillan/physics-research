#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void** malloc2d(size_t numRows, size_t rowSize);
void free2d(void **a);

void loadNextFile(int, int, int, float, char*, float**, int**);
int possiblePairs(int);

int main(int _argc, char **_argv) {

	/*
		Description:
			Finds the average values of the local clustering coefficients by degree in the ensemble of sequences and writes them to the file.
		Inputs:
			_argv[1]: pathtofiles (string)
			_argv[2]: numfiles (int)
			_argv[3]: gamma (float)
			_argv[4]: ntwsize (int)
			_argv[5]: enssize (int)
			_argv[6]: pathtostorage (string)
		Outputs:
	*/

	/* MARK: Variable Declaration */
	char pathtofiles[200]={}, pathtostorage[200]={}, filename[200]={}, filepath[200]={};
	int numfiles, ntwsize, enssize, i, j, nextval;
	int **lccGenrDegreeCnt;
	float gamma;
	float **lccGenrDegree;
	FILE *genrlcc;

	/* MARK: Help Display */
	if(_argc < 2) {
		printf("\n\tHELP:\n");
		printf("\t_argv[1]: pathtofiles (string)\n\t_argv[2]: numfiles (int)\n\t_argv[3]: gamma (float)\n\t_argv[4]: ntwsize (int)\n\t_argv[5]: enssize (int)\n\t_argv[6]: pathtostorage (string)\n\n");
		return -1;
	}

	/* MARK: Read Command-Line Inputs */
	strcpy(pathtofiles, _argv[1]);
	numfiles = atoi(_argv[2]);
	gamma = atof(_argv[3]);
	ntwsize = atoi(_argv[4]);
	enssize = atoi(_argv[5]);
	strcpy(pathtostorage, _argv[6]);

	/* MARK: Allocate Memory */
	lccGenrDegree = (float**) malloc2d(ntwsize+1, numfiles * sizeof(float));
	lccGenrDegreeCnt = (int**) malloc2d(ntwsize+1, numfiles * sizeof(int));
	for(i = 0; i < ntwsize+1; i++) {
		for(j = 0; j < numfiles; j++) {
			lccGenrDegree[i][j] = 0.0;
			lccGenrDegreeCnt[i][j] = 0.0;
		}
	}

	//printf("%d\n", possiblePairs(316));
	/* MARK: Begin Parsing Data */
	for(i = 1; i <= numfiles; i++) {
		loadNextFile(i, ntwsize, enssize, gamma, pathtofiles, lccGenrDegree, lccGenrDegreeCnt);
	}

	/* MARK: Write data to file */
	strcpy(filepath, pathtostorage);
	sprintf(filename, "HistoLccGenr_%d_%.4d_%.3f_%d.txt", ntwsize, numfiles, gamma, enssize);
	strcat(filepath, filename);
	printf("%s\n", filepath);

	if((genrlcc = fopen(filepath, "w")) == NULL) {
		printf("Error opening %s exiting to system\n", filepath);
		return 1;
	}

	//{avg degree by network number}
	for(i = 0; i < numfiles; i++) {
		for(j = 0; j < ntwsize+1; j++) {
			printf("%d\t%d\t%f\n",j, lccGenrDegreeCnt[j][i], lccGenrDegree[j][i]);
			fprintf(genrlcc, "%f\t", lccGenrDegree[j][i]/lccGenrDegreeCnt[j][i]);
		}
		fprintf(genrlcc,"\n");
	}

	fclose(genrlcc);

	return 0;
}


void loadNextFile(int _i, int _ntwsize, int _enssize, float _gamma, char *_pathgenr, float **_gnLcc, int **_geCnt) {
	// Declare Variables
	char filepath[200]={}, filename[200]={};
	FILE *genr;
	int i = 0, j = 0, generated, numtris;
	float templcc;

	// Open files
	strcpy(filepath, _pathgenr);
	sprintf(filename, "DegreeLCC_%d_%.4d_%.3f_%d.txt", _ntwsize, _i, _gamma, _enssize);
	strcat(filepath,filename);
	printf("%s\n", filepath);
	if((genr = fopen(filepath, "r")) == NULL) {
	 	printf("Error opening %s exiting to system\n", filepath);
	 	return;
	}

	// Begin Parsing Data
	for(i = 0; i < 10000; i++) {
		for(j = 0; j < _ntwsize; j++) {
			fscanf(genr, "%d", &generated);
			fscanf(genr, "%d", &numtris);

			if(generated > 1) {
				templcc = numtris/(float)possiblePairs(generated);
			} else {
				templcc = 0;
			}
			if(templcc < 0 || templcc > 1) {
				printf("templcc\t%d\t%d\t%d\t%d\t%f\t%d\n", i, j, generated, numtris, templcc, possiblePairs(generated));
			}

			_gnLcc[generated][_i-1] += templcc;
			_geCnt[generated][_i-1]++;
		}
	}

	fclose(genr);

	// End Return
	return;
}
int possiblePairs(int _degree) {
	/*
		Finds the number of possible pairs between a node and its neighbours using a C[2,_degree] = _degree!/(degree!(_degree-2)!).
	*/

	long int nfact=1, nrfact=1, i;

	if(i > 2) {
		for(i = _degree-1; i <= _degree; i++) {
			//printf("%ld\n", i);
			nfact *= i;

		}
	} else if(i == 2) {
		return 1;
	} else {
		return 0;
	}
	//printf("improve program execution by just storing all the possibilities from 1 to 316")
	//printf("\t%ld\n", nfact);
	return nfact/2;

}
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
void free2d(void **a) {
		void **row;

		/* first free rows */
		for(row = a; *row != 0; row++) {
				free(*row);
		}

		/* then free array of rows */
		free(a);
}
