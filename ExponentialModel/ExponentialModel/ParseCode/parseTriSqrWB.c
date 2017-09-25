#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int _argc, char **_argv) {

	/*
		Description:
			Finds the histogram of the #ofTriangles and #ofSquares in the ensemble of sequences and writes them to the file.
		Inputs:
			_argv[1]: pathtofiles (string)
			_argv[2]: numfiles (int)
			_argv[3]: gamma (float)
			_argv[4]: ntwsize (int)
			_argv[5]: enssize (int)
			_argv[6]: pathtostorage (string)
		Outputs:
			Will read in all the files in the directory specified which conform to the format: TriSqrGcc_$ntwsize_$numfiles_$gamma_$enssize.txt and write the histogram for tri and sqr to a file called HistoTriSqr_$ntwsize_$numfiles_$gamma_$enssize.txt in two columns the first for the triangle distribution, the second for the square distribution.
	*/

	/* MARK: Variable Declaration */
	char pathtofiles[200]={}, pathtostorage[200]={}, filename[200]={}, filepath[200]={}, filenametri[200]={}, filepathtri[200]={}, filenamesqr[200]={}, filepathsqr[200]={};
	int numfiles, ntwsize, enssize, largestval=0, i, nexttrival, nextsqrval;
	int *histoSqr, *histoTri;
	float gamma;
	FILE *curr_Ftri, *curr_Fsqr, *currfile;

	/* MARK: Read Command-Line Inputs */
	strcpy(pathtofiles, _argv[1]);
	numfiles = atoi(_argv[2]);
	gamma = atof(_argv[3]);
	ntwsize = atoi(_argv[4]);
	enssize = atoi(_argv[5]);
	strcpy(pathtostorage,_argv[6]);

	/* MARK: Allocate Memory */
	histoSqr = (int*) calloc(50000, sizeof(int));
	histoTri = (int*) calloc(50000, sizeof(int));

	/* MARK: Begin Parsing Data */
	for(i = 1; i <= numfiles; i++) {
		strcpy(filepathtri, pathtofiles);
		sprintf(filenametri, "Coolen_L3_N%dr%.1fSq%d.txt", ntwsize, -1*gamma, i);
		strcat(filepathtri,filenametri);
		printf("%s\n", filepathtri);
		//Coolen_L3_N100r2.0Sq821

		strcpy(filepathsqr, pathtofiles);
		sprintf(filenamesqr, "Coolen_L4_N%dr%.1fSq%d.txt", ntwsize, -1*gamma, i);
		strcat(filepathsqr,filenamesqr);
		printf("%s\n", filepathsqr);

		if((curr_Ftri = fopen(filepathtri, "r")) == NULL) {
			printf("Error opening %s exiting to system\n", filepathtri);
			return 1;
		}

		if((curr_Fsqr = fopen(filepathsqr, "r")) == NULL) {
			printf("Error opening %s exiting to system\n", filepathsqr);
			return 1;
		}

		while(fscanf(curr_Ftri, "%d", &nexttrival) != EOF && fscanf(curr_Fsqr, "%d", &nextsqrval) != EOF) {
			histoTri[nexttrival]++;
			histoSqr[nextsqrval]++;

			if(nexttrival > largestval) {
				largestval = nexttrival;
			} else if(nextsqrval > largestval) {
				largestval = nextsqrval;
			} else {

			}
		}

		fclose(curr_Ftri);
		fclose(curr_Fsqr);
	}

	strcpy(filepath, pathtostorage);
	sprintf(filename, "HistoTriSqrWB_%d_%.4d_%.3f_%d.txt", ntwsize, numfiles, gamma, enssize);
	strcat(filepath, filename);
	printf("%s\n", filepath);

	if((currfile = fopen(filepath, "w")) == NULL) {
		printf("Error opening %s exiting to system\n", filepath);
		return 1;
	}

	for(i = 0; i < largestval; i++) {
		fprintf(currfile, "%d\t%d\n", histoTri[i], histoSqr[i]);
	}

	fclose(currfile);

	return 0;
}
