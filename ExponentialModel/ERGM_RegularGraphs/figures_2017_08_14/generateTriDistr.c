#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ran2.h"

int main(int argc, char **argv) {

	FILE *file;
	int i, j, k, numtri, bytadd;
	long _seed = atol(argv[1]);
	char nameresult[200]={}, pathStrg[200]={};
	strcpy(pathStrg, argv[2]);

	for(i = 0; i < 316; i++) {
		//create filename
		bytadd = sprintf(nameresult, pathStrg);
		sprintf(bytadd+nameresult, "tridist_%d.txt", i);

		printf("%s\n", nameresult);
		if((file = fopen(nameresult,"w")) == NULL) {
			printf("opening file failed\n");
			return 1;
		}
		for(j = 0; j < 1000; j++) {
			numtri = 0;
			for(k = 0; k < 3*5209260; k++) {
				if(ran2(&_seed) < (double)(i*(i-1)*(i-1))/(double)(315*314*314)) {
					numtri++;
				}
			}
			fprintf(file, "%f\n", (double)numtri/3.0);
		}
		fclose(file);
	}

	return 0;
}
