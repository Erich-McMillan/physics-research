// Erich McMillan
// 10 June 2016

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "InvTrans.h"

int main(int argc, char **argv) {
  if(argc != 6) {
    printf("Error testInvTrans. Improper number of arguments. Exiting to System:\n");
    return 1;
  }

  int _size, _seqNum, i, j, *seq, *dist;
  float _gamma, C;
  long _seed;
  char *_filepath, filename[100];
  FILE *fp;

  _filepath = (char*) calloc(200, sizeof(char));

  _size = atoi(argv[1]);
  _gamma = atof(argv[2]);
  _seqNum = atoi(argv[3]);
  _seed = atol(argv[4]);
  strcpy(_filepath, argv[5]);

  /* create filename */
  sprintf(filename, "/degseq_N%dr%.3lfSq%d.txt", _size, _gamma, _seqNum);
  strcat(_filepath, filename);

  /* open file */
  if((fp = fopen(_filepath, "w")) == NULL) {
    printf("File could not be opened. Exiting program:\n");
    return 1;
  }

  seq = (int*) malloc(_size*sizeof(int));
  dist = (int*) calloc(_size, sizeof(int));


  /* main process */
    C = InvTrans(_size, _gamma, &_seed, seq);
    //C = InvTransNoCheck(_size, _gamma, &_seed, seq);
    //printf("hello\n");
		/* print to file */
		for(i = 0; i < _size; i++) {
			fprintf(fp, "%d\n", seq[i]);
		}

  //printf("hello\n");

	fclose(fp);

  return 0;
}
