/* Headers to include */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


/* Global variables and structure definition */
typedef struct {
    int *degree;
    int *label;
    int *forbidden;
} sequence;
typedef struct {
	int **list;
	double weight;
} graph;
static int orlen, len, offset, *allowed, *xk, *newseq;
static double minposw = 0;
static sequence orig, indseq;
static graph sample;

static int mt;//z_my_dd4m
/* Prototypes */
              void   gsaminit(const int *seq, const int n);
              void   gsamclean(void);
              graph  gsam    (double (*rng)(void), double (*AnyWeight)(double doffset, double dj, double m));
              //double MyWeight(double doffset, double dj);//z_my_dd4m
              //double MyWeight_dd2m(double doffset, double dj);//z_my_dd2m
static inline void   condsort(const int t);
static        void   build   (double (*rng)(void), double (*AnyWeight)(double doffset, double dj, double m));
static        void   allow   (int *alll, unsigned long int *restubs);


/* Inlined functions and macros */
static inline void condsort(const int t)
{
	int min, max, found=0, seek, ind, dummy;

	if (t<len-1) {															// If the node is not the last
		seek = indseq.degree[offset+t+1];
		if (indseq.degree[offset+t]<seek) {									// then if its degree is less than the degree of the next node
			if (indseq.degree[offset+t]==0) {								// check if it exhausted its stubs:
				ind=len-1;													// If so, erase it and move the last node in its place;
				len--;
				indseq.degree   [offset+t] = 1;
				indseq.label    [offset+t] = indseq.label    [offset+ind];
				indseq.forbidden[offset+t] = indseq.forbidden[offset+ind];
			} else {														// otherwise (if it still has stubs)
				min = t+1;													// find the last node whose degree is the same of the one of the next node...
				max = len-1;
				ind = min;
				do {
					if (max==min) found = 1;
					else {
						if (indseq.degree[offset+ind]==seek) {
							if (ind==max) found = 1;
							else {
								min = ind;
								ind = ceil((max+min)/2.0);
							}
						} else {
							max = ind-1;
							ind = ceil((max+min)/2.0);
						}
					}
				} while (!found);
				dummy = indseq.degree[offset+t];							// ... swap them...
				indseq.degree[offset+t] = seek;
				indseq.degree[offset+ind] = dummy;
				dummy = indseq.label[offset+t];
				indseq.label[offset+t] = indseq.label[offset+ind];
				indseq.label[offset+ind] = dummy;
				indseq.forbidden[offset+t] = indseq.forbidden[offset+ind];
				indseq.forbidden[offset+ind] = 1;							// and mark the node forbidden.
			}
		} else indseq.forbidden[offset+t] = 1;								// If instead the degree is not less than the degree of the next node, just mark the node forbidden.
	} else {																// If instead the node was the last,
		if (indseq.degree[offset+t]==0) len--;								// if it has exhausted its stubs, just delete it;
		else indseq.forbidden[offset+t] = 1;								// otherwise, mark it forbidden.
	}

	return;
}
