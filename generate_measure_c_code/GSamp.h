/* Structure definition */
typedef struct {
	int **list;
	double weight;
} graph;


/* Prototypes */
void  gsaminit (const int *seq, const int n);
void  gsamclean(void);
graph gsam     (double (*rng)(void), double (*AnyWeight)(double doffset, double dj, double m));
