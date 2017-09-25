#include"GSampler.h"


/* Graph sampler initializer */
void gsaminit(const int *seq, const int nodes)
{
	int i, n=0;
	unsigned long int stubs = 0;
	
	do n++;
	while (n<nodes && seq[n]!=0);
	
	/* Memory allocations */
	orig.degree      = calloc(n,sizeof(int));								// Original sequence
	orig.label       = calloc(n,sizeof(int));								// Original node labels
	orig.forbidden   = calloc(n,sizeof(int));								// Original forbidden status labels
	indseq.degree    = calloc(n,sizeof(int));								// Work sequence
	indseq.label     = calloc(n,sizeof(int));								// Work sequence labels
	indseq.forbidden = calloc(n,sizeof(int));								// Work sequence forbidden status labels
	allowed          = calloc(n,sizeof(int));								// Allowed nodes
	newseq           = calloc(n,sizeof(int));								// D' sequence
	xk               = calloc(n,sizeof(int));								// Crossing indices hash table
	
	/* Initializtions of sequences */
	for (i=0; i<n; i++) {
		orig.degree[i]     = seq[i];
		orig.label[i]      = i;
		orig.forbidden[i]  = 0;
		stubs             += seq[i];
	}
	
	/* Initialization of sample */
	sample.list    = malloc(n*sizeof(int*));								// Adjacency list
	sample.list[0] = calloc(stubs,sizeof(int));
	for (i=1; i<n; i++) sample.list[i] = sample.list[i-1]+orig.degree[i-1];
	orlen = n;
	
	return;
}



/* Graph sampler cleanup */
void gsamclean(void)
{
	
	free(orig.degree);
	free(orig.label);
	free(orig.forbidden);
	free(indseq.degree);
	free(indseq.label);
	free(indseq.forbidden);
	free(allowed);
	free(newseq);
	free(xk);
	free(sample.list[0]);
	free(sample.list);
	
	return;
}



/* Graph sampler */
graph gsam(double (*rng)(void), double (*AnyWeight)(double doffset, double dj, double m))
{
	int i;
	mt=0;//z_my_dd4m
	len = orlen;
	for (i=0; i<len; i++) {
		indseq.degree[i]    = orig.degree[i];
		indseq.label[i]     = i;
		indseq.forbidden[i] = 0;
		mt                 += orig.degree[i];//z_my_dd4m
	}
	offset = 0;
	sample.weight = 0.0;
	
	mt=mt/2;//z_my_dd4m

	build(rng,AnyWeight);
	
	return sample;
}


/* Recursively places all the possible links in the graph */
static void build(double (*rng)(void), double (*AnyWeight)(double doffset, double dj, double m))
{
	int i, j, alll;
	unsigned long int ext, count=0, restubs=0;
	double Saberi_ext,Saberi_count=0, Saberi_Sum=0;//Saberi

	for (j=1; j<len; j++) Saberi_Sum += AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+j], (double)mt);//z_my_dd4m
	Saberi_ext = Saberi_Sum*rng();
	j = 0;
	do {
		j++;
		Saberi_count += AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+j], (double)mt);//z_my_dd4m
	} while (Saberi_count<=Saberi_ext);
	sample.weight += log(Saberi_Sum)-lgamma(indseq.degree[offset]+1)-log(AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+j], (double)mt));//z_my_dd4m									// Start computing the weight
	sample.list[indseq.label[offset]][orig.degree[indseq.label[offset]]-indseq.degree[offset]] = indseq.label[offset+j];			// Connect the node to the hub
	sample.list[indseq.label[offset+j]][orig.degree[indseq.label[offset+j]]-indseq.degree[offset+j]] = indseq.label[offset];
	indseq.degree[offset]--;																										// Reduce the degree of the hub
	indseq.degree[offset+j]--;														 												// Reduce tuhe degree of the target node
	condsort(j);																													// Reorder the sequence if needed.
	mt--;//z_my_dd4m

	while (indseq.degree[offset]>0) {																								// If the hub node still has stubs
		alll=0;																														// Reset the size of the allowed nodes set
		restubs=0;
		Saberi_Sum=0;//Saberi
		allow(&alll,&restubs);																										// Build the set of allowed nodes
		for (i=0; i<alll; i++) Saberi_Sum += AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+allowed[i]],(double)mt);//z_my_dd4m
		Saberi_ext = Saberi_Sum*rng();
		Saberi_count = 0;
		i = -1;
		do {
			i++;
			Saberi_count += AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+allowed[i]],(double)mt);//z_my_dd4m
		} while (Saberi_count<=Saberi_ext);
		j = allowed[i];
		sample.weight += log(Saberi_Sum)-log(AnyWeight((double)indseq.degree[offset], (double)indseq.degree[offset+j], (double)mt));//z_my_dd4m																// Update the weight
		sample.list[indseq.label[offset]][orig.degree[indseq.label[offset]]-indseq.degree[offset]] = indseq.label[offset+j];		// Connect the node to the hub
		sample.list[indseq.label[offset+j]][orig.degree[indseq.label[offset+j]]-indseq.degree[offset+j]] = indseq.label[offset];
		indseq.degree[offset]--;																									// Reduce the degree of the hub
		indseq.degree[offset+j]--;														 											// Reduce tuhe degree of the target node
		condsort(j);																												// Reorder the sequence if needed.
		mt--;//z_my_dd4m
	}
	
	if (len>1) {																													// If there are more nodes to be connected
		offset++;																													// Update the memory offset
		len--;																														// Reduce the length of the sequence
		for (j=0; j<len; j++) indseq.forbidden[offset+j] = 0;																		// Clear the forbidden flags
		build(rng,AnyWeight);																											// and keep sampling.
	}
	
	return;
}



/* Builds the set of the nodes which we are allowed to connect to */
static void allow(int *alll, unsigned long int *restubs)
{
	int stop=0, faildeg, count=0, lastdeg, flag=0, k=0, j=1, t, dumdeg, firstdeg, curdeg=-1, moven=0, tarpos, oripos, check, degchk, index, indexcopy;
	long int degsum, minsum;
	
	
    /* Start making a copy of the sequence.
       Add the leftmost adjacency set to the allowed node set.
       Also decrease the degrees of the nodes in the copy sequence
       which are added to the set and make them forbidden */
	while (count<indseq.degree[offset] && j<len) {
		*(newseq+j-1) = indseq.degree[offset+j];
		if (!indseq.forbidden[offset+j]) {
			(*(newseq+j-1))--;
			(*alll)++;
			*(allowed+(*alll)-1) = j;
			(*restubs) += indseq.degree[offset+j];
			count++;
			if (curdeg!=*(newseq+j-1)) {
				curdeg = (*(newseq+j-1));
				firstdeg = j-1;
			}
		} else {
			if (j>1) {
				if (*(newseq+j-1)>*(newseq+j-2)) {
					dumdeg = *(newseq+firstdeg);
					*(newseq+firstdeg) = *(newseq+j-1);
					*(newseq+j-1) = dumdeg;
					firstdeg++;
				}
			} else {
				firstdeg = 0;
				curdeg = *(newseq);
			}
		}
		j++;
	}
	lastdeg = indseq.degree[offset+(*(allowed+(*alll)-1))];												// Degree of the last node in the leftmost adjacency set
	
	/* Revert the degree (in the copy sequence) of the last node added to the set.
	   If needed, also add to the allowed set all the nodes that are not restricted and have
	   the same degree of the last node in the leftmost adjacency set and keep copying the sequence */
	check = (*(allowed+(*alll)-1))-1;
	(*(newseq+check))++;
	if (check>0) degchk = *(newseq+check-1);
	else degchk=-1;
	while (j<len && (indseq.degree[offset+j]==lastdeg || indseq.forbidden[offset+j])) {
		*(newseq+j-1) = indseq.degree[offset+j];
		if (!indseq.forbidden[offset+j]) {
			(*alll)++;
			*(allowed+(*alll)-1) = j;
			(*restubs) += indseq.degree[offset+j];
		}
		j++;
	}
	
	if (j<len) {																						// If there are more nodes left
	    for (index=j; index<len; index++) *(newseq+index-1) = indseq.degree[offset+index];				// Finish copying the sequence
	    *(newseq+len-1) = 1;																			// Hub at the end with degree 1
	    
	    /* If needed, quickly reorder the sequence */
		if (degchk!=-1) {
			oripos = check;
			while (check<len && degchk<*(newseq+check)) {moven++;check++;}
			if (moven>0) {
				degchk = *(newseq+oripos);
				tarpos = oripos-1;
				while (tarpos>=0 && degchk>*(newseq+(tarpos))) tarpos--;
				tarpos++;
				for (t=oripos-1; t>=tarpos; t--) *(newseq+t+moven)=*(newseq+t);
				for (t=0; t<moven; t++) *(newseq+tarpos+t) = degchk;
			}
		}
		
	    /* Build the hash table for the crossing indices */
	    indexcopy=-1;
	    for (index=len-1; index>=(*newseq); index--) *(xk+index) = 0;
	    for (index=1; index<len; index++) for (indexcopy=(*(newseq+index-1))-1; indexcopy>=(*(newseq+index)); indexcopy--) *(xk+indexcopy) = index;
	    for (index=indexcopy; index>=0; index--) *(xk+indexcopy) = len;
	    
	    /* Start the E-G test */
		degsum  = (*newseq);
		minsum  = len-1;
		faildeg = 0;																					// If L=R, then the maximum fail degree is the degree of the next node or 1 less than the
		if (degsum==minsum) {																			// degree of the last node in the leftmost adjacency set of the hub, whichever is smaller.
			if ((*(newseq+1))<lastdeg) faildeg = (*(newseq+1));											// If instead L=R-1, then the fail degree for this value of k is either the degree of the
			else faildeg = lastdeg-1;																	// node whose index is the crossing index for (k+1) or the degree of the next node,
			stop = 1;																					// depending on whether the crossing index for (k+1) is greater than k or not, respectively.
		} else if (degsum==(minsum-1)) if (lastdeg>1) faildeg = 1;										// For k=0, the case L=R-1 reduces to a fail degree of 1.
		
		if (len>=3) {																					// If there are more nodes, keep testing as long as necessary, as described above
			while ((!stop) && (k<(len-2))) {
				k++;
				degsum += (*(newseq+k));
				if ((*(xk+k))<k+1) flag = 1;
				if (!flag) minsum += (*(xk+k))-1;
				else minsum += (2*k-(*(newseq+k)));
				if (degsum==minsum) {
					if ((*(newseq+k+1))<lastdeg) faildeg = (*(newseq+k+1));
					else faildeg = lastdeg-1;
					stop = 1;
				} else if (degsum==(minsum-1)) {
					if ((*(xk+k+1))>k) {
						if ( (*(newseq+(*(xk+k+1))))>=faildeg ) faildeg = (*(newseq+(*(xk+k+1))));
						else stop = 1;
					} else {
						if ( (*(newseq+k+1))>=faildeg ) faildeg = (*(newseq+k+1));
						else stop = 1;
					}
				}
			}
		}
		
		/* Add to the allowed set all the non-forbidden nodes
		   whose degree is greater than the maximum fail degree */
		while (j<len) {
			if ( !indseq.forbidden[offset+j] && indseq.degree[offset+j]>faildeg ) {
				(*alll)++;
				*(allowed+(*alll)-1) = j;
				(*restubs) += indseq.degree[offset+j];
			}
			j++;
		}
	}
	
	return;
}
