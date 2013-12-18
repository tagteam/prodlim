#include <R.h>
void iindex(int *iindex,
	    int *strata,
	    double *L,
	    double *R,
	    double *U,
	    int *N,
	    int *NS){
  int s, i, k;
  k=0;
  for (s=0;s<(*NS-1);s++){
    i=0;
    for (i=0; i<*N;i++){
      if ((L[i]==R[i] && L[i]==U[s+1]) /* exact obs */
	  ||
	  (L[i]<U[s+1] && R[i]>U[s])){   /* [U[s],U[s+1]] intersects [L[i],R[i]] */
	iindex[k] = i+1;
	k++;
      }
    }
    strata[s]=k;
  }
}
