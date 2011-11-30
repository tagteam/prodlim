#include <R.h>
void predict_individual_survival(double *pred,
				 double *surv,
				 double *jumptime,
				 double *Y,
				 int *first,
				 int *size,
				 int *n,
				 int *lag){
  int j,i; /* start at index 0 */

  /* predicted survival probabilities at or just before the
     individual event times Y[i] */  
  for (i=0;i<(*n);i++){
    j=0;
    /* index j is in stratum i if j < size[i] */
    while(j < size[i] - 1 &&
	  jumptime[first[i] - 1 + j] != Y[i])
      j++;
    if (j - *lag < 0)
      pred[i]=1;
    else
      pred[i] = surv[first[i] - 1 + j - *lag];
  }

}
