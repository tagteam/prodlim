/*
  (2011) Thomas A. Gerds 
  --------------------------------------------------------------------
  distributed under the terms of the GNU public license 
*/
	      
#include <math.h>
#include <R.h>

void loo_surv(double *Y,
	      double *D,
	      double *time,
	      double *obsT,
	      double *status,
	      double *S, 
	      int *N,
	      int *NT){
  int k, t;
  double na,pl;
  for (k=0; k<*N;k++){
    /* compute the Nelson-Aalen estimate */
    pl=1;
    for (t=0; t<*NT;t++){
      if (obsT[k]>time[t]){
	/* decrease the number at risk */
	na = D[t]/(Y[t]-1);
      }
      else{
	if (obsT[k]==time[t]){
	  /*
	    maybe decrease the number of events,
	    but not the number at risk
	  */
	  na = (D[t]-status[k])/Y[t];
	}
	else{
	  /* do nothing */
	  na = D[t]/Y[t];
	}
      }
      /* compute the product-limit estimate */
      pl *= (1-na);
      S[k+(*N)*t]=pl;
    }
  }
}

void loo_comprisk(double *Y,
		  double *D,
		  double *time,
		  double *obsT,
		  double *status,
		  double *lagSurv,
		  double *F, 
		  int *N,
		  int *NT){
  int k, t;
  double na,aj;
  for (k=0; k<*N;k++){
    /* compute the Nelson-Aalen estimate */
    aj=0;
    for (t=0; t<*NT;t++){
      if (obsT[k]>time[t]){
	/* decrease the number at risk */
	na = D[t]/(Y[t]-1);
      }
      else{
	if (obsT[k]==time[t]){
	  /*
	    maybe decrease the number of events,
	    but not the number at risk
	  */
	  na = (D[t]-status[k])/Y[t];
	}
	else{
	  /* do nothing */
	  na = D[t]/Y[t];
	}
      }
      /* compute the Aalen-Johansen estimate */
      aj += lagSurv[t * (*N) + k] * na;
      F[k+(*N)*t]=aj;
    }
  }
}

