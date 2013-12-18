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
    /* Rprintf("\n"); */
    /* compute the Nelson-Aalen estimate */
    pl=1;
    for (t=0; t<*NT;t++){
      if (obsT[k]>time[t]){
	/* decrease the number at risk
	   because individual k was at risk
	   at time[t]
	 */
	na = D[t]/(Y[t]-1);
      }
      else{
	if (obsT[k]==time[t]){
	  /*
	    decrease the number of events
	    if k was an event,
	    and decrease the number at risk
	    because k was in the risk set at
	    time[t]
	  */
	  na = (D[t]-status[k])/(Y[t]-1);
	}
	else{
	  /* do nothing */
	  na = D[t]/Y[t];
	}
      }
      /* compute the product-limit estimate */
      pl *= (1-na);
      S[k+(*N)*t]=pl;
      /* Rprintf("t=%d\tk=%d\tD[t]=%1.2f\tY[t]=%1.2f\tna=%1.2f\tS[k](t)=%1.2f\n",t,k,D[t],Y[t],na,S[k+(*N)*t]); */
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
	/* decrease the number at risk
	   because k was in the risk set at time[t]
	*/
	na = D[t]/(Y[t]-1);
      }
      else{
	if (obsT[k]==time[t]){
	  /*
	    decrease the number of events
	    if k was an event,
	    and decrease the number at risk
	    because k was in the risk set at
	    time[t]
	  */
	  na = (D[t]-status[k])/(Y[t]-1);
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

