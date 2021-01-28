/*
  (2011, 2020) Thomas A. Gerds 
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
	      double *loo,
	      int *N,
	      int *NT,
	      int *NP,
	      int *pos,
	      int *lag){
  int k, t, p;
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
      S[t]=pl;
    }
    for (p=0; p<*NP;p++){
    if (*lag==1){
      if (pos[p]<=1) loo[k+(*N)*p]=1; else loo[k+(*N)*p] = S[pos[p]-2];
    }else{
      if (pos[p]==0) loo[k+(*N)*p]=1; else loo[k+(*N)*p] = S[pos[p]-1];
    }
    }
  }
}


void loo_comprisk(double *Y,
		  double *D,
		  double *Dall,
		  double *time,
		  double *obsT,
		  double *status,
		  double *event,
		  double *S,
		  double *F,
		  double *loo,
		  int *N,
		  int *NT,
		  int *NP,
		  int *pos){
  int k, t, p;
  double na,naall,aj,pl;
  for (k=0; k<*N;k++){
    /* compute the Nelson-Aalen estimate */
    aj=0;
    pl=1;
    for (t=0; t<*NT;t++){
      if (obsT[k]>time[t]){
	/* decrease the number at risk
	   because k was in the risk set at time[t]
	*/
	naall = Dall[t]/(Y[t]-1);
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
	  naall = (Dall[t]-status[k])/(Y[t]-1);
	  na = (D[t]-event[k])/(Y[t]-1);
	}
	else{
	  /* do nothing */
	  naall = Dall[t]/Y[t];
	  na = D[t]/Y[t];
	}
      }
      /* calculate the event-free Kaplan-Meier estimate*/
      pl *= (1-naall);
      S[t]=pl;
      /* compute the Aalen-Johansen estimate */
      if (t==0)
	aj += na;
      else
	aj += S[t-1] * na;
      /* if (k==0) Rprintf("t=%d\tD[t]=%1.2f\tY[t]=%1.2f\tS[t]=%1.2f\tna=%1.2f\tnaall=%1.2f\t\n",t,D[t],Y[t],S[t],na,naall); */
      F[t]=aj;
    }
    for (p=0; p<*NP;p++){
      if (pos[p]==0) loo[k+(*N)*p]=1; else loo[k+(*N)*p] = F[pos[p]-1];
    }
  }
}

