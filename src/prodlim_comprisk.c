#include <math.h>
#include <R.h>
#include "prodlim.h"

/*

Compute the Aalen-Johannsen estimate in a loop over "NS" causes.

Important: the vector "cause" has code "-1" for censored obs

*/

/* {{{ Header */

void prodlim_comprisk(double* y, 
		      double* status,
		      int* cause,
		      int* NS, /* number of causes (states) */
		      double* time,
		      double* nrisk,
		      double* event, 
		      double* loss, 
		      double* surv,
		      double* cuminc,
		      double* cause_hazard,
		      double* varcuminc,
		      double* I, /* current cumulative incidence */ 
		      double*I_lag, /* time lagged cumulative incidence */ 
		      double* v1,
		      double* v2,
		      int *t,
		      int start,
		      int stop) {
  
  
  int i,j,s,d,d1,d2;
  double S,S_lag,H,varH,n;

  /* }}} */

  /* {{{ initialization */
  s=(*t);
  S=1;
  H=0;
  for(j=0; j < (*NS); ++j) {
    I[j]=0;
    I_lag[j]=0;
    v1[j]=0;
    v2[j]=0;
  }
  varH=0;
  n=(double) stop-start; /* (sub-)sample size */ 

  
  if (status[start]>0)
    event[s *(*NS) + cause[start]]=1;
  else
    loss[s]=1;
  /* }}} */
  
  for (i=(1+start);i<=stop;i++){
    /* {{{ if tie then wait */
    if (i<stop && y[i]==y[i-1]){
      if (status[i]>0)
	event[s * (*NS) + cause[i]] +=1;
      else
	loss[s]+=1;
    }
    /* }}} */
    else {
      /* {{{ at s: set time, atrisk; reset d */
      time[s]=y[i-1];
      nrisk[s]=n;
      d = 0;
      /* }}} */
      /* {{{ loop over causes: compute cuminc */
      for(j=0; j < (*NS); ++j) {
	cause_hazard[s * (*NS) + j] = (event[s * (*NS) + j] / n);
	I_lag[j] = I[j];
	I[j] += S * cause_hazard[s * (*NS) + j];
	cuminc[s * (*NS) + j] = I[j];
	d += event[s * (*NS) + j];
      }
      /* }}} */
      /* {{{ compute survival */
      S_lag = S;
      pl_step(&S, &H, &varH, n, d, 0);
      surv[s] = S;
      /* }}} */
      /* {{{ variance estimate Marubini & Valsecchi (1995), Wiley, chapter 10, page 341 */
      for (j=0; j < (*NS); ++j){
	d1 = event[s * (*NS) + j];
	d2 = d - d1;
	v1[j] += I[j] * (d / (n * (n - d))) + (S_lag * d1) / (n * n);
	v2[j] += (I[j] * I[j]) * (d / (n * (n - d)))
	  + ((S_lag * S_lag) * (n - d1) * d1) / (n * n * n)
	  + (2 * I[j] * S_lag * d1) / (n * n);
	varcuminc[s * (*NS) + j] = (I[j] * I[j]) *  varH - 2 * I[j] * v1[j] + v2[j];
	/* varH is greenwood's formula */
	/* variance estimate Korn & Dorey (1992), Stat in Med, Vol 11, page 815 */
	/* I1 = (I[j] - I_lag[j]) / 2; */
      }
      /* }}} */
      /* {{{ update atrisk, set n.event, loss, for the next time point */
      if (i<stop){
	n -= (d + loss[s]);
	s++;
	if (status[i]>0){
	  event[s *(*NS) + cause[i]]=1;
	}
	else
	  loss[s]=1;
      }
      /* }}} */
    }
  }
  *t=(s+1); /* for the next strata  */
}


void prodlimCompriskPlus(double* y, 
			 double* status,
			 int* cause,
			 double *entrytime,
			 double *caseweights,
			 int* NS, /* number of causes (states) */
			 double* time,
			 double* nrisk,
			 double* event, 
			 double* loss, 
			 double* surv,
			 double* cuminc,
			 double* cause_hazard,
			 double* varcuminc,
			 double* I, /* current cumulative incidence */ 
			 double* I_lag, /* time lagged cumulative incidence */ 
			 double* v1,
			 double* v2,
			 int *t,
			 int start,
			 int stop,
			 int *delayed,
			 int *weighted
			 ) {
  
  
  int i,e,j,s,d,d1,d2,entered;
  double S,S_lag,H,varH,atrisk;

  /* }}} */

  /* {{{ initialization */
  s=(*t);
  S=1;
  H=0;
  for(j=0; j < (*NS); ++j) {
    I[j]=0;
    I_lag[j]=0;
    v1[j]=0;
    v2[j]=0;
  }
  varH=0;
  if (*weighted==1){
    atrisk=0;
    for (i=start;i<stop;i++) atrisk += caseweights[i];
  } else{
    if (*delayed==1){
      atrisk=0;
      /* sort the delayed entry times */
      qsort(entrytime+start,
	    (stop-start),
	    (size_t) sizeof(double),
	    (int (*)(const void *, const void *))(doubleNewOrder));
    }else{
      atrisk=(double) stop-start; /* (sub-)sample size */
    }
  }
  i=0;
  e=0; /* index for delayed entry */
  if (*weighted==1){
    if (status[start]>0){
      event[s *(*NS) + cause[start]]=caseweights[start];
    } else{
      loss[s]=caseweights[start];
    }
  }
  else{
    if (status[start]>0){
      event[s *(*NS) + cause[start]]=1;
    } else{
      loss[s]=1;
    }
  }
  
  /* }}} */
  
  for (i=(1+start);i<=stop;i++){
    /* {{{ if tie then wait */
    if (i<stop && y[i]==y[i-1]){
      if (*weighted==1){
	if (status[i]>0)
	  event[s * (*NS) + cause[i]] +=caseweights[i];
	else
	  loss[s]+=caseweights[i];
      }
      else{
	if (status[i]>0)
	  event[s * (*NS) + cause[i]] ++;
	else
	  loss[s]++;
      }
    }
    /* }}} */
    else{
      /* {{{ at s: set time, atrisk; reset d */
      time[s]=y[i-1];
      if (*delayed==1){
	/* delayed entry: find number of subjects that
	   entered at time[s] */
	entered=0;
	while(e<stop && entrytime[e]<time[s]){
	  e++;
	  entered++;
	}
	atrisk += (double) entered;
	/* Rprintf("D1=%d\tD2=%d\tF3=%1.2f\tF4=%1.2f\t\n",entered,e,atrisk,entrytime[e]); */
      }
      nrisk[s]=atrisk;
      d = 0;
      /* }}} */
      /* {{{ loop over causes: compute cuminc */
      for(j=0; j < (*NS); ++j) {
	cause_hazard[s * (*NS) + j] = (event[s * (*NS) + j] / atrisk);
	I_lag[j] = I[j];
	I[j] += S * cause_hazard[s * (*NS) + j];
	cuminc[s * (*NS) + j] = I[j];
	d += event[s * (*NS) + j];
      }
      /* }}} */
      /* {{{ compute survival */
      S_lag = S;
      pl_step(&S, &H, &varH, atrisk, d, 0);
      surv[s] = S;
      /* }}} */
      /* {{{ variance estimate Marubini & Valsecchi (1995), Wiley, chapter 10, page 341 */
      for (j=0; j < (*NS); ++j){
	d1 = event[s * (*NS) + j];
	d2 = d - d1;
	v1[j] += I[j] * (d / (atrisk * (atrisk - d))) + (S_lag * d1) / (atrisk * atrisk);
	v2[j] += (I[j] * I[j]) * (d / (atrisk * (atrisk - d)))
	  + ((S_lag * S_lag) * (atrisk - d1) * d1) / (atrisk * atrisk * atrisk)
	  + (2 * I[j] * S_lag * d1) / (atrisk * atrisk);
	varcuminc[s * (*NS) + j] = (I[j] * I[j]) *  varH - 2 * I[j] * v1[j] + v2[j];
	/* varH is greenwood's formula */
	/* variance estimate Korn & Dorey (1992), Stat in Med, Vol 11, page 815 */
	/* I1 = (I[j] - I_lag[j]) / 2; */
      }
      /* }}} */
      /* {{{ update atrisk, set n.event, loss, for the next time point */
      if (i<stop){
	atrisk -= (d + loss[s]);
	s++;
	if (*weighted==1){
	  if (status[i]>0){
	    event[s *(*NS) + cause[i]]=caseweights[i];
	  }
	  else
	    loss[s]=caseweights[i];
	}
	else{
	  if (status[i]>0){
	    event[s *(*NS) + cause[i]]=1;
	  }
	  else
	    loss[s]=1;
	}
      }
      /* }}} */
    }
  }
  *t=(s+1); /* for the next strata  */
}


