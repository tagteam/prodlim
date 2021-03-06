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
		      double* risk,
		      double* cause_hazard,
		      double* varrisk,
		      double* I, /* current cumulative incidence */ 
		      double*I_lag, /* time lagged cumulative incidence */ 
		      double* v1,
		      double* v2,
		      int *t,
		      int start,
		      int stop) {
  
  
  int i,j,s,d,d1;
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
      /* {{{ loop over causes: compute risk */
      for(j=0; j < (*NS); ++j) {
	cause_hazard[s * (*NS) + j] = (event[s * (*NS) + j] / n);
	I_lag[j] = I[j];
	I[j] += S * cause_hazard[s * (*NS) + j];
	risk[s * (*NS) + j] = I[j];
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
	/* d2 = d - d1; */
	v1[j] += I[j] * (d / (n * (n - d))) + (S_lag * d1) / (n * n);
	v2[j] += (I[j] * I[j]) * (d / (n * (n - d)))
	  + ((S_lag * S_lag) * (n - d1) * d1) / (n * n * n)
	  + (2 * I[j] * S_lag * d1) / (n * n);
	varrisk[s * (*NS) + j] = (I[j] * I[j]) *  varH - 2 * I[j] * v1[j] + v2[j];
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
			 double* risk,
			 double* cause_hazard,
			 double* varrisk,
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
  
  
  int i,e,j,s,d,d1,entered;
  double S,S_lag,H,varH,atrisk;

  /* }}} */

  /* {{{ initialization */
  s=(*t);
  e=0;
  S=1;
  S_lag=1;
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
      /* do not initialize with those that have entry times=0 because by convention in case of ties 
	 entry happens after events and after censoring  */
      /* sort the delayed entry times */
      qsort(entrytime+start,
	    (stop-start),
	    (size_t) sizeof(double),
	    (int (*)(const void *, const void *))(doubleNewOrder));
      e=start; /* index for delayed entry */
      /* initialize the number at risk */

    }else{
      atrisk=(double) stop-start; /* (sub-)sample size */
    }
  }
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
      if (*delayed==1){
	/* delayed entry: find number of subjects that
	   entered at time[s] */
	entered=0;
	while(e<stop && entrytime[e]< y[i-1]){ /*entry happens at t+ and events at t*/
	  entered++;
	  /* FIXED: 17 Nov 2019 (12:50) */
	  if ((e+1>=stop) || entrytime[e] < entrytime[e+1]){ /* it has to be the last tie of the current entry time */
	  /* WRONG: 17 Nov 2019 (12:50) if (e==start || entrytime[e]>entrytime[e-1]){ */
	    /* unless there is a tie between the current
	       and the next entry-time, add time to list of times, increase s
	       and move the values of event, loss etc. to the next event time */
	    nrisk[s]=atrisk+entered;
	    if (s==0 || entrytime[e]!=time[s-1]){
	      /* if entrytime[e]==time[s-1] then only increase
		 the number at risk (done two lines above)
		 but dont change the time counter or the values
		 of event, etc.
	      */
	      for(j=0; j < (*NS); ++j) { 
		event[(s+1) * (*NS) + j]=event[s * (*NS) + j];
		event[s * (*NS) + j]=0;
	      }
	      loss[s+1]=loss[s];
	      loss[s]=0;
	      if (entrytime[e]<y[start]){
		surv[s]=1;
		for(j=0; j < (*NS); ++j) {
		  risk[s * (*NS) + j]=0;
		  varrisk[s * (*NS) + j]=0;
		}
	      } else{
		surv[s]=S_lag;
		for(j=0; j < (*NS); ++j) {
		  risk[s * (*NS) + j]=risk[(s-1) * (*NS) + j];
		  varrisk[s * (*NS) + j]=varrisk[(s-1) * (*NS) + j];
		}
	      }
	      time[s]=entrytime[e];
	      s++;
	    } 
	  }
	  e++;/* increase cumulative counter  */
	}
	atrisk += (double) entered;
      }
      time[s]=y[i-1];
      nrisk[s]=atrisk;
      d = 0;
      /* }}} */
      /* {{{ loop over causes: compute risk */
      for(j=0; j < (*NS); ++j) {
	cause_hazard[s * (*NS) + j] = (event[s * (*NS) + j] / atrisk);
	I_lag[j] = I[j];
	I[j] += S * cause_hazard[s * (*NS) + j];
	risk[s * (*NS) + j] = I[j];
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
	/* d2 = d - d1; */
	v1[j] += I[j] * (d / (atrisk * (atrisk - d))) + (S_lag * d1) / (atrisk * atrisk);
	v2[j] += (I[j] * I[j]) * (d / (atrisk * (atrisk - d)))
	  + ((S_lag * S_lag) * (atrisk - d1) * d1) / (atrisk * atrisk * atrisk)
	  + (2 * I[j] * S_lag * d1) / (atrisk * atrisk);
	varrisk[s * (*NS) + j] = (I[j] * I[j]) *  varH - 2 * I[j] * v1[j] + v2[j];
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


