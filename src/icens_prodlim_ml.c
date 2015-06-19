#include <math.h>
#include <R.h>

#define max(A,B) ((A) > (B) ? (A):(B))
#define min(A,B) ((A) < (B) ? (A):(B))

void icens_prodlim_ml(double *L,
		      double *R,
		      double *petoL,
		      double *petoR,
		      int *indexL,
		      int *indexR,
		      int *status,
		      double *N,
		      double *NS,
		      double *nrisk,
		      double *nevent,
		      double *ncens,
		      double *hazard,
		      double *var_hazard,
		      double *surv,
		      double *oldsurv,
		      double *tol,
		      int *maxstep,
		      int *educate,
		      int *niter) {
  
  int i, s, done=0, step=0;
  double atrisk, pl, haz, varhaz, diff, tmpR, tmpL ,survL, survR, lenOBS;
  
    while (done==0 && step < *maxstep){
      /*       Rprintf("Step %d\n",step);  */
      diff=0;
      atrisk = *N;
      pl=1; 
      haz=0;
      varhaz=0;
      nevent[0] = 0;
      ncens[0] = 0;
      
      for (s=0; s < *NS; s++){ /* loop over peto intervals */ 
	nrisk[s]=atrisk;
	for (i=0; i < *N; i++){
	  /* loop only over those intervals */
	  /* that touch the current peto interval */
   	  if (L[i]<=petoR[s] && R[i]>=petoL[s]){   
/* 	    /\* educated first step *\/ */
 	    if (step==0){ 
/* 	      if (*educate==0){ */
		
/* 	      } */
/* 	    else */
	      if (status[i]==0 && L[i] <= petoL[s]) ncens[s]++; /* right censored at L[i] before JL*/
	      if (status[i]==1){ 
		lenOBS = R[i] - L[i];
		if (lenOBS==0 && L[i] == petoL[s]) {
		  nevent[s] ++; /* exact observations */
		}
		if (lenOBS > 0){ /* interval censored  */
		  if (s==0 && L[i]<petoL[s])
		    tmpL=L[i];
		  else if(L[i]>petoL[s])
		    tmpL=L[i];
		  else
		    tmpL=petoL[s];
		  
		  if (s==(*NS-1) && R[i]>petoR[s])
		    tmpR=R[i];
		  else if (R[i]<petoL[s+1])
		    tmpR=R[i];
		  else
		    tmpR=petoL[s+1];
		  nevent[s] += max(0,tmpR - tmpL)/lenOBS;
		}
	      }
	      /* Rprintf("L[i]=%1.2f\tR[i]=%1.2f\tpetoL[s]=%1.2f\tPetoR[s]=%1.2f\tnevent[s]=%1.2f\ttmpL=%1.2f\ttmpR=%1.2f\n",L[i],R[i],petoL[s],petoR[s],nevent[s],tmpL,tmpR);  */
 	    } 
	    else{
	      if (indexL[i]<=1)
		survL=1;
	      else
		survL=surv[indexL[i]-2];
	      if (indexR[i]>=(*NS-1)) 
		survR=0; 
	      survR=surv[indexR[i]-1]; 
	      if (s==0) tmpL=1;
	      else tmpL=surv[s-1];
	      if (s==(*NS-1))
		tmpR=0;
	      else
		tmpR=surv[s];
	      nevent[s] += (tmpL - tmpR)/(survL - survR); 
	      /* Rprintf("i=(%1.0f,%1.0f)\ts=[%1.0f,%1.0f]\tnevent[s]=%1.2f\tsurv[s-1]=%1.2f\tsurv[s]=%1.2f\tsurvL=%1.2f\tsurvR=%1.2f\n",L[i],R[i],petoL[s],petoR[s],nevent[s],tmpL,tmpR,survL,survR);   */
	    }
	  }
	}
	if (nevent[s]>0){
	  haz = nevent[s] / atrisk;
	  pl*=(1 - (nevent[s] / atrisk)); 
	  varhaz += nevent[s] / (atrisk * (atrisk - nevent[s])); 
	}
	if (step>0) oldsurv[s]= surv[s];
	surv[s]=pl;
	/* Rprintf("\ns=%d\tatrisk=%1.8f\tnevent[s]=%1.8f\tsurv[s]=%1.2f\n\n",s,atrisk,nevent[s],surv[s]);   */
	hazard[s] = haz;
	var_hazard[s] = varhaz;
	atrisk-=(nevent[s]+ncens[s]);
	nevent[s+1] = 0;
	ncens[s+1] = 0;
      }
      for (s=0;s<*NS;s++){
	diff=max(max(surv[s]-oldsurv[s],oldsurv[s]-surv[s]),diff);
      }
      if (diff < *tol) done=1;
      step++;
    }
    /*     Rprintf("Step %d\n",step); */
    niter[0]=step; 
}
