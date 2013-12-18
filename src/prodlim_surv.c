#include <math.h>
#include <R.h>
#include "prodlim.h"


void prodlim_surv(double *y,
		  double *status,
		  double *time,
		  double *nrisk,
		  double *event,
		  double *loss,
		  double *surv,
		  double *hazard,
		  double *varhazard,
		  int *reverse,
		  int *t,
		  int start,
		  int stop
		  ){
  int i,s;
  double surv_temp,hazard_temp,varhazard_temp,atrisk;
  
  s=(*t);
  surv_temp=1; 
  hazard_temp=0;
  varhazard_temp=0;
  atrisk=(double) stop-start;
  
  event[s] = status[start];
  loss[s] = (1-status[start]);
  
  for (i=(1+start);i<=stop;i++){

    if (i<stop && y[i]==y[i-1]){
      event[s] += status[i];
      loss[s]  += (1-status[i]);
    }
    else {
      time[s]=y[i-1];
      nrisk[s]=atrisk;
      
      if (*reverse==1)
	pl_step(&surv_temp, &hazard_temp, &varhazard_temp, atrisk, loss[s], event[s]);
      else
	pl_step(&surv_temp, &hazard_temp, &varhazard_temp, atrisk, event[s], 0);
      surv[s]=surv_temp;
      hazard[s]=hazard_temp;
      varhazard[s] = varhazard_temp;
      
      if (i<stop){
	atrisk-=(event[s]+loss[s]);
	s++;
	event[s]=status[i];
	loss[s]=(1-status[i]);
      }
    }
  }
  *t=(s+1);			/* for the next strata and finally for R */
}


void prodlim_surv_weighted(double *y,
			   double *status,
			   double *caseweights,
			   double *time,
			   double *nrisk,
			   double *event,
			   double *loss,
			   double *surv,
			   double *hazard,
			   double *varhazard,
			   int *reverse,
			   int *t,
			   int start,
			   int stop
			   ){
  int i,s;
  double surv_temp,hazard_temp,varhazard_temp,atrisk;
  
  s=(*t);
  surv_temp=1; 
  hazard_temp=0;
  varhazard_temp=0;
  atrisk=0;
  for (i=start;i<stop;i++) atrisk += caseweights[i];
  i=0;
  /* atrisk=(double) stop-start; */
  event[s] = caseweights[start] * status[start];
  loss[s] = caseweights[start] * (1-status[start]);
  
  for (i=(1+start);i<=stop;i++){

    if (i<stop && y[i]==y[i-1]){
      event[s] += caseweights[i] * status[i];
      loss[s]  += caseweights[i] * (1-status[i]);
    }
    else {
      time[s]=y[i-1];
      nrisk[s]=atrisk;
      if (*reverse==1)
	pl_step(&surv_temp, &hazard_temp, &varhazard_temp, atrisk, loss[s], event[s]);
      else
	pl_step(&surv_temp, &hazard_temp, &varhazard_temp, atrisk, event[s], 0);
      surv[s]=surv_temp;
      hazard[s]=hazard_temp;
      varhazard[s] = varhazard_temp;
      
      if (i<stop){
	atrisk-=(event[s]+loss[s]);
	s++;
	event[s]=caseweights[i] * status[i];
	loss[s]=caseweights[i] * (1-status[i]);
	/* Rprintf("\ntime=%1.2f\tnrisk=%1.2f\tnev=%1.2f\tstatus=%1.2f\tcasew=%1.2f\n",y[i],atrisk,event[s],status[i],caseweights[i]); */
      }
    }
  }
  *t=(s+1);			/* for the next strata and finally for R */
}
