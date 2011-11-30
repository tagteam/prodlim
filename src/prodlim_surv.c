#include <math.h>
#include <R.h>
#include "prodlim.h"


void prodlim_surv(double *y,
		 int *status,
		 double *time,
		 double *nrisk,
		 int *event,
		 int *loss,
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

    if (y[i]==y[i-1] && i<stop){
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
