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
  *t=(s+1); /* for the next strata and finally for R */
}

int doubleNewOrder (double *a, double *b){
  if (*a < *b) return -1; else return 1;}

void prodlimSurvPlus(double *y,
		     double *status,
		     double *entrytime,
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
		     int stop,
		     int *delayed,
		     int *weighted
		     ){
  int i,e,s,entered;
  double surv_temp,hazard_temp,varhazard_temp,atrisk;
  
  e=0;
  s=(*t);
  surv_temp=1; 
  hazard_temp=0;
  varhazard_temp=0;
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
        e=start; /* index for delayed entry */
    }else{
      atrisk=(double) stop-start;
    }
  }
  
  s=(*t);
  if (*weighted==1){
    event[s] = caseweights[start] * status[start];
    loss[s] = caseweights[start] * (1-status[start]);
  }else{  
    event[s] = status[start];
    loss[s] = (1-status[start]);
  }
  
  for (i=(1+start);i<=stop;i++){
    if (i<stop && y[i]==y[i-1]){ /* for ties  */
      if (*weighted==1){
	event[s] += caseweights[i] * status[i];
	loss[s]  += caseweights[i] * (1-status[i]);
      }else{
	event[s] += status[i];
	loss[s]  += (1-status[i]);
      }
    }
    else {
      if (*delayed==1){
	/* delayed entry: find number of subjects that
	   entered at time[s] */
	entered=0;
	while(e<stop && entrytime[e] < y[i-1]){ /*entry happens at t+ events at t*/
	  /* unless there is a tie between the current
	     and the next entry-time, add time to list of times, increase s
	     and move the values of event, loss etc. to the next event time */
	  entered++;
	  if (e==start || entrytime[e]>entrytime[e-1]){
	    nrisk[s]=atrisk+entered;
	    if (entrytime[e]!=time[s-1]){ /* if entrytime[e]==y[i] then only increase
					the number at risk but not move the
					time counter or the values of event, etc.
				      */
	      /* Rprintf("e=%d\ts=%d\tentrytime[e]=%1.2f\ty[i-1]=%1.2f\ttime[s]=%1.2f\ti=%d\t\n",e,s,entrytime[e],y[i-1],time[s],i); */
	      event[s+1]=event[s];
	      event[s]=0;
	      loss[s+1]=loss[s];
	      loss[s]=0; 
	      surv[s]=surv_temp; 
	      hazard[s]=0; 
	      varhazard[s]=varhazard_temp;
	      time[s]=entrytime[e];
	      s++;
	    } 
	  }
	  e++; /* increase cumulative counter  */
	}
	atrisk += (double) entered;
      }
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
	if (*weighted==1){
	  event[s] = caseweights[i] * status[i];
	  loss[s] = caseweights[i] * (1-status[i]);
	}else{  
	  event[s] = status[i];
	  loss[s] = (1-status[i]);
	}
	/* Rprintf("e=%d\tstop=%d\ti=%d\ts=%d\ttime=%1.2f\tstatus=%1.2f\tevent=%1.2f\t\n",e,stop,i,s,time[s],status[i],event[s]);  */
      }
    }
  }
  *t=(s+1); /* for the next strata and finally for R */
}
