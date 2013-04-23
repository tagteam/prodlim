#include <math.h>
#include <R.h>
#include "prodlim.h"

void prodlim_clustersurv(double *y,
			 int *status,
			 int *cluster,
			 int *NC,
			 double *time,
			 double *nrisk,
			 double *cluster_nrisk,
			 int *nevent,
			 int *lost,
			 int *ncluster_with_event,
			 int *ncluster_lost,
			 int *sizeof_cluster,
			 int *nevent_in_cluster,
			 double *surv,
			 double *hazard,
			 double *varhazard,
			 double *adj1,
			 double *adj2,
			 double *adjvarhazard,
			 int *t,
			 int start,
			 int stop){
  
  int s,i,l,k;
  double surv_step, hazard_step, V1, V2, atrisk, cluster_atrisk;
  /*   Rprintf("Call clustersurv\n\n");  */
  /* initialize the time counter */
  s = (*t); 
    
  /*
    cluster is an indicator of the cluster number.
    for example if the individual (tooth) 'i'
    belongs to patient 'k' then 'cluster[i]=k'
    First we need to re-initialize sizeof_cluster, nevent_in_cluster, etc
    are set to zero.    
  */
  for (k=0;k<*NC;k++) {
    sizeof_cluster[k]=0;
    nevent_in_cluster[k]=0;
    adj1[k]=0;
    adj2[k]=0;
  }
  /*
    Then, the vector "sizeof_cluster" is
    initialized with the current number of individuals 
    in the cluster.
  */

  for (i=start;i<stop;i++) sizeof_cluster[cluster[i]-1]++;
  
  /* initialize  */
  surv_step=1; hazard_step=0; V1=0; V2=0;
  atrisk=(double) stop-start;
  cluster_atrisk= (double) *NC;
  nevent[s] = status[start];
  ncluster_with_event[s] = status[start];
  ncluster_lost[s] = 0;
  nevent_in_cluster[cluster[start]-1] = status[start];
  lost[s] = (1-status[start]);
  
  for(i=(1+start);i <=stop;i++){
    /*
      start at i=1 to deal with ties.
      first check if current time is equal
      to the previous time.
    */
    if (i<stop && y[i]==y[i-1]){
      nevent[s] += status[i];
      lost[s] += (1 - status[i]);
      nevent_in_cluster[cluster[i]-1] += status[i];
      if (cluster[i]!=cluster[i-1]){
	ncluster_with_event[s]+= status[i];
      }
    }
    else {
      time[s] = y[i-1];
      nrisk[s] = atrisk;
      cluster_nrisk[s] = cluster_atrisk;

      /* marginal Kaplan-Meier and naive variance estimator */
      pl_step(&surv_step, &hazard_step, &V1, atrisk, nevent[s],0);
	
      surv[s]=surv_step;
      hazard[s]=hazard_step;
      varhazard[s] = V1;
	
      /* adjusted variance estimator of Ying and Wei (1994)  */
	
      V2=0;
      for (k=0;k<*NC;k++) {
	adj1[k] += nevent_in_cluster[k] / (double) atrisk;
	adj2[k] += sizeof_cluster[k] * nevent[s] / (double) (atrisk * atrisk);	
	V2 += (adj1[k]-adj2[k]) * (adj1[k]-adj2[k]);
      }	
      /* collect the results for unique time points */
      surv[s] = surv_step; 
      varhazard[s]=V1;
      adjvarhazard[s]=V2; 
      
      /* initialize the next time point */
      if (i < stop) {
	atrisk-=(nevent[s]+lost[s]);
	/*
	  looking back at the individuals with tie at time s
	  this makes sense as presently: y[i]!=y[i-1]
	*/
	for (l=1;l<=(nevent[s]+lost[s]);l++) {
	  /*
	    decrease the size of corresponding clusters
	  */
	  sizeof_cluster[cluster[i-l]-1]--;
	  /*
	    if the last obs in a cluster is gone, then
	    the number of clusters at risk is decreased
	    by 1.
	  */	  
	  if (sizeof_cluster[(cluster[i-l]-1)]==0) {
	    cluster_atrisk--;
	    ncluster_lost[s] += (1-status[i-l]);
	  }
	  nevent_in_cluster[cluster[i-l]-1]=0; /* reset for next time point  */
	}
	s++;
	nevent_in_cluster[cluster[i]-1] = status[i];
	nevent[s] = status[i];
	ncluster_with_event[s] = status[i];
	lost[s] = (1-status[i]);
      }
    }
  }
  *t=(s+1); /* for the next stratum and finally for R */
}

