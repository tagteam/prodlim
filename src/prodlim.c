/*
  (2006--2010) Thomas A. Gerds 
  --------------------------------------------------------------------
  distributed under the terms of the GNU public license 
  y the SORTED failure times with ties 
  status is 1 if the individual has failed (from any cause), 0 otherwise 
  cause indicates the cause 
  N is the length of Y  
  NC is the number of different clusters 
  NS is the number of states (aka causes) 
  cluster indicates the cluster
  size is a vector with the number of individuals in strata
*/
	      
#include <math.h>
#include <R.h>
#include "prodlim.h"

void prodlim(double *y,
	     int *status,
	     int *cause,
	     int *cluster,
	     int *N,
	     int *NS,
	     int *NC,
	     int *NU,
	     int *size,
	     double *time,
	     double *nrisk,
	     int *event,
	     int *lost,
	     double *surv,
	     double *cuminc,
	     double *hazard,
	     double *varhazard,
	     double *extra_double,
	     int *extra_int,
	     int *len,
	     int *size_strata,
	     int *first_strata,
	     int *reverse,
	     int *model,
	     int *independent) {
  int t, u, start, stop, size_temp;
  t=0;
  start=0;
  size_temp=0;
  for (u=0;u<*NU;u++){
    stop=start+size[u];
    if (*model==0)
      if (*independent==1) 
	/* shake(int start, int stop, int back, int m, int **res ) */
	  prodlim_surv(y,status,time,nrisk,event,lost,surv,hazard,varhazard,reverse,&t,start,stop);
      else{
	double *cluster_nrisk, *adj1, *adj2, *adjvarhazard;
	int *ncluster_lost, *ncluster_with_event, *sizeof_cluster, *nevent_in_cluster, *max_nc;
	max_nc = extra_int;
	/*
	  tag: 12 Nov 2010 (18:41)
	  
	the length of nrisk, nevent and lost is 2 * N
	the first half is used for the individual level
	the second for the cluster level.

	the function is thus still restricted to a single cluster variable
	*/
	  cluster_nrisk = nrisk + *N;
	ncluster_with_event = event + *N;
	ncluster_lost = lost + *N;
	adjvarhazard = varhazard + *N;
	adj1 = extra_double;
	adj2 = extra_double + *max_nc;
	nevent_in_cluster = extra_int + 1;
	sizeof_cluster = 1 + extra_int + *max_nc;
	prodlim_clustersurv(y,status,cluster,NC + u,time,nrisk,cluster_nrisk,event,lost,ncluster_with_event,ncluster_lost,sizeof_cluster,nevent_in_cluster,surv,hazard,varhazard,adj1,adj2,adjvarhazard,&t,start,stop);}
    else
      if (*model==1){
	double *cuminc_temp, *cuminc_lag, *v1, *v2;
	cuminc_temp = extra_double;
	cuminc_lag = extra_double + *NS;
	v1 = extra_double + *NS + *NS;
	v2 = extra_double + *NS + *NS + *NS;
	prodlim_comprisk(y,status,cause,NS,time,nrisk,event,lost,surv,cuminc,hazard,varhazard,cuminc_temp,cuminc_lag,v1,v2,&t,start,stop);
      }
      /* else */
	/* Rprintf("Unknown model"); */
    start+=size[u];
    size_strata[u] = t - size_temp;
    first_strata[u] = t + 1 - size_strata[u];
    size_temp += size_strata[u];
  }
  *len=t;
}


void pl_step(double *pl,double *aj,double *v,int n,int d,int rev)
{
  if (d > 0){	    
    *aj = (d / (double) (n - rev));	/* nelson-aalen */
    *pl *= (1 - *aj); /* product limit */
    *v += (double) d / ((double) (n - rev) * (double) (n - rev - d)); /* greenwood variance */
  }
}











