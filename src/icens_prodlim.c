/*

  The product limit method for interval censored data

  Copyright 2007-2009 Department of Biostatistics, University of Copenhagen

  Written by Thomas Alexander Gerds <tag@biostat.ku.dk>
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public
  License along with this program; if not, write to the Free
  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA  02110-1301, USA.


  The structure of the algorithm:
  
  looping until convergence or maxstep
  over all grid points
  starting with the interval
  [grid[0] ; grid[1]]
      
  the first time s=0 is a dummy time
  used to catch exact events at 0.
  to compute the hazard and the survival
  probability at the END of a grid interval
      
  [grid[s] ; grid[s+1]]
  
  first count events and censored between
  grid[s] and grid[s+1], then devide
  by the number at risk at grid[s].
  Note: nevent[s+1] is the number of
  subjects at risk at time grid[s].

  use only the observed intervals

  [L[i],R[i]]

  that overlap the currentn
  grid interval:

  [grid[s] ; grid[s+1]]
  
  whether or not an interval overlaps is determined by
  iindex, a vector of indices where the part
  from imax[x] to imax[x+1] identifies observations that
  overlap grid interval x

  Exact and right censored observations are handled
  as for the usual Kaplan-Meier method.
  Real interval censored observations contribute to
  the number of events by the relative to the overlap
  with the current grid-interval.

  To compute the relative event count at the very first step
  assume a uniform distribution, in subsequent steps
  use the survival probability of in the previous step

*/


#include <math.h>
#include <R.h>

#define max(A,B) ((A) > (B) ? (A):(B))
#define min(A,B) ((A) < (B) ? (A):(B))

void icens_prodlim(double *L,
		   double *R,
		   double *grid,
		   int *indexL,
		   int *indexR,
		   int *iindex,
		   int *imax,
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
		   int *niter) {
  
  int i, j, s, done=0, step=0, n, ns, start, stop;
  /*   int verbose; */
  double atrisk, pl, haz, varhaz, diff, survL, survR, lenOBS, nom;
  
  n = (int) *N; /* number of interval censored observations */
  ns = (int) *NS;  /* number of grid points + 1 */

  while (done==0 && step < *maxstep){
    
    surv[0]=1;
    oldsurv[0]=1; 
    diff=0;
    atrisk = *N;
    nrisk[0]= *N;
    varhaz=0;
    haz=0;
    pl=1;
    start=0;
    stop=max(0,imax[0]);

    /* LOOP OVER GRID INTERVALS */
    for (s=0; s < (ns-2); s++){
      nrisk[s+1]=atrisk;
      nevent[s+1] = 0;
      ncens[s+1] = 0;
      
      /* LOOP OVER OBSERVED INTERVALS */
      for (j=start; j < stop; j++){
	i=iindex[j]-1; /* R starts counting at 1 */
	if (status[i]==0 && L[i] == grid[s+1]) ncens[s+1]++; /* right censored */
	if (status[i]>0){
	  lenOBS = R[i] - L[i];
	  if (lenOBS==0 && L[i] == grid[s+1]) nevent[s+1] ++; /* exact observation */
	  if (lenOBS > 0){
	    if (L[i] < grid[s+1] && R[i]>grid[s]){
	      if (step==0){ 
		nevent[s+1] += max(0,min(R[i],grid[s+1]) -max(grid[s],L[i]))/lenOBS;
	      }
	      else{
		survL = surv[indexL[i]-1]; 
		survR = surv[indexR[i]-1];
		nom = (min(survL,surv[s]) - max(surv[s+1],survR)); /* overlap */
		if (nom>=*tol) nevent[s+1] += nom/(survL-survR); 
	      }
	    }
	  }
	}
      }
      start=max(0,imax[s]);
      stop=max(imax[s+1],0);
      if (nevent[s+1]>0){
	haz = nevent[s+1] / (double) atrisk;
	pl*=(1 - (nevent[s+1] / (double) atrisk)); 
	varhaz += nevent[s+1] / (double) (atrisk * (atrisk - nevent[s+1])); 
      }
      if (step>0) oldsurv[s+1]= surv[s+1]; /* move the current estimate to oldsurv */
      surv[s+1]=pl; /* update the survival probability */
      hazard[s+1] = haz;
      var_hazard[s+1] = varhaz;
      atrisk-=(nevent[s+1]+ncens[s+1]);  /* update the number at risk */
    }
    for (s=0;s<(ns-2);s++){  /*  check if the algorithm converged */
      diff=max(max(surv[s+1]-oldsurv[s+1],oldsurv[s+1]-surv[s+1]),diff);
    }
    if (diff < *tol) done=1;
    step++;
  }
  niter[0]=step; 
}



/*   verbose=-2; */
/* THE CURRENT SURVIVAL ESTIMATE
  if (verbose>=0){ 
  Rprintf("\nStep %d\n",step); 
  for (s=0; s < (ns-2); s++)
  Rprintf("s(%1.2f)=%1.2f\n",grid[s],surv[s]); 
  Rprintf("\n\n",step); 
  } 
*/

/* THE GRID INTERVAL
  if (step<=verbose){  
  Rprintf("\n");  
  Rprintf("grid=[%1.3f,%1.3f]\n",grid[s],grid[s+1]);}
*/


/*
   THE OBSERVED INTERVAL
  if (step<=verbose){  
  Rprintf("\n");  
  Rprintf("Obs=[%1.3f,%1.3f]\n",L[i],R[i]);  
  }
*/

/* THE EVENT COUNT IN STEPS >0  
   if (step<=verbose){  
   Rprintf("survGrid=[%1.2f,%1.2f]\tsurvObs=[%1.2f,%1.2f]\tzaehl=%1.2f\tnenn=%1.2f\tjump=%1.2f\n",surv[s],surv[s+1],survL,survR,nom,(survL-survR),nevent[s+1]);  
   }
*/

      
/* EVENTS, ATRISK, SURVPROB
  if (step<=verbose){ 
  Rprintf("nevent=%1.2f\tnrisk=%1.2f\tsurv=%1.2f\t\n",nevent[s+1],atrisk,pl); 
  }
*/

