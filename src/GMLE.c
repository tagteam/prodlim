#include <math.h>
#include <R.h>

#define max(A,B) ((A) > (B) ? (A):(B))
#define min(A,B) ((A) < (B) ? (A):(B))

void GMLE(int *Mstrata,
	  int *Istrata,
	  int *Mindex,
	  int *Iindex,
	  int *N,
	  int *M,
	  double *z,
	  double *oldZ,
	  double *tol,
	  int *maxstep,
	  int *niter){
  
  int i,j,k,l,m,step,done;
  double newZ,nom, denom, diff;
  step=0;
  done=0;
  while (done==0 && step < *maxstep){
    /*     Rprintf("\n\nStep=%d\t\n",step);  */
    diff=0;
    for(k=0;k<*M;k++) oldZ[k]= z[k];
    for(k=0;k<*M;k++){
      nom=0;
      newZ=0;
      for(j=Mstrata[k]; j< Mstrata[k+1];j++){
	i=Mindex[j]-1;
 	denom=0; 
	for(l=Istrata[i]; l < Istrata[i+1];l++){
	  m=Iindex[l]-1;
	  denom += oldZ[m];
	}
	nom = oldZ[k];
	newZ += nom/denom;
      }
      z[k]=newZ/(*N);
    }
    for (k=0;k<*M;k++){
      /*       Rprintf("k=%d\toldZ[k]=%1.2f\tz[k]=%1.2f\tdiff=%1.2f\t\n",k,oldZ[k],z[k],diff);   */
      diff=max(max(z[k]-oldZ[k],oldZ[k]-z[k]),diff);
    }
    if (diff < *tol) done=1;
    step++;
  }
  niter[0]=step; 
}
