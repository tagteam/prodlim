/*

  define symmetric neighborhoods for unique values u in x

  input
  =====
  n: the sample size
  nu: number of unique x values 
  cumtabu: n times the cumulative empirical df at u
  cumtabx: n times the cumulative empirical df at x
  tabx: frequency of x
  radius: n times the bandwidth

  output specific to neighborhood's 
  =================================
  first: the first neighbor
  size:  the size the neighborhood
  neighbors sorted from the first to last neighborhood
  
*/

#include <stdlib.h> 
void neighborhood(int *first,
		  int *size,
		  int *cumtabu,
		  int *cumtabx,
		  int *tabx,
		  int *radius,
		  int *nu,
		  int *n){
  int u,last;
  
  for (u=0;u<*nu;u++){
    
    /* make a first guess */
    
    first[u]=cumtabu[u]-*radius;
    last=cumtabu[u]+*radius;
    
    /* if x[first[u]] is tied, move
       to the first[u] member in the bin */
    
    if (first[u]<=0) first[u]=1;
    else first[u] = cumtabx[first[u]-1]-tabx[first[u]-1]+1;
    
    /* if x[last] is tied and not the last
       in its bin, move to the previous bin */
    
    if (last>*n) last=*n;
    else if (cumtabx[last-1] > last) last=cumtabx[last-1]-tabx[last-1];
    
    size[u]=last-first[u]+1;
  }
}

int neworder (int *a, int *b){
  if (*a < *b) return -1; else return 1;}

void neighbors(int *first,
	       int *size,
	       int *orderx,
	       int *neighbors,
	       int *nu){
  int u,i,new,start=0;

  /* fill the neighborhoods */
  new=0;
  for (u=0;u<*nu;u++){
    for (i=0;i<size[u];i++){
      neighbors[new] = orderx[first[u]-1+i];
      new++;
    }

    /*     sort the neighborhood */
    qsort(neighbors+start,size[u],(size_t) sizeof(int),(int (*)(const void *, const void *))(neworder));  

    start+=size[u];
  }
}
  


