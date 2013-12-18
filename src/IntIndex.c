#include <math.h>
#include <R.h>
void IntIndex(double *L,
	    double *R,
	    double *p,
	    double *q,
	    int *N,
	    int *M,
	    int *Iindex,
	    int *Mindex,
	    int *Istrata,
	    int *Mstrata){
  
  int i,m,k,l;
  k=0;
  
  for (i=0; i<*N;i++){
    for (m=0; m<*M;m++){
      if ((L[i]==R[i] && p[m]==q[m] && L[i]==q[m])  /* point */
	  ||
	  (L[i]<q[m] && L[i]<=p[m] && R[i]>=q[m] && R[i]>p[m]))  /* interval */
	{
	  Iindex[k]=m+1;
	  k++;
	}
    }
    Istrata[i]=k;
  }
  l=0;
  for (m=0; m<*M;m++){
    for (i=0; i<*N;i++){
      if ((L[i]==R[i] && p[m]==q[m] && L[i]==q[m])  /* point */
	  ||
	  (L[i]<q[m] && L[i]<=p[m] && R[i]>=q[m] && R[i]>p[m]))  /* interval */
	{
	  Mindex[l]=i+1;
	  l++;
	}
    }
    Mstrata[m]=l;
  }
}




void Turnb(int *Mstrata,
	   int *Istrata,
	   int *Mindex,
	   int *Iindex,
	   int *N,
	   int *M,
	   double *Z,
	   double *nplme){
  
  int i,l,u,j,Iind, Mind;
  double Ilast, ZI, ZM, Mlast, Zlast, ZMI;
  
  for(i=0;i<*M;i++){

  Zlast=0;
  ZMI=0;

  for(l=0;l<*N; l++){
    
    Mlast=0;
    ZM=0;
    Mind=0;
    
    for(u=Mstrata[l];u<Mstrata[l+1];u++){
      
      Mind=Mindex[u];
      
      Ilast=0;
      ZI=0;
      Iind=0;
  
      for(j=Istrata[l]; j<Istrata[l+1];j++)
    { 
      Iind=Iindex[j];
      ZI=Z[Iind-1]+Ilast;
      Ilast=ZI;
    }
  
  ZM=Z[Mind-1]/Ilast + Mlast;
  Mlast=ZM;
 
 }

  ZMI=Mlast+Zlast;
  Zlast=ZMI;
  }

  nplme[i]=Mlast;
  }
}
