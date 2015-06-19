#include <math.h>


/*********************************************************************/
/* declaration of some functions called by 'trans'                    */
/*********************************************************************/
void init_start_risk(int t, int nt, int ns, int u, int* nrisk, int* nstart);

void init_next_risk(int t, int nt, int ns, int* nrisk);

void init_aj(int ns, double* aj);

void set_event(int i, int t, int nt, int ns, 
	       int* tra_from, int* tra_to, int* trow,
	       int* cens_in, int* cpos,
	       int* nevent, int* ncens, int* status, int* nrisk);

void multi_state(int t, int ntr, int ns, int* tra_from, int* tra_to, 		    
		 int* nrisk, int* nevent, double* hazard, 
		 double* aj, double* prob);
 
void compute_hazard(int t, int ntr, int ns, int* tra_from, int* tra_to, 
		    int* nrisk, int* nevent, double* hazard);

void compute_diag(int t, int ns, double* hazard);

void compute_aj(int t, int ns, double* hazard, double* aj);
  
void store_aj(int t, int ns, double* aj, double* prob);



/*********************************************************************/
/* function 'prodlim_multistates' called by C-function 'trans'               */
/*********************************************************************/
void prodlim_multistates(int* n,
			 int* nstates,
			 int* nobserv,
			 int* size,
			 int* ntra,
			 int* tra_from,
			 int* tra_to,
			 int* trow,
			 int* nci,
			 int* cens_in,
			 int* cpos,  
			 double* y,
			 int* status,          
			 int* nstart,
			 double* time,
			 double* hazard,
			 double* prob,
			 int* nevent,
			 int* ncens,
			 int* nrisk,
			 int *first_strata,
			 int *ntimes_strata) {

  int i=0;
  int k=0;
  int s=0;
  int u=0;
  int t=0;


  int nt   = *n;         /* N */
  int ns   = *nstates;   /* number of states, if censoring -1 is included */
  int no   = *nobserv;   /* number of observations */
  int ntr  = *ntra;      /* number of (unique) possible transitions */


  double aj[(ns*ns)];  /* matrix for the aalen-johansen */


  for(i=0; i < no; ++i) { /* loop over the observations (jumps) */

    if( s == 0 ) {      
      /* initialize nrisk with the start distribution for the strata*/
      init_start_risk(t, nt, ns, u, nrisk, nstart);
      
      /* initialize aj */	
      init_aj(ns, aj);
    }
    
    set_event(i, t, nt, ns, tra_from, tra_to, trow,
	      cens_in, cpos, nevent, ncens, status, nrisk);
  
    if( (s < size[u]-1 && y[i] != y[i+1]) || s == size[u]-1 ) {
      /* compute the hazards and aalen */
      multi_state(t, ntr, ns, tra_from, tra_to, nrisk, nevent, hazard, aj, prob);
      
      /* store the time-point */
      time[t] = y[i];
      
      ++t;
      ++k;

      if(s < size[u]-1 ){
	/* initialize nrisk for the next time-point */
	init_next_risk(t, nt, ns, nrisk);
      }
    }
    

    if(s == size[u]-1) {
      first_strata[u]  = t-k+1;
      ntimes_strata[u] = k;
      s=0;
      k=0;
      ++u;
    }
    else {
      ++s;      
    }
  }
}

	   	   
/*********************************************************************/
/* implementation of the functions called by 'trans_multi'           */
/*********************************************************************/
void init_start_risk(int t, int nt, int ns, int u, int* nrisk, int* nstart) {      
  int j = 0;

  nrisk[t*ns + j] = nstart[u];
         
  for(j=1; j < ns; ++j) {	
    nrisk[t*ns + j] = 0;     
  }

  init_next_risk(t, nt, ns, nrisk);
}

void init_next_risk(int t, int nt, int ns, int* nrisk) {
  int j;

  if(t < (nt - 1) ) {	
    for(j=0; j < ns; ++j) {	  
      nrisk[(t+1)*ns + j] = nrisk[t*ns + j];	
    }      
  }  
}

void init_aj(int ns, double* aj) {
  int i,j;

  for(i=0; i < ns; ++i){
    for(j=0; j < ns; ++j) {
      aj[i*ns+j] = 0;
      if( i == j ) {
	aj[i*ns+j] = 1;	  	  
      }
    }
  }
}

void set_event(int i, int t, int nt, int ns,
	       int* tra_from, int* tra_to, int* trow,
	       int* cens_in, int* cpos,
	       int* nevent, int* ncens, int* status, int* nrisk) {
 
     
  if( status[i] == 1 ) {		
    /* add the transition */	
    nevent[ (t*ns*ns) + (tra_from[trow[i]]*ns + tra_to[trow[i]]) ] += 1;
 
    /* risk */
    if(t < (nt - 1) ) {
      nrisk[ (t+1)*ns + tra_from[trow[i]] ] = nrisk[ (t+1)*ns + tra_from[trow[i]] ] - 1;	
      nrisk[ (t+1)*ns + tra_to[trow[i]] ]   = nrisk[ (t+1)*ns + tra_to[trow[i]] ]   + 1;
    } 
  }
  else {

    /* add censoring */
    ncens[ (t*ns) + cens_in[cpos[i]] ] += 1;

    /* risk */
    if(t < (nt - 1) ) {
      nrisk[ (t+1)*ns + cens_in[cpos[i]] ]   = nrisk[ (t+1)*ns + cens_in[cpos[i]] ] - 1;
    }    
  }
}
          
void multi_state(int t, int ntr, int ns, int* tra_from, int* tra_to, 		    
		 int* nrisk, int* nevent, double* hazard, 
		 double* aj, double* prob) {
  
  /* compute the hazards */		  
  compute_hazard(t, ntr, ns, tra_from, tra_to, nrisk, nevent, hazard);

  /* compute the aalen-johansen */	
  compute_aj(t, ns, hazard, aj);

  /* store the aalen-johansen for time-point t */
  store_aj(t, ns, aj, prob);
}

void compute_hazard(int t, int ntr, int ns, int* tra_from, int* tra_to, 
		    int* nrisk, int* nevent, double* hazard) {
  int j;
      
  /* compute the hazards */    
  for(j=0; j < ntr; ++j) {	
    if(nevent[(t*ns*ns) + (tra_from[j]*ns + tra_to[j])] > 0 ) {

	
      hazard[(t*ns*ns) + (tra_from[j]*ns + tra_to[j])] = 
	(double) nevent[(t*ns*ns) + (tra_from[j]*ns + tra_to[j])] / nrisk[t*ns + tra_from[j]];	        
    }	    
  }     

  /* compute the diagonal of the matrix hazard[(t*ns*ns)] */
  compute_diag(t, ns, hazard);
}


void compute_diag(int t, int ns, double* hazard) { 
  int r,c;
  double sumrow;
  
  /* compute the diagonal elements: the sum over each row must be 1 */
  for(r=0; r < ns; ++r ) {              
    sumrow = 0.;       
      
    for( c = 0; c < ns; ++c ) {	  
      if( c != r ) {	    
	sumrow += hazard[(t*ns*ns) + (r*ns+c)];	  
      }	         
    }
	
    hazard[(t*ns*ns)+ (r*ns+r)] = (double)(1 - sumrow);           
  } 
}     

void compute_aj(int t, int ns, double* hazard, double* aj) {    
  int r,c,i;

  double m[ns*ns];

  for(r=0; r < ns; ++r) {		
    for(c=0; c < ns; ++c) {
      m[r*ns+c] = 0.0;
      for(i=0; i < ns; ++i) {
	m[r*ns+c] += aj[r*ns+i] * hazard[(t*ns*ns) + (i*ns+c)];
      }
    }	
  }

  for(i=0; i < (ns*ns); ++i) {
    aj[i] = m[i];
  }
}
   
void store_aj(int t, int ns, double* aj, double* prob) {
  int i;

  for(i=0; i < (ns*ns); ++i) {
    prob[(t*ns*ns) + i] = aj[i];
  }
}
