/* compute the values of a step function, 
   ie how many of the jumps are smaller or
   equal to the eval points  */

void sindex(int *index,
	    double *jump,
	    double *eval,
	    int *N,
	    int *NT,
	    int *strict){
  int i,t;
  index[0] = 0;
  i = 0;
  if (*strict==0){
    for (t=0;t<*NT;t++){
      while(i<*N && jump[i]<=eval[t]) i++;
      index[t] = i;
    }
  }
  else{
    for (t=0;t<*NT;t++){
      while(i<*N && jump[i] < eval[t]) i++;
      index[t] = i;
    }
  }
}
