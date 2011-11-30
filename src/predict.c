#include <R.h>

void findex(int *findex,
	    int *type,
	    int *S,
	    int *freq_strata,
	    double *Z,
	    double *NN,
	    int *NR,
	    int *NT){
  
  int i,x,last;
  
  for (i=0;i<*NR;i++){
    
    /* goto strata of subject i */
    if (S[i]==1)
      x=0;
    else
      x = freq_strata[S[i]-2];
    last = freq_strata[S[i]-1] -1;
    
    /* find the closest neighbor */
    if (*type==0)
      findex[i]=x;
    else{
      if (Z[i] <= NN[x]) /* <= first */
	findex[i] = x; 
      else{
	if (Z[i] >= NN[last]){/* >= last */
	  findex[i] = last;
	}
	else { /* sitting between two neighbors*/
	  while (Z[i] >= NN[x]) x++;
	  if ((NN[x] - Z[i]) < (Z[i] - NN[x-1]))
	    findex[i] = x;
	  else
	    findex[i] = x-1;
	}
      }
    }
    findex[i]+=1; /* in `R' counting starts at 1 */
  }
}

void pred_index(int *pindex,
		double *Y,
		double *time,
		int *first,
		int *size,
		int *NR,
		int *NT){
  
  int i,t,f;
    
  for (i=0;i<*NR;i++){    
    f=0;
    for (t=0;t<(*NT);t++){
      
      if (Y[t] < time[first[i]-1]){ /* < first */
	pindex[t + i * (*NT)] = 0;
      }
      else{
	if (Y[t] > time[first[i]-1 + size[i]-1]){ /* > last */
	  while(t<(*NT)){ 
	    pindex[t + i * (*NT)] = -1;
	    t++; 
	  } 
	}
	else{ /* sitting between to jump times */
	  
	  while (Y[t] >= time[first[i]-1 + f]
		 && f <= size[i]-1)
	    f++;
	  pindex[t + i * (*NT)] = first[i] -1 + f;
	  /* do NOT reset f because the next requested time
	     is greater or equal to the current time */
	}
      }
    }
  }
}




