#include <R.h>
void life_table(int *pred_nrisk,
		int *pred_nevent,
		int *pred_nlost,
		int *nrisk,
		int *nevent,
		int *nlost,
		double *lower,
		double *upper,
		double *eventTime,
		int *first,
		int *size,
		int *NR,
		int *NT,
		int *delayed){

  
  int i,t,s,count_e,count_l,First,Last;
  double min_eventTime, max_eventTime;

  /*
    Aim:

     life table intervals are given by
    
    [lower[t] ; upper[t])

    NOTE: the intervals are closed on the right and open on the left

    in a loop across covariate strata    
    find the
    
    a) the number at risk just before lower[t]
    b) the number of uncensored events in interval
    c) the number of censored in interval
    
    Notation:
    
    i: runs through covariate strata
    t: runs through lower and upper
    s: runs through intervals between eventTimes

    the covariate stratum starts at

    First=first[i]-1
    
    and stops at
    
    Last=first[i]-1 + size[i]-1
    
    the censored event times are in `eventTime'

    There are three cases:
    
    (1) the interval lays before the first event time
    (2) the interval lays includes one event time
    (3) the interval lays behind the last event time
    
  */
  for (i=0;i<*NR;i++){
    First=first[i]-1;
    Last=first[i]-1 + size[i]-1;
    min_eventTime = eventTime[First];
    max_eventTime = eventTime[Last];
    s=0;
    for (t=0;t<(*NT);t++){
      count_e =0;
      count_l =0;
      if (upper[t] < min_eventTime){
	/*
	  case (1) interval before the first event time:
	  
	  [)....

	  with delayed entry no one is at risk before the first entry time
	*/
	if (*delayed)
	  pred_nrisk[t + i *(*NT)] = 0;
	else
	  pred_nrisk[t + i *(*NT)] = nrisk[First];
	pred_nevent[t + i *(*NT)] = 0;
	pred_nlost[t + i *(*NT)] = 0;
      }
      else{
	if (lower[t] > max_eventTime){ /* the left side of the interval is larger than max_eventTime.*/
	  /*
	    case (3) after the last eventTime: ....[)
	  */
	  while(t<(*NT)){
	    pred_nrisk[t + i *(*NT)] = 0;
	    pred_nevent[t + i *(*NT)] = 0;
	    pred_nlost[t + i *(*NT)] = 0;
	    t++; 
	  } 
	}
	else{
	  /*
	    case (2) between .[..)..
	    here
	    upper[t] >= min_eventTime
	    and
	    lower[t] <= max_eventTime
	  */

	  /*
	    first find number at risk just before lower[t] ...
	  */
	  /* Rprintf("s=%d\tt=%d\ti=%d\tupper[t]=%1.2f\tlower[t]=%1.2f\teventTime[i]=%1.2f\t\n",s,t,i,upper[t],lower[t],eventTime[i]); */
	  if (s==0){
	    if (*delayed)
	      	    pred_nrisk[t + i *(*NT)] = 0;
	    else
         	    pred_nrisk[t + i *(*NT)] = nrisk[First];
	  }
	  else{
	    pred_nrisk[t + i *(*NT)] = nrisk[First+s];
	  }
	  /* ... then count events and censored in interval [lower[t],upper[t]) */
	  
	  while ((s <= size[i]-1) && (eventTime[First + s] < upper[t])){
	    count_e +=nevent[First+s];
	    count_l +=nlost[First+s];
	    s++;
	  }
	  pred_nevent[t + i *(*NT)] = count_e;
	  pred_nlost[t + i *(*NT)] = count_l;
	  /*
	    now s is such that either
	    eventTime[First + s] >= upper[t] =lower[t+1]
	    or
	    s==size[i]
	  */
	}
      }
    }
    /* do NOT reset s because the
       next event Time is greater
       or equal to the current.
    */
  }
}
