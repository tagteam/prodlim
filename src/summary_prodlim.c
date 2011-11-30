#include <R.h>
void summary_prodlim(int *pred_nrisk,
		     int *pred_nevent,
		     int *pred_nlost,
		     int *nrisk,
		     int *nevent,
		     int *nlost,
		     double *evalTime,
		     double *eventTime,
		     int *first,
		     int *size,
		     int *NR,
		     int *NT){
  
  int i,t,s,First,Last;
  double min_eventTime, max_eventTime;

  /*
    in a loop across covariate strata, count events,
    right censored (lost) and numbers at risk
    at the eval time points:

    we aim to find the

    a) number at risk just before evalTime[t]
    b) the number of uncensored events at evalTime[t]
    c) the number of censored at evalTime[t] 
    
    i: covariate strata
    t: runs through evalTime
    s: runs through intervals between eventTimes
    
    the requested time points are in `evalTime'
    the censored event times are in `eventTime'

    There are three cases:
    
    (1) before the first event time
    (2) between event times
    (3) after the last event time
    
    the covariate stratum starts at

    First=first[i]-1
    
    and stops at
    
    Last=first[i]-1 + size[i]-1
  */
  
  for (i=0;i<*NR;i++){
    First=first[i]-1;
    Last=first[i]-1 + size[i]-1;
    min_eventTime = eventTime[First];
    max_eventTime = eventTime[Last];
    s=0;
    for (t=0;t<(*NT);t++){
      if (evalTime[t] < min_eventTime){
	pred_nrisk[t + i *(*NT)] = nrisk[First];
	pred_nevent[t + i *(*NT)] = 0;
	pred_nlost[t + i *(*NT)] = 0;
      }
      else{
	if (evalTime[t] > max_eventTime){
	  while(t<(*NT)){
	    pred_nrisk[t + i *(*NT)] = 0;
	    pred_nevent[t + i *(*NT)] = 0;
	    pred_nlost[t + i *(*NT)] = 0;
	    t++; 
	  }
	}
	else{
	  /* move to the largest event time before the eval time  */
	  while ((eventTime[First + s] < evalTime[t]) && (s <= size[i]-1)){
	    s++;
	  }
	  /* Rprintf("s=%d\tevalTime=%1.2f\teventTime[First+s]=%1.2f\tFirst=%d\tnrisk=%d\n",s,evalTime[t],eventTime[First+s],First,nrisk[First+s]);  */
	  pred_nrisk[t + i *(*NT)] = nrisk[First+s];
	  if  (eventTime[First + s] == evalTime[t]){
	    pred_nevent[t + i *(*NT)] = nevent[First+s];
	    pred_nlost[t + i *(*NT)] = nlost[First+s];
	  }
	  else{
	    pred_nevent[t + i *(*NT)] = 0;
	    pred_nlost[t + i *(*NT)] = 0;
	  }
	}
      }
      /* do NOT reset s because the
	 next evalTime is greater
	 or equal to the current.
      */
    }
  }
}
