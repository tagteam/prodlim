void pl_step(double *pl,double *aj,double *v,double n,double d,int rev);
void prodlim_surv(double *y,double *status,double *time,double *nrisk,double *event,double *loss,double *surv,double *hazard,double *varhazard,int *reverse,int *t,int start,int stop);
void prodlimSurvPlus(double *y,double *status,double *entrytime,double *caseweights,double *time,double *nrisk,double *event,double *loss,double *surv,double *hazard,double *varhazard,int *reverse,int *t,int start,int stop,int *delayed,int *weighted);
void prodlim_clustersurv(double *y,double *status,int *cluster,int *NC,double *time,double *nrisk,double *cluster_nrisk,double *nevent,double *loss,double *ncluster_with_event,double *ncluster_lost,double *sizeof_cluster,double *nevent_in_cluster,double *surv,double *hazard,double *varhazard,double *adj1,double *adj2,double *adjvarhazard,int *t,int start,int stop);
void prodlim_comprisk(double* y,double* status,int* cause,int* NS,double* time,double* nrisk,double* event,double* loss,double* surv,double* cuminc,double* cause_hazard,double* varcuminc,double* cuminc_temp,double* cuminc_lag,double* v1,double* v2,int *t,int start,int stop);
void prodlimCompriskPlus(double* y,double* status,int* cause,double *entrytime,double *caseweights,int* NS,double* time,double* nrisk,double* event,double* loss,double* surv,double* cuminc,double* cause_hazard,double* varcuminc,double* cuminc_temp,double* cuminc_lag,double* v1,double* v2,int *t,int start,int stop,int *delayed,int *weighted);
int neworder (int *a, int *b);
int doubleNewOrder (double *a, double *b);
  
