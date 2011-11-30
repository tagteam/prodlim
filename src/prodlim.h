void pl_step(double *pl,double *aj,double *v,int n,int d,int rev);
void prodlim_surv(double *y,int *status,double *time,double *nrisk,int *event,int *loss,double *surv,double *hazard,double *varhazard,int *reverse,int *t,int start,int stop);
void prodlim_clustersurv(double *y,int *status,int *cluster,int *NC,double *time,double *nrisk,double *cluster_nrisk,int *nevent,int *loss,int *ncluster_with_event,int *ncluster_lost,int *sizeof_cluster,int *nevent_in_cluster,double *surv,double *hazard,double *varhazard,double *adj1,double *adj2,double *adjvarhazard,int *t,int start,int stop);
void prodlim_comprisk(double* y,int* status,int* cause,int* NS,double* time,double* nrisk,int* event,int* loss,double* surv,double* cuminc,double* cause_hazard,double* varcuminc,double* cuminc_temp,double* cuminc_lag,double* v1,double* v2,int *t,int start,int stop);
float randF(float f,float l);
int randI(int number);
int neworder (int *a, int *b);
  
