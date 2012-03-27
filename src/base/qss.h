#ifndef __QSS_H__
#define __QSS_H__

typedef struct {
	int num_eqn;
	double *y0;
	double *y1;
	double *rs; 
	double *a0; 
	double *a1;
	double *w0;
	double *w1;
	double *buf;
	void (* rates)(double t, double *y, void *params, double *w, double *a);
	void (* adjust)(double t, double *y, void *params);
} qss_system;

qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ); 
void qss_free( qss_system *sys );
void qs1_step( qss_system *sys, double t, double dt, double yf[], void *params );
void qsn_step( qss_system *sys, double t, double dt, double yf[], void *params );
void qss_solve( qss_system *sys, double t_begin, double delta_t, double y[], const double err[], void *params );

#endif
