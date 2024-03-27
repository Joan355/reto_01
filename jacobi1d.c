#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "timing.h"
#include <math.h>
#include <sys/time.h>
#include <time.h>
#define THREADS 6

typedef struct {
    int nsweeps;
    int n;
    double *u;
    double *f; 
    int start; 
    int end;
} ThreadData;
/* --
 * Do nsweeps sweeps of Jacobi iteration on a 1D Poisson problem
 * 
 *    -u'' = f
 *
 * discretized by n+1 equally spaced mesh points on [0,1].
 * u is subject to Dirichlet boundary conditions specified in
 * the u[0] and u[n] entries of the initial vector.
 */
void *jacobi(void *d)
{
    ThreadData *data = (ThreadData *)d;
    int i, sweep;
    double h  = 1.0 / data->n;
    double h2 = h*h;
    double* utmp = (double*) malloc( (data->n+1) * sizeof(double) );

    /* Fill boundary conditions into utmp */
    utmp[0] = data->u[0];
    utmp[data->n] = data->u[data->n];

    for (sweep = 0; sweep < data->nsweeps; sweep += 2) {
        
        /* Old data in u; new data in utmp */
        for (i = 1 + data->start; i < data->end + 1 ; ++i)
            utmp[i] = (data->u[i-1] + data->u[i+1] + h2*data->f[i])/2;
        
        /* Old data in utmp; new data in u */
        for (i = 1 + data->start; i < data->end + 1; ++i)
            data->u[i] = (utmp[i-1] + utmp[i+1] + h2*data->f[i])/2;
    }

    free(utmp);
}


void write_solution(int n, double* u, const char* fname)
{
    int i;
    double h = 1.0 / n;
    FILE* fp = fopen(fname, "w+");
    for (i = 0; i <= n; ++i)
        fprintf(fp, "%g %g\n", i*h, u[i]);
    fclose(fp);
}


int main(int argc, char** argv)
{
    int i;
    int n, nsteps, hilos;
    double* u;
    double* f;
    double h;
    timing_t tstart, tend;
    char* fname;

    /* Process arguments */
    n      = (argc > 1) ? atoi(argv[1]) : 100;
    nsteps = (argc > 2) ? atoi(argv[2]) : 100;
    hilos  = (argc > 3) ? atoi(argv[3]) : 1;
    fname  = (argc > 4) ? argv[3] : NULL;
    h      = 1.0/n;

    /* Allocate and initialize arrays */
    u = (double*) malloc( (n+1) * sizeof(double) );
    f = (double*) malloc( (n+1) * sizeof(double) );
    memset(u, 0, (n+1) * sizeof(double));
    for (i = 0; i <= n; ++i)
        f[i] = i * h;

    pthread_t threads[hilos];
    ThreadData thread_data[hilos];

    int t_r = (int)floor(n / hilos);
    int step = t_r == 0 ? 1 : t_r;
    int r_step = n % hilos;
    int lim = (t_r == 0) ? r_step : hilos;

    /* Run the solver */
    get_time(&tstart);
    for(int i = 0; i < hilos; i++){

        thread_data[i].nsweeps = nsteps;
        thread_data[i].n = n;
        thread_data[i].u = u;
        thread_data[i].f = f;
        thread_data[i].start = i * step;
        thread_data[i].end = ((hilos - 1) == i) ? (i*step + r_step + step): (i*step + step);
        pthread_create(&threads[i], NULL, jacobi, (void *)&thread_data[i]);
        
    }
    for (int i = 0; i < hilos; i++) {
        pthread_join(threads[i], NULL);
    }
    get_time(&tend);

    /* Run the solver */    
    printf(
           "%g\n", 
           timespec_diff(tstart, tend));

    /* Write the results */
    if (fname)
        write_solution(n, u, fname);

    free(f);
    free(u);
    return 0;
}
