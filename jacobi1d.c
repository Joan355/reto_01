#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/shm.h>
#include <unistd.h>
#include <sys/types.h> 
#include <math.h>
#include <sys/wait.h>
#include "timing.h"
#define CORES (int)(sysconf(_SC_NPROCESSORS_CONF) / 2)
/* --
 * Do nsweeps sweeps of Jacobi iteration on a 1D Poisson problem
 * 
 *    -u'' = f
 *
 * discretized by n+1 equally spaced mesh points on [0,1].
 * u is subject to Dirichlet boundary conditions specified in
 * the u[0] and u[n] entries of the initial vector.
 */
void jacobi(int nsweeps, int n, double* u, double* f, int start, int end)
{
    //printf("start -> %i end -> %i\n", start, end);
    //exit(0);
    int i, sweep;
    double h  = 1.0 / n;
    double h2 = h*h;
    double* utmp = (double*) malloc( (n+1) * sizeof(double) );

    /* Fill boundary conditions into utmp */
    utmp[0] = u[0];
    utmp[n] = u[n];

    for (sweep = 0; sweep < nsweeps; sweep += 2) {
        
        /* Old data in u; new data in utmp */
        for (i = start + 1; i < end + 1; ++i){
            utmp[i] = (u[i-1] + u[i+1] + h2*f[i])/2;
            
        }
        
        /* Old data in utmp; new data in u */
        for (i = start + 1; i < end + 1; ++i)
            u[i] = (utmp[i-1] + utmp[i+1] + h2*f[i])/2;
            //printf("%f\n",u[i]);
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
    int n, nsteps;
    double* u;
    double* f;
    double h;
    timing_t tstart, tend;
    char* fname;

    /* Process arguments */
    n      = (argc > 1) ? atoi(argv[1]) : 100;
    nsteps = (argc > 2) ? atoi(argv[2]) : 100;
    fname  = (argc > 3) ? argv[3] : NULL;
    h      = 1.0/n;


    //Crear memoria compartida
    int shm_id_u = shmget(IPC_PRIVATE, (n + 1) * sizeof(double), IPC_CREAT | 0666);
    if(shm_id_u < 0){
        perror("Error en shmget");
        return 1;
    }


    /* Allocate and initialize arrays */
    //u = (double*) malloc( (n+1) * sizeof(double) );
    u = (double *)shmat(shm_id_u, NULL, 0);
    f = (double*) malloc( (n+1) * sizeof(double) );
    memset(u, 0, (n+1) * sizeof(double));
    for (i = 0; i <= n; ++i)
        f[i] = i * h;


    //Calculando segementos de vector
    //Todo: Aqui se calcula los segementos que se repartiran entre los diferentes procesos
    int p_r = (int)floor(n / CORES);
    int step = p_r == 0 ? 1 : p_r;
    int r_step = n % CORES;
    int lim = (p_r == 0) ? r_step : CORES;


    /* Run the solver */
    //Obtener tiempo de inicio de proceso
    get_time(&tstart);

    //calculo de jacobiano

    pid_t pid = 1;

    for(int i = 0; i < CORES; i++){
        if(pid){
            pid = fork();
        }else if(pid < 0){
            return 1;
        }

        if(pid == 0){
            clock_t child_ts = clock();
            jacobi(nsteps, n, u, f,
            i*step,
            ((CORES - 1) == i) ? (i*step + r_step + step): (i*step + step));
            clock_t child_fs = clock();
            double tt = (double)(child_fs - child_ts) / CLOCKS_PER_SEC;
            printf("Tiempo -> %f Process -> %d\n",tt,getpid());
            exit(0);
        }

    }

    for(int i = 0; i < CORES; i++){
        wait(NULL);
    }
    //Obtener tiempo de finalizacion de proceso
    get_time(&tend);

    /* Run the solver */    
    printf("n: %d\n"
           "nsteps: %d\n"
           "Elapsed time: %g s\n", 
           n, nsteps, timespec_diff(tstart, tend));

    /* Write the results */
    if (fname)
        write_solution(n, u, fname);

    free(f);
    shmdt(u);
    shmctl(shm_id_u, IPC_RMID, NULL);
    return 0;
}
