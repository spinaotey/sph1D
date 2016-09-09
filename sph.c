#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "utility.h"
#include "sphlib.h"

/* 1D SMOOTHED PARTICLE HYDRODYNAMICS
 * Computes the evolution of particles in a
 * smoothed 1D hydrodynamic system and prints
 * their position and energy and density for
 * certain time values into files. For the
 * dynamics of the system, one considers an
 * effective volume which is determined by
 * the closeness of the neighbour particles
 * given the parameter NSHP. Npart determines
 * the number of particles involved, Nsteps,
 * the number of time steps; and t_max, the
 * the time to evolve the system.
 */

int main(int argc, char *argv[]){
    /* VARIABLES */
    long Npart, Nsteps, NSHP;   //Input variables: num.  particles, time steps, neighbours
    double t_max;               //Maximum time
    double dt, dx;              //time and space steps
    long i,j;                   //dummy variables
    particle_t *P;              //particle structur vector
    
    /* INPUT READ */
    if (argc < 5
       || sscanf(argv[1],"%ld",  &Npart)!=1
       || sscanf(argv[2],"%ld",  &Nsteps)!=1
       || sscanf(argv[3],"%ld",  &NSHP)!=1
       || sscanf(argv[4],"%lf",  &t_max)!=1
       ) {
           fprintf(stderr,"%s  Npart Nsteps NSHP t_max\n",argv[0]);
           return FALSE;
    }

    /* MEMORY ALLOCATION */

    P = (particle_t*) malloc(Npart*sizeof(particle_t)); assert(P != NULL);

    /* INITIAL CONDITION FOR THE SYSTEM */
    dx = 1./(Npart-1);          //equally distributed 
    for(i=0;i<Npart;i++){
        P[i].m = 1.;
        P[i].x = dx*i;
        P[i].v = 0.;
        P[i].rho = P[i].m/dx;
        P[i].e = 1e-5;
    }
    P[(Npart-1)/2].e = 1.;      //thermalisation of central particle

    /* EVOLUTON */
    dt = t_max/Nsteps;
    
    print_posdens(P,Npart,0.);    //print initial state
    
    /*Time evolution*/
    for(i=1;i<=Nsteps;i++){
        
#pragma omp parallel for private(j) shared(Npart,NSHP,P) default(none) schedule(dynamic,50)
        for(j=0;j<Npart;j++)
            get_coeff(P,j,Npart,NSHP);  //Get coefficients for euler step for each particle
#pragma omp parallel for private(j) shared(dt,P,Npart) default(none)
        for(j=0;j<Npart;j++)
            evolve(P,j,dt);             //Compute euler step for each particle

        if(i==Nsteps/2 || i==Nsteps/4 || i==Nsteps/4*3) //print some intermediate states
            print_posdens(P,Npart,dt*i);
    }
    
    print_posdens(P,Npart,t_max);   //print final state
    
    /* Free memory */
    free(P);
    
    return TRUE;
}

