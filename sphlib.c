#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "utility.h"
#include "sphlib.h"

/* Computes the W function given h and the distance r */
double dW_r (double r, double h){

    if((0.<=r) && (r<=h))
        return r*(3.*r-4.*h)/(2.*pow4(h));
    else if((h<r) && (r<=2*h))
      return -pow2(r-2.*h)/(2.*pow4(h));
    else
        return 0;
}



/* Computes the coefficients of particle i into the particle structure to
 * perform later euler step evolution.
 * Variables:
 *  -p = array of particles.
 *  -i = index of particle we want to compute the coefficients of.
 *  -N = total number of particles.
 *  -NSHP = particles to consider to define the distance h.
 */ 
#define P(j)    ((GAMMA-1.)*p[j].rho*p[j].e)
#define Rho(j)  (p[j].rho)

void get_coeff(particle_t *p, long i, long N, long NSHP){
    long j,k;           //dummy variables
    long *dist_index;   //list to contain index of ordered distances
    double *distance;   //distance |x_i-x_j| from particle i to particle j
    double h,h2;        //h distance and two time h parameter
    double sign;        //sign variable for x_i-x_j
    double PiPkrho, vik; // minimize number of calculations

    /* MEMORY ALLOCATION */
    distance = (double*) malloc(N*sizeof(double)); assert(distance!=NULL);
    dist_index = (long*) malloc(N*sizeof(long)); assert(dist_index!=NULL);

    /* Compute the distance |x_i-x_j| */
    for(j=0;j<N;j++)
        distance[j] = fabs(p[i].x-p[j].x);

    /* Index of ordered distance computation */
    indexx(N,distance-1,dist_index-1);
    
    /* Initialization of loop for coefficients */
    j=1;                                //loop variable
    k=dist_index[j]-1;                  //index of variable j-closest to particle i
    p[i].v_coeff=0.;
    p[i].rho_coeff=0.;
    p[i].e_coeff=0.;
    h = distance[dist_index[NSHP]-1];     //h distance 
    h2 = 2.*h;
    
    /* While we have particles and we are inside h2, particles contribute to the coefficients */
    while(j<N && distance[k]<h2){
        sign = SIGN(p[i].x-p[k].x); //sign computation needed for performing derivative
      
        PiPkrho = P(k)/pow2(Rho(k))+P(i)/pow2(Rho(i));
        vik     = p[i].v-p[k].v;
      
        /* Coefficient modification for each loop*/
        p[i].v_coeff -= p[k].m*PiPkrho*dW_r(distance[k],h)*sign;
        p[i].rho_coeff += p[k].m*vik*dW_r(distance[k],h)*sign;
        p[i].e_coeff += p[k].m*PiPkrho*vik*dW_r(distance[k],h)*sign;
        /* Loop variables update */
        j++;
        k=dist_index[j]-1;
    }
    
    p[i].e_coeff /= 2.;
    
    /* Free memory */
    free(distance);
    free(dist_index);
}

#undef P
#undef Rho

/* Evolves the i-th particle using the coefficients provided
 * in the structure. 
 * Variables:
 *  -p = array of particle structures.
 *  -i = index of the particle to compute its evolution.
 *  -dt = time differential.
 */
void evolve(particle_t *p, long i, double dt){
    p[i].x += dt*p[i].v;
    p[i].v += dt*p[i].v_coeff;
    p[i].rho += dt*p[i].rho_coeff;
    p[i].e += dt*p[i].e_coeff;
}

/* Prints into file positon and density of all particles.
 * Variables:
 *  -P = array of particle structures.
 *  -N = size of the array.
 *  -t = current time of the system.
 */
void print_posdens(particle_t *P, long N, double t){
    char name[30];
    int j;
    FILE* fout;

    sprintf(name,"t_%6.4lf.dat",t);
    
    fout = fopen(name,"w");
    
    for(j=0;j<N;j++)
        fprintf(fout,"%lf %lf\n",P[j].x,P[j].rho);

    fclose(fout);
}
