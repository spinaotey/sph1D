
/*STRUCTURE FOR 1D PARTICLE*/
typedef struct {
double m;           //mass
double x;           //position
double v;           //velocity
double rho;         //mass density
double e;           //thermal energy
/* Coefficients for time evolution */
double v_coeff;
double rho_coeff;
double e_coeff;
} particle_t;

#define GAMMA 5./3.


/* Computes the W function given h and the distance r */
double dW_r (double r, double h);


/* Computes the coefficients of particle i into the particle structure to
 * perform later euler step evolution.
 * Variables:
 *  -p = array of particles.
 *  -i = index of particle we want to compute the coefficients of.
 *  -N = total number of particles.
 *  -NSHP = particles to consider to define the distance h.
 */ 
void get_coeff(particle_t *p, long i, long N, long NSHP);


/* Evolves the i-th particle using the coefficients provided
 * in the structure. 
 * Variables:
 *  -p = array of particle structures.
 *  -i = index of the particle to compute its evolution.
 *  -dt = time differential.
 */
void evolve(particle_t *p, long i, double dt);


/* Prints into file positon and density of all particles.
 * Variables:
 *  -P = array of particle structures.
 *  -N = size of the array.
 *  -t = current time of the system.
 */
void print_posdens(particle_t *P, long N, double t);
