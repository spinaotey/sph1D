# SPH1D (Smoothed Particle Hydrodynamics in one spatial dimension)

## Introduction

This softwers computes the evolution of particles in a
smoothed 1-dimensional hydrodynamic system and prints
their position and energy and density for
certain time values into files. For the
dynamics of the system, one considers an
effective volume which is determined by
the closeness of the neighbour particles
given the parameter `NSHP`. 

OpenMP has been used to parallelize some of the steps performed.

## Compilation

Just use the command `make`.

## Usage

By executing the program `./sph` without any additional arguments, 
we are shown that it is used as `./sph  Npart Nsteps NSHP t_max`.


The arguments are the following:
* `Npart`: number of gas particles.
* `Nsteps`: number of integration steps.
* `NSHP`: number of SPH neighbour particles used with the solver.
* `t_max`: end-time of the integration, Tstart=0.
