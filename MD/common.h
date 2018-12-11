#ifndef _COMMON_H_
#define _COMMON_H_

#include "spherical.h"
#include <fstream>
#include <vector>
#include <set>
#include <utility>

using namespace std;

#define N 30000 // total number of particles: active + passive
#define N_inactive 0 // number of inactive particles

extern double lx, ly, lz;
extern double x_0, y_0, z_0;
extern double radii_max;

extern unsigned int no_of_particles;
extern double Time,timestep;
extern ofstream fphase, fenergy;

extern int nstep, nprint, nenergy;
typedef vector<spherical> ParticleList;
extern ParticleList particle;

void step();
void make_forces();
void integrate();
void init_algorithm();
void phase_plot(ostream & os);
void init_system();
void make_forces();


#endif
