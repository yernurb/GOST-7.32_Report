#include "common.h"
#include "spherical.h"
#include <cmath>
#include <cstdlib>

using namespace std;

extern ParticleList particle;

ParticleList safe;
double Timesafe;

double verlet_ratio = 0.6, verlet_distance = 0.0025, verlet_grid = 30.0;
double verlet_increase = 1.1;
int vnx,vny;
double dx,dy;

typedef vector<vector<vector<int> > > CellListType;
typedef vector<set<int> > VerletList;

CellListType celllist;
VerletList verlet;

bool verlet_needs_update();
bool make_verlet();
bool do_touch(int i, int k);

void make_forces()
{
        for(unsigned int i=0; i<particle.size()-1; i++){
		set<int>::iterator iter;
		for(iter=verlet[i].begin(); iter!=verlet[i].end(); iter++){
			force(particle[i], particle[*iter], lx, ly, lz);
		}
	}	
}

bool make_verlet()
{
	bool ok = true;

	verlet.resize(no_of_particles);
	for (int ix=0; ix<vnx; ix++){
		for (int iy=0; iy<vny; iy++){
			celllist[ix][iy].clear();
		}
	}
	for (unsigned int i=0; i<no_of_particles; i++){
		int ix = int((particle[i].x()-x_0)/dx);
		int iy = int((particle[i].y()-y_0)/dy);
		celllist[ix][iy].push_back(i);
	}
	for (unsigned int i=0; i<no_of_particles; i++){
		set<int> oldverlet = verlet[i];
		verlet[i].clear();
		int ix = int((particle[i].x()-x_0)/dx);
                int iy = int((particle[i].y()-y_0)/dy);
		for (int iix=ix-1; iix<=ix+1; iix++){
			for (int iiy=iy-1; iiy<=iy+1; iiy++){
				int wx = (iix+vnx)%vnx;
				int wy = (iiy+vny)%vny;
				for (unsigned int k=0; k<celllist[wx][wy].size(); k++){
					int pk = celllist[wx][wy][k];
					if (pk<(int)i){
						if (Distance(particle[i],particle[pk],lx,ly,lz)<
						    particle[i].r()+particle[pk].r()+verlet_distance){
							if (particle[i].ptype()==0||particle[pk].ptype()==0){
								verlet[i].insert(pk);
								if (oldverlet.find(pk)==oldverlet.end()){
									if (do_touch(i,pk)){
										ok = false;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return ok;
}

void init_algorithm()
{
	safe = particle;
	Timesafe = Time;
	verlet_grid = 5.0;
	vnx = int(lx/verlet_grid);
	vny = int(ly/verlet_grid);
	if (vnx==0) vnx = 1;
	if (vny==0) vny = 1;
	dx = lx/vnx;
	dy = ly/vny;
	celllist.resize(vnx);
	for(int i=0; i<vnx; i++) celllist[i].resize(vny);
	make_verlet();
}

bool do_touch(int i, int k)
{
	return (Distance(particle[i],particle[k],lx,ly,lz)<
		particle[i].r()+particle[k].r());
}

bool verlet_needs_update()
{
	for (unsigned int i=0; i<no_of_particles; i++){
		if (Distance(particle[i],safe[i],lx,ly,lz)>=verlet_ratio*verlet_distance){
			return true;
		}
	}
	return false;
}

void step()
{
	bool ok=true, newverlet=false;
	
	if (verlet_needs_update()){
		ok = make_verlet();
		newverlet = true;
	}
	if (!ok){
		cout << "fail: going back from " << Time << " to " << Timesafe << endl;
		particle = safe;
		Time = Timesafe;
		verlet_distance *= verlet_increase;
		make_verlet();
	}
	if (newverlet && ok){
		safe = particle;
		Timesafe = Time;
	}

        integrate();
}















