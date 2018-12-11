#include "common.h"
#include "spherical.h"
#include <cmath>
#include <cstdlib>

using namespace std;

extern ParticleList particle;
typedef vector<vector<vector<int> > > LinkList;
typedef vector<vector<vector<pair<int,int> > > > NeighbourList;

int nx,ny;

LinkList lc;
NeighbourList neighbours;

void make_link_cell();
bool is_valid_neighbour(int ix, int iy, int iix, int iiy);
void init_neighbours();

void make_link_cell()
{
	for (unsigned int ix=0; ix<lc.size(); ix++){
		for (unsigned int iy=0; iy<lc[ix].size(); iy++){
			lc[ix][iy].clear();
		}
	}
	for (unsigned int i=0; i<no_of_particles; i++){
		int ix = int(nx*(particle[i].x()-x_0)/lx);
		int iy = int(ny*(particle[i].y()-y_0)/ly);
		if ((ix>=0) && (ix<nx) && (iy>=0) && (iy<ny)){
			lc[ix][iy].push_back(i);
		} else {
			cout << "Particle "<< i << " outside of the simulation area" << endl;
			exit(0);
		}
	}
}

void init_neighbours()
{
	int iix, iiy;

	neighbours.resize(nx);
	for (int ix=0; ix<nx; ix++){
		neighbours[ix].resize(ny);
	}
	for (int ix=0; ix<nx; ix++){
		for (int iy=0; iy<ny; iy++){
			for (int dx=-1; dx<=1; dx++){
				for (int dy=-1; dy<=1; dy++){
					iix = (ix+dx+nx)%nx;
					iiy = (iy+dy+ny)%ny;
					if (is_valid_neighbour(ix,iy,iix,iiy)){
						neighbours[ix][iy].push_back(pair<int,int>(iix,iiy));
					}
				}
			}
		}
	}
}

bool is_valid_neighbour(int ix, int iy, int iix, int iiy)
{
	if ((iix==(ix-1+nx)%nx) && (iiy==(iy+1+ny)%ny)) return true;
	if ((iix==(ix+nx)%nx) && (iiy==(iy+1+ny)%ny)) return true;
	if ((iix==(ix+1+nx)%nx) && (iiy==(iy+1+ny)%ny)) return true;
	if ((iix==(ix+1+nx)%nx) && (iiy==(iy+ny)%ny)) return true;
	return false;
}

void make_forces()
{
	for (unsigned int ix=0; ix<lc.size(); ix++){
		for (unsigned int iy=0; iy<lc[ix].size(); iy++){
			for (unsigned int j=0; j<lc[ix][iy].size(); j++){
				int pj = lc[ix][iy][j];
				for (unsigned int k=j+1; k<lc[ix][iy].size(); k++){
					int pk = lc[ix][iy][k];
					force(particle[pj],particle[pk],lx,ly,lz);
				}
				for (unsigned int n=0; n<neighbours[ix][iy].size(); n++){
					int iix = neighbours[ix][iy][n].first;
					int iiy = neighbours[ix][iy][n].second;
					for (unsigned int k=0; k<lc[iix][iiy].size(); k++){
						int pk = lc[iix][iiy][k];
						force(particle[pj],particle[pk],lx,ly,lz);
					}
				}
			}
		}
	}
}

void step()
{
	make_link_cell();
	integrate();
}

void init_algorithm()
{
	nx = int(lx/(4*radii_max));
	ny = int(ly/(4*radii_max));
	lc.resize(nx);
	for (unsigned int ix=0; ix<lc.size(); ix++) lc[ix].resize(ny);
	init_neighbours();
	make_link_cell();
}




















