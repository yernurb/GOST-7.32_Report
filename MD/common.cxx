#include "common.h"
#include "spherical.h"
#include <cmath>
#include <cstdlib>

using namespace std;

void init_system()
{
	x_0 = 0.0; y_0 = 0.0; z_0 = 0.0;
	lx = 2000; ly = 2000; lz = 0;
	timestep = 0.01;
	no_of_particles = N;
	radii_max = 2.5;
	particle.resize(no_of_particles);
	cout << "Filling factor: " << no_of_particles*pi*radii_max*radii_max/(lx*ly) << endl;
	cout << "Number of particles: " << no_of_particles << endl;

	for (unsigned int i=0; i<N; i++)
	{
		particle[i].set_x((lx-x_0)*drand48());
		particle[i].set_y((ly-y_0)*drand48());
		particle[i].set_z(0.0);
		particle[i].set_vx(50*(2*drand48()-1));
                particle[i].set_vy(50*(2*drand48()-1));
                particle[i].set_vz(double(0));
                particle[i].set_r(radii_max);
                particle[i].set_m(2.0);
                particle[i].set_Y(1000.0);
                particle[i].set_A(0.01);
		particle[i].set_nu(0.25);
                particle[i].set_ptype(0);
		particle[i].set_gamma(0.2);
	}

/*		
	for (unsigned int i=0; i<100; i++)
	{
		particle[i].set_x(100);
		particle[i].set_y(10+5*i);
		particle[i].set_z(0.0);
		particle[i].set_vx(double(0));
                particle[i].set_vy(double(0));
                particle[i].set_vz(double(0));
                particle[i].set_r(2.5);
                particle[i].set_m(2.0);
                particle[i].set_Y(500.0);
                particle[i].set_A(0.05);
                particle[i].set_ptype(1);
	}
	for (unsigned int i=100; i<200; i++)
        {
                particle[i].set_x(400);
                particle[i].set_y(10+5*(i-100));
                particle[i].set_z(0.0);
                particle[i].set_vx(double(0));
                particle[i].set_vy(double(0));
                particle[i].set_vz(double(0));
                particle[i].set_r(2.5);
                particle[i].set_m(2.0);
                particle[i].set_Y(500.0);
                particle[i].set_A(0.05);
                particle[i].set_ptype(1);
        }
	for (unsigned int i=200; i<N_inactive; i++)
        {
                particle[i].set_x(100+5*(i-200));
                particle[i].set_y(5);
                particle[i].set_z(0.0);
                particle[i].set_vx(double(0));
                particle[i].set_vy(double(0));
                particle[i].set_vz(double(0));
                particle[i].set_r(2.5);
                particle[i].set_m(2.0);
                particle[i].set_Y(500.0);
                particle[i].set_A(0.05);
                particle[i].set_ptype(1);
        }
	
	for (unsigned int i=N_inactive; i<N-1; i+=30)
		for (unsigned int k=0; k<30; k++)
		{
			particle[i+k].set_x(double(110+8*k));
			particle[i+k].set_y(double(150+8*(int((i-N_inactive)/30))));
			particle[i+k].set_z(double(0));
			particle[i+k].set_vx(double(0));
			particle[i+k].set_vy(double(0));
			particle[i+k].set_vz(double(0));
			particle[i+k].set_r(2.0);
			particle[i+k].set_m(2.0);
			particle[i+k].set_Y(500.0);
			particle[i+k].set_A(0.05);
			particle[i+k].set_ptype(0);		
		}
	
*/
/*
// Test of the restitution coeffitient (only two particles)
	particle[0].set_x(double(100));
	particle[0].set_y(double(250));
	particle[0].set_z(double(0));
	particle[0].set_vx(double(5));
	particle[0].set_vy(double(0));
	particle[0].set_vz(double(0));
	particle[0].set_r(15.0);
	particle[0].set_m(30.0);
	particle[0].set_Y(1000.0);
	particle[0].set_A(0.01);
	particle[0].set_nu(0.25);
	particle[0].set_ptype(0);
	particle[0].set_gamma(6.0);

	particle[1].set_x(double(135));
	particle[1].set_y(double(250));
	particle[1].set_z(double(0));
	particle[1].set_vx(double(-5));
	particle[1].set_vy(double(0));
	particle[1].set_vz(double(0));
	particle[1].set_r(15.0);
	particle[1].set_m(30.0);
	particle[1].set_Y(1000.0);
	particle[1].set_A(0.01);
	particle[1].set_nu(0.25);
	particle[1].set_ptype(0);
	particle[1].set_gamma(6.0);
*/
/*
	for (unsigned int i=0; i<particle.size(); i++)
        {
                particle[i].set_x(100+4.99*i);
                particle[i].set_y(100);
                particle[i].set_z(0.0);
                particle[i].set_vx(double(0));
                particle[i].set_vy(double(0));
                particle[i].set_vz(double(0));
                particle[i].set_r(2.5);
                particle[i].set_m(2.0);
                particle[i].set_Y(1000.0);
                particle[i].set_A(0.01);
                particle[i].set_ptype(0);
		particle[i].set_nu(0.25);
		particle[i].set_gamma(15.0);
        }

	particle[particle.size()-1].set_x(225);
	particle[particle.size()-1].set_y(400);
	particle[particle.size()-1].set_vy(-double(10));

*/

}

void integrate()
{
        for (unsigned int i=0; i<particle.size(); i++)
        {
                if (particle[i].ptype()==0){
                        particle[i].set_force_to_zero();
                        particle[i].predict(timestep);
                } else {
                        particle[i].boundary_conditions(i,timestep,Time);
                }
        }
        make_forces();
        for (unsigned int i=0; i<particle.size(); i++){
                if (particle[i].ptype()==0) particle[i].correct(timestep);
        }
        for (unsigned int i=0; i<particle.size(); i++){
                particle[i].periodic_bc(x_0,y_0,z_0,lx,ly,lz);
		//              particle[i].box_bc(x_0,y_0,z_0,lx,ly,lz);
        }
        Time += timestep;
}






