#include "spherical.h"
#include <cstdlib>

using namespace std;

extern Vector3d G;

void spherical::predict(double dt)
{
	double a1 = dt;
	double a2 = a1*dt/2; 
	double a3 = a2*dt/3;
	double a4 = a3*dt/4;

	rtd0 += a1*rtd1 + a2*rtd2 + a3*rtd3 + a4*rtd4;
	rtd1 += a1*rtd2 + a2*rtd3 + a3*rtd4;
	rtd2 += a1*rtd3 + a2*rtd4;
	rtd3 += a1*rtd4;
}

void spherical::correct(double dt)
{
	static Vector3d accel, corr;
	Vector3d shear;
	double dtrez = 1/dt;
	const double coeff0 = double(19)/double(90)*(dt*dt/double(2));
	const double coeff1 = double(3)/double(4)*(dt/double(2));
	const double coeff3 = double(1)/double(2)*(double(3)*dtrez);
	const double coeff4 = double(1)/double(12)*(double(12)*(dtrez*dtrez));

	accel = (1/_m)*_force+G;
	
	corr = accel-rtd2;

	rtd0 += coeff0*corr;
	shear.set_x(-1.5*0.0001*(rtd0.y()-1000)); shear.set_y(0.0); shear.set_z(0.0);
	rtd1 += coeff1*corr;
	rtd2 = accel;
	rtd3 += coeff3*corr;
	rtd4 += coeff4*corr;
}

void force(spherical & p1, spherical & p2, double lx, double ly, double lz)
{
	double dx = normalize(p1.x()-p2.x(),lx);
	double dy = normalize(p1.y()-p2.y(),ly);
	double dz = normalize(p1.z()-p2.z(),lz);
	double rr = sqrt(dx*dx+dy*dy+dz*dz);
	double r1 = p1.r();
	double r2 = p2.r();
	double xi = r1+r2-rr;

	if (xi>0)
	{
		Vector3d dv = relvel(p1,p2);
	        Vector3d edd(dx/rr,dy/rr,dz/rr);
        	double xidot = -scalprod3d(edd,dv);
	        double D1 = (double(1)-p1.nu*p1.nu)/p1.Y;
        	double D2 = (double(1)-p2.nu*p2.nu)/p2.Y;
	        double D = 0.5*(D1+D2);
        	double A = 0.5*(p1.A+p2.A);
	        double gamma = 0.5*(p1.gamma+p2.gamma);
        	double reff = (r1*r2)/(r1+r2);
	        double alpha = double(2)*sqrt(reff)/(double(3)*D);
        	double beta = double(4)*pi*gamma*reff*sqrt(reff)/D;
//	        double asep = exp(log(double(3)*pi*D*gamma*reff*reff/double(2))/double(3));
//        	double xisep = asep*asep/reff-sqrt(double(8)*pi*gamma*D*asep/double(3));
		double Fh = alpha*xi*sqrt(xi);
		double Fb = sqrt(beta*xi*sqrt(xi));
		double Fd = double(3)*A*xidot*alpha*sqrt(xi)/double(2);
		double fn = Fh-Fb+Fd;
		Vector3d fnormplus = fn*edd;
		Vector3d fnormmin  = -fn*edd;
		
		if (p1.ptype() == 0)
			p1.add_force(fnormplus);
		if (p2.ptype() == 0)
			p2.add_force(fnormmin);
	}
//	else if ((xi>xisep)&&(xidot>0))
}

double spherical::kinetic_energy() const
{
	return (_m/double(2))*(rtd1.x()*rtd1.x()+rtd1.y()*rtd1.y()+rtd1.z()*rtd1.z());
}

void spherical::boundary_conditions(int n, double timestep, double Time)
{
	switch(ptype()){
	case(0): break;
	case(1): break;
//	case(2): {
		
//	} break;
	default: {
		cerr  << "ptype: " << ptype() << " not implemented" << endl;
		exit(1);
	}
	}
}

void spherical::periodic_bc(double x_0, double y_0, double z_0, double lx, double ly, double lz)
{
	double x,y,z;
	x = rtd0.x(); y = rtd0.y(); z = rtd0.z();
	while (x < x_0) x += lx;
	while (x > x_0+lx) x -= lx;
	while (y < y_0) y += ly;
	while (y > y_0+ly) y -= ly;
	while (z < z_0) z += lz;
	while (z > z_0+lz) z -= lz;
	rtd0.set_x(x);
	rtd0.set_y(y);
	rtd0.set_z(z);
}

void spherical::box_bc(double x_0, double y_0, double z_0, double lx, double ly, double lz)
{
	double x,y,z;
	x = rtd0.x(); y = rtd0.y(); z = rtd0.z();
//	vx = rtd1.x(); vy = rtd1.y(); vz = rtd1.vz();
//	while (x < x_0) vx = -vx;
	while (x > x_0+lx) x -= lx;
//	while (y < y_0) vy = -vy;
	while (y > y_0+ly) y -= ly;
//	while (z < z_0) vz = -vz;
	while (z > z_0+lz) z -= lz;
//	rtd1.set_x(vx);
//	rtd1.set_y(vy);
//	rtd1.set_z(vz);
}




















