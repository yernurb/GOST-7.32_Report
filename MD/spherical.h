#ifndef _SPHERICAL_H_
#define _SPHERICAL_H_

#include "Vector3d.h"

using namespace std;

inline double normalize(double dx, double L)
{
	if (fabs(dx)<0.00001) return dx;
	while (dx < -L/2) dx += L;
	while (dx >= L/2) dx -= L;
	return dx;
}

class spherical{
	friend istream & operator >> (istream & is, spherical & p){
		is >> p.rtd0 >> p.rtd1 >> p._r >> p._m >> p._ptype >> p.Y >> p.A >> p._force
                   >> p.rtd2 >> p.rtd3 >> p.rtd4;
		return is;
	}
	friend ostream & operator << (ostream & os, const spherical & p){
		os << p.rtd0 << " " << p.rtd1 << " " << p._r << " " << p._m << " " << p._ptype << " ";
	        os << p.Y << " " << p.A << " " << p._force << " ";
	        os << p.rtd2 << " " << p.rtd3 << " " << p.rtd4 << endl;
	        return os;
	}
	friend double Distance(const spherical & p1, const spherical & p2, double lx, double ly, double lz){
		double dx = normalize(p1.rtd0.x()-p2.rtd0.x(),lx);
		double dy = normalize(p1.rtd0.y()-p2.rtd0.y(),ly);
		double dz = normalize(p1.rtd0.z()-p2.rtd0.z(),lz);
		return sqrt(dx*dx+dy*dy+dz*dz);
	}
	friend Vector3d distvect(const spherical & p1, const spherical & p2, double lx, double ly, double lz){
		return Vector3d(p1.rtd0.x()-p2.rtd0.x(), p1.rtd0.y()-p2.rtd0.y(), p1.rtd0.z()-p2.rtd0.z());
	}
	friend Vector3d relvel(const spherical & p1, const spherical & p2){
		return Vector3d(p1.vx()-p2.vx(), p1.vy()-p2.vy(), p1.vz()-p2.vz());
	}
	friend void force(spherical & p1, spherical & p2, double lx, double ly, double lz);

public:
	spherical(): rtd0(null), rtd1(null), rtd2(null), rtd3(null), rtd4(null) {}
	Vector3d position() const { return rtd0; } // returns the position of the particle
	double x() const { return rtd0.x(); }
	double y() const { return rtd0.y(); }
	double z() const { return rtd0.z(); }
	Vector3d velocity() const { return rtd1; } // returns the velocity of the particle
	double vx() const { return rtd1.x(); }
	double vy() const { return rtd1.y(); }
	double vz() const { return rtd1.z(); }
	double r() const { return _r; }
	double m() const { return _m; }
	int ptype() const { return _ptype; }
	int species() const { return _species; }
	void add_force(const Vector3d & f) { _force += f; }
	void predict(double dt);
	void correct(double dt);
	void set_force_to_zero() { _force = null; }
	double kinetic_energy() const;
	void boundary_conditions(int n, double timestep, double Time);
	void periodic_bc(double x_0, double y_0, double z_0, double lx, double ly, double lz);
	void box_bc(double x_0, double y_0, double z_0, double lx, double ly, double lz);
	void set_x(double xval) { rtd0.set_x(xval); }
	void set_y(double yval) { rtd0.set_y(yval); }
	void set_z(double zval) { rtd0.set_z(zval); }
	void set_vx(double vxval) { rtd1.set_x(vxval); }
	void set_vy(double vyval) { rtd1.set_y(vyval); }
	void set_vz(double vzval) { rtd1.set_z(vzval); }
	void set_r(double rval) { _r = rval; }
	void set_m(double mval) { _m = mval; }
	void set_ptype(int type) { _ptype = type; }
	void set_Y(double Yval) { Y = Yval; }
	void set_A(double Aval) { A = Aval; }
	void set_nu(double nuval) { nu = nuval; }
	void set_gamma(double gammaval) { gamma = gammaval; }
	void set_species(int spec) { _species = spec; }	

private:
	double _r, _m;
	int _ptype, _species;
	double Y,A,nu,gamma;
	Vector3d rtd0, rtd1, rtd2, rtd3,rtd4;
	Vector3d _force;	
};

#endif
