/************************************************/
/* Provides class Vector3d			*/
/* This is an implementation of 3 dimentinal	*/
/* vector v=(x,y,z) 				*/
/* The class provides 13 methods:		*/
/* consider a,b,c to be vectors and alpha to be */
/* a double variable				*/
/* 1)  c = a + b; Summation			*/
/* 2)  c = a - b; Subtraction			*/
/* 3)  c += a; Increment			*/
/* 4)  c -= a; Decrement			*/
/* 5)  c = alpha*a; Multiplication		*/
/* 6)  c *= alpha; Renormalization		*/
/* 7)  alpha = norm3d(a); Calculate norm	*/
/* 8)  alpha = scalprod3d(a,b); Scalar product	*/
/* 9)  c = vecprod3d(a,b); Vector product	*/
/* 10) istream >> a; Writes from istream	*/
/* 11) ostream << a; Wrties into ostream	*/
/* 12) c = null; Zero vector			*/
/* 13) c = unitvec3d(a); Unit vector of a	*/
/************************************************/


#ifndef _VECTOR3D_H_
#define _VECTOR3D_h_

#include <cmath>
#include <iostream>

using namespace std;

#define pi M_PI

class Vector3d{
	friend istream & operator >> (istream & is, Vector3d & v){
		is >> v._x >> v._y >> v._z;
		return is;
	}
	friend ostream & operator << (ostream & os, const Vector3d & v){
		os << v.x() <<" "<< v.y() <<" "<< v.z();
		return os;
	}
	friend Vector3d operator + (const Vector3d & v1, const Vector3d & v2){
		return Vector3d(v1._x+v2._x, v1._y+v2._y, v1._z+v2._z);
	}
	friend Vector3d operator - (const Vector3d & v1, const Vector3d & v2){
		return Vector3d(v1._x-v2._x, v1._y-v2._y, v1._z-v2._z);
	}
	friend Vector3d operator * (const double val, const Vector3d & v){
		return Vector3d(val*v._x, val*v._y, val*v._z);
	}
	friend double norm3d(const Vector3d & v) { // returns the norm of the vector v
		return sqrt(v._x*v._x+v._y*v._y+v._z*v._z);
	}
	friend Vector3d unitvec3d(const Vector3d v){ // returns the unit vector of v
		return Vector3d(v._x/norm3d(v), v._y/norm3d(v), v._z/norm3d(v));
	}
	friend double scalprod3d(const Vector3d & v1, const Vector3d & v2) { // scalar product of v1 and v2
		return v1._x*v2._x+v1._y*v2._y+v1._z*v2._z;
	}
	friend Vector3d vecprod3d(const Vector3d & v1, const Vector3d & v2) { // vector product of v1 and v2
		return Vector3d(v1._y*v2._z-v1._z*v2._y, v1._z*v2._x-v1._x*v2._z, v1._x*v2._y-v1._y*v2._x);
	}

public:
	explicit Vector3d(double x=0, double y=0, double z=0): _x(x), _y(y), _z(z){};

	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
	void set_x(double xval) { _x=xval; }
	void set_y(double yval) { _y=yval; }
	void set_z(double zval) { _z=zval; }	

	const Vector3d & operator += (const Vector3d & v){
		_x += v._x; _y += v._y; _z += v._z;
		return *this;
	}
	const Vector3d & operator -= (const Vector3d & v){
		_x -= v._x; _y -= v._y; _z -= v._z;
		return *this;
	}
	const Vector3d & operator *= (const double val){
		_x *= val; _y *= val; _z *= val;
		return *this;
	}
	

private:
	double _x;
	double _y;
	double _z;
};


const Vector3d null(0,0,0);

#endif
