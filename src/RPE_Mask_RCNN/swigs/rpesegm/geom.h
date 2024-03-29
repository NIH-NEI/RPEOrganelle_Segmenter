#ifndef geom_h
#define geom_h

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#ifndef uint64
typedef unsigned long long uint64;
#endif

#define AUNITMAX 252
#define AUNITHALF 126		// AUNITMAX/2
#define AUNIT90 63			// AUNITMAX/4
#define AUNIT270 189		// AUNITMAX*3/4

const int HOOD_SIZE_NEUMANN = 5;
const int HOOD_SIZE_MOORE = 9;
const int HOOD_SIZE_FATCROSS = 21;
const int HOOD_SIZE_RAD3 = 37;

const int HOOD3D_6 = 7;
const int HOOD3D_18 = 19;
const int HOOD3D_26 = 27;

struct Point
{
	int x, y;
	Point() {}
	Point(int _x, int _y) : x(_x), y(_y) {}
	// Point(Point &p) : x(p.x), y(p.y) {}
	double dist(Point &p1) {
		double dx(x - p1.x);
		double dy(y - p1.y);
		return sqrt(dx*dx + dy*dy);
	}
	double dist(int _x, int _y) {
		double dx(x - _x);
		double dy(y - _y);
		return sqrt(dx*dx + dy*dy);
	}
	int idist(Point &p1) {
		int dx(x - p1.x); if (dx < 0) dx = -dx;
		int dy(y - p1.y); if (dy < 0) dy = -dy;
		return dx > dy ? dx : dy;
	}
	void swap(Point &p1) {
		x ^= p1.x; p1.x ^= x; x ^= p1.x;
		y ^= p1.y; p1.y ^= y; y ^= p1.y;
	}
	bool equals(Point &p1) {
		return x==p1.x && y==p1.y;
	}
};

struct Boundary
{
	int xmin, ymin, xmax, ymax;
	Boundary() {}
	Boundary(int _xmin, int _ymin, int _xmax, int _ymax) :
		xmin(_xmin), ymin(_ymin), xmax(_xmax), ymax(_ymax) {}
	//Boundary(const Boundary &bnd) :
	//	xmin(bnd.xmin), ymin(bnd.ymin), xmax(bnd.xmax), ymax(bnd.ymax) {}
	void expand(int bord=1) {
		xmin -= bord; ymin -= bord;
		xmax += bord; ymax += bord;
	}
	void combo(Boundary &bnd) {
		if (xmin > bnd.xmin) xmin = bnd.xmin;
		if (ymin > bnd.ymin) ymin = bnd.ymin;
		if (xmax < bnd.xmax) xmax = bnd.xmax;
		if (ymax < bnd.ymax) ymax = bnd.ymax;
	}
	bool intersects(Boundary &bnd) {
		return xmin < bnd.xmax && xmax > bnd.xmin &&
			ymin < bnd.ymax && ymax > bnd.ymin;
	}
	bool IsInside(int x, int y) {
		return x>=xmin && x<=xmax && y>=ymin && y<=ymax;
	}
	bool IsInside(Point p) {
		return p.x>=xmin && p.x<=xmax && p.y>=ymin && p.y<=ymax;
	}
	long long area() { return ((long long)(ymax-ymin+1)) * (xmax-xmin+1); }
};

struct WalkAngle
{
	int ssin;
	double asin;
	int scos;
	double acos;
};

struct NbrPoint
{
	int dx, dy;
};

struct HSeg
{
	int y, xl, xr;
	HSeg() {}
	HSeg(int _y, int _xl, int _xr) :
		y(_y), xl(_xl), xr(_xr) {}
	bool less_than(HSeg &other) {
		if (y != other.y) return y < other.y;
		return xl < other.xl;
	}
};

struct ColorCounter
{
	int c;
	int cnt;
	ColorCounter() {}
	ColorCounter(int _c, int _cnt) : c(_c), cnt(_cnt) {}
};

struct Contour
{
	std::vector<Point> path;
	bool is_closed;
	Boundary bnd;
	double plen;
	Contour() {}
	Contour(std::vector<Point> &_path) {
		set_path(_path);
	}
	Contour(const std::vector<int> & flat_cont, int iStart=-1, int len=-1);
	void set_path(std::vector<Point> &_path);
	void _update_from_path();
	void optimize(double mindist=3., double maxdist=5.);

	bool low_trend(Point *pp1, Point *pp2, double mindist=3.);
	bool hi_trend(Point *pp3, Point *pp4, double mindist=3.);
	bool can_close(double maxdist=15., double mindist=4.);
	
	int isLeft(Point P0, Point P1, Point P2);
	bool IsInside(Point P);
};

struct Particle
{
	int x0, y0;
	Boundary bnd;
	std::vector<HSeg> fill;
	//
	Particle() {}
	//Particle(Particle &ptc) {
	//	copyfrom(ptc);
	//}
	Particle(Contour &cont) { fromContour(cont); }
	void clear() {
		x0 = y0 = 0;
		bnd = Boundary(0,0,0,0);
		fill.clear();
	}
	void copyfrom(Particle& ptc) {
		fill.resize(ptc.fill.size());
		for (size_t i=0; i<ptc.fill.size(); i++)
			fill[i] = ptc.fill[i];
		x0 = ptc.x0;
		y0 = ptc.y0;
		bnd = ptc.bnd;
	}
	//
	void fromContour(Contour &cont);
	int update_from_fill();
	bool IsInside(int x, int y);
	Point center_mass();
	int overlay_area(std::vector<HSeg> &other_fill);
	int overlay_area(Particle &other);
};

struct Slice : public Particle
{
	int area;
	Slice() : Particle() {}
};

//--- 3D stuff

struct Point3D
{
	int x, y, z;
	Point3D() {}
	Point3D(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}
	double dist(int _x, int _y, int _z) {
		double dx(x - _x);
		double dy(y - _y);
		double dz(z - _z);
		return sqrt(dx*dx + dy*dy + dz*dz);
	}
	bool equals(Point3D &p1) {
		return x==p1.x && y==p1.y && z==p1.z;
	}
};

struct Boundary3D
{
	int xmin, ymin, zmin, xmax, ymax, zmax;
	Boundary3D() {}
	Boundary3D(int _xmin, int _ymin, int _zmin, int _xmax, int _ymax, int _zmax) :
		xmin(_xmin), ymin(_ymin), zmin(_zmin), xmax(_xmax), ymax(_ymax), zmax(_zmax) {}
	//Boundary3D(Boundary3D &bnd) :
	//	xmin(bnd.xmin), ymin(bnd.ymin), zmin(bnd.zmin), xmax(bnd.xmax), ymax(bnd.ymax), zmax(bnd.zmax) {}
	Boundary boundary2d() { return Boundary(xmin, ymin, xmax, ymax); }
	void set2d(Boundary &bnd) {
		xmin = bnd.xmin;
		ymin = bnd.ymin;
		xmax = bnd.xmax;
		ymax = bnd.ymax;
	}
	void expand(int bord=1) {
		xmin -= bord; ymin -= bord; zmin -= bord;
		xmax += bord; ymax += bord; zmax += bord;
	}
	void combo(Boundary3D &bnd) {
		if (xmin > bnd.xmin) xmin = bnd.xmin;
		if (ymin > bnd.ymin) ymin = bnd.ymin;
		if (zmin > bnd.zmin) zmin = bnd.zmin;
		if (xmax < bnd.xmax) xmax = bnd.xmax;
		if (ymax < bnd.ymax) ymax = bnd.ymax;
		if (zmax < bnd.zmax) zmax = bnd.zmax;
	}
	void combo2d(Boundary &bnd) {
		if (xmin > bnd.xmin) xmin = bnd.xmin;
		if (ymin > bnd.ymin) ymin = bnd.ymin;
		if (xmax < bnd.xmax) xmax = bnd.xmax;
		if (ymax < bnd.ymax) ymax = bnd.ymax;
	}
	bool intersects(Boundary3D &bnd) {
		return xmin < bnd.xmax && xmax > bnd.xmin &&
			ymin < bnd.ymax && ymax > bnd.ymin &&
			zmin < bnd.zmax && zmax > bnd.zmin;
	}
	bool intersects2d(Boundary &bnd) {
		return xmin < bnd.xmax && xmax > bnd.xmin &&
			ymin < bnd.ymax && ymax > bnd.ymin;
	}
	bool IsInside(int x, int y, int z) {
		return x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax;
	}
};

struct NbrPoint3D
{
	int dx, dy, dz;
};

struct Particle3D
{
	int x0, y0, z0;
	Boundary3D bnd;
	std::vector<std::vector<HSeg>> fills;
	//
	Particle3D() {}
	void clear() {
		x0 = y0 = z0 = 0;
		bnd = Boundary3D(0,0,0,0,0,0);
		for (std::vector<HSeg> &fill : fills)
			fill.clear();
	}
	bool empty() {
		for (std::vector<HSeg> &fill : fills) {
			if (!fill.empty()) return false;
		}
		return true;
	}
	long long update_from_fill();
	long long overlay_volume(Particle3D &other);
	int overlay_area2d(Particle &other, int z) {
		if (!bnd.intersects2d(other.bnd) || fills[z].empty()) return 0;
		return other.overlay_area(fills[z]);
	}
	void add_slice(Particle &ptc, int z) {
		std::vector<HSeg>& my_fill = fills[z];
		for (HSeg& hs : ptc.fill)
			my_fill.push_back(hs);
		if (z < bnd.zmin) bnd.zmin = z;
		if (z > bnd.zmax) bnd.zmax = z;
		bnd.combo2d(ptc.bnd);
	}
	double iou_score(int max_gap=3);
};

struct Nucleus : public Particle3D
{
	int idx;
	long long vol;
	int z_lo, z_hi, zdelta;
	Nucleus() : Particle3D() {}
	long long update_from_fill(double afac=0.4);
};

struct Cell : public Particle3D
{
	int idx;
	long long vol;
	int nnucs;
	int z_lo, z_hi;
	int zmin, zmax;
	Cell() : Particle3D() {}
	long long update_from_fill();
};

struct CsvCellDataReader
{
	std::ifstream fin;
	bool eof = false;
	bool is_2d, is_3d;
	char *sbuf;
	int line = 0;
	
	CsvCellDataReader(const char *filename) : fin(filename, std::ios::in)
	{
		sbuf = new char[256];
		is_2d = is_3d = false;
		if (fin.is_open()) {
			fin.getline(sbuf, 256);
			if (fin.fail()) {
				eof = true;
			} else {
				if (!strncmp(sbuf, "ID,Frame,y,xL,xR", 16)) {
					is_3d = true;
				} else if (!strncmp(sbuf, "ID,y,xL,xR", 10)) {
					is_2d = true;
				} else eof = true;
				if (!eof) ++line;
			}
		} else {
			eof = true;
		}
	}
	virtual ~CsvCellDataReader() { delete [] sbuf; }
	int read_hs(int *id, int *fr, int *y, int *xl, int *xr) {
		if (eof) return 0;
		if (!fin.eof()) {
			fin.getline(sbuf, 256);
			if (fin.fail()) eof = true;
		} else eof=true;
		if (eof) return 0;
		if (is_3d) {
			if (sscanf(sbuf, "%d,%d,%d,%d,%d", id, fr, y, xl, xr) == 5) {
				return 5;
			}
		} else if (is_2d) {
			if (sscanf(sbuf, "%d,%d,%d,%d", id, y, xl, xr) == 4) {
				*fr = 0;
				return 4;
			}
		}
		eof = true;
		return 0;
	}
	void close() { fin.close(); }
};

inline int au_diff(int a1, int a2) {
	int d = a1 < a2 ? a2 - a1 : a1 - a2;
	return d <= AUNITHALF ? d : AUNITMAX - d;
}

inline void sort_fill(std::vector<HSeg> & fill) {
	std::sort(fill.begin(), fill.end(), [](HSeg &a, HSeg &b) {
        return a.less_than(b);   
    });
}
inline int fill_area(std::vector<HSeg> & fill) {
	int a = 0;
	for (HSeg &s : fill) a += (s.xr-s.xl+1);
	return a;
}

void reverse_path(std::vector<Point> & path);
void remove_lead_from_path(std::vector<Point> & path, size_t k);
double path_length(std::vector<Point> & path);
std::vector<Point> path_from_flat(const std::vector<int>& flat_cont, int iStart, int len);
Boundary boundary_around(int x, int y, int dist, int w, int h);

Boundary fill_boundary(std::vector<HSeg> & fill);
double fill_centroid(std::vector<HSeg> & fill, double *px, double *py);
double fill_circularity(std::vector<HSeg> & fill);
int fill_overlay_area(std::vector<HSeg> &fill, std::vector<HSeg> &other_fill);
void fill_from_contour(std::vector<HSeg> &fill, Contour &cont);

void write_particle_data(std::vector<Slice> particles, const char *outfn);
void write_cell_data(std::vector<Particle3D> &cells, const char *outfn);
void write_cell_data(std::vector<Cell> &cells, const char *outfn);
int read_cell_data(const char *infn, std::vector<Particle3D> &cells, int w, int h, int d);
int read_particle_data(const char *infn, std::vector<Particle> &particles, int w, int h);

#ifndef __geom_main__
extern NbrPoint hood_pts[37];
extern WalkAngle wangs[AUNITMAX];
extern NbrPoint3D hood3d_pts[27];
#endif

#endif
