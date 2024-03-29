#define __geom_main__
#include "geom.h"

//----------------------------- Contour --------------------------------

Contour::Contour(const std::vector<int> & flat_cont, int iStart, int len)
{
	if (iStart < 0) {
		len = flat_cont[0];
		iStart = 1;
	}
	is_closed = false;
	path = path_from_flat(flat_cont, iStart, len);
	if (path.size() > 3) {
		path.push_back(path[0]);
		is_closed = true;
	}
	_update_from_path();
}

void Contour::set_path(std::vector<Point> &_path)
{
	is_closed = false;
	path.resize(_path.size());
	if (path.size() > 0) {
		for (size_t j=0; j<path.size(); j++)
			path[j] = _path[j];
		if (path.size() > 3 && path[0].equals(path[path.size()-1]))
			is_closed = true;
	}
	_update_from_path();
}

void Contour::_update_from_path()
{
	plen = 0.;
	if (path.size() > 0) {
		for(size_t j=1; j<path.size(); j++) {
			plen += path[j-1].dist(path[j]);
		}
		bnd.xmin = bnd.xmax = path[0].x;
		bnd.ymin = bnd.ymax = path[0].y;
		for (Point &p : path) {
			if (bnd.xmin > p.x) bnd.xmin = p.x;
			if (bnd.xmax < p.x) bnd.xmax = p.x;
			if (bnd.ymin > p.y) bnd.ymin = p.y;
			if (bnd.ymax < p.y) bnd.ymax = p.y;
		}
	} else {
		bnd.xmin = bnd.xmax = bnd.ymin = bnd.ymax = 0;
	}
}

void Contour::optimize(double mindist, double maxdist)
{
	if (path.size() < 4) return;
	std::vector<Point> npath;
	npath.reserve(path.size());
	int iLow=0, iHigh=1;
	while (size_t(iHigh) < path.size()) {
		Point &p0 = path[iLow];
		npath.push_back(p0);
		for (; size_t(iHigh) < path.size(); iHigh++) {
			if (p0.dist(path[iHigh]) >= mindist) break;
		}
		while (p0.dist(path[iHigh]) > maxdist) {
			if (iHigh-1 <= iLow) break;
			--iHigh;
		}
		if (iLow == iHigh) break;
		iLow = iHigh;
		++iHigh;
	}
	if (iLow < iHigh && size_t(iHigh) < path.size()) {
		npath.push_back(path[iHigh]);
	}
	if (is_closed && !npath[0].equals(npath[npath.size()-1])) {
		npath.push_back(npath[0]);
	}
	set_path(npath);
}

bool Contour::low_trend(Point *pp1, Point *pp2, double mindist)
{
	if (path.size() == 0) return false;
	Point &p0 = path[0];
	Point p1;
	for (size_t j=1; j<path.size(); j++) {
		p1 = path[j];
		if (p0.dist(p1) >= mindist) break;
	}
	*pp1 = p0;
	*pp2 = p1;
	return true;
}

bool Contour::hi_trend(Point *pp3, Point *pp4, double mindist)
{
	if (path.size() == 0) return false;
	Point &p0 = path[path.size()-1];
	Point p1;
	for (size_t j=path.size()-1; j>0; j--) {
		p1 = path[j-1];
		if (p0.dist(p1) >= mindist) break;
	}
	*pp3 = p1;
	*pp4 = p0;
	return true;
}

/*
         P5  P6
        :      :
      P2         P3
     /             \
   P1               P4
    +---------------+
*/

bool Contour::can_close(double maxdist, double mindist)
{
	if (path.size() < 4) return false;
	double dist = path[0].dist(path[path.size()-1]);
	if (dist <= mindist) return true;
	if (dist > maxdist) return false;
	Point p1, p2, p3, p4;
	if (!low_trend(&p1, &p2)) return false;
	if (!hi_trend(&p3, &p4)) return false;
	double gap = p2.dist(p3);
	if (p1.dist(p4) <= gap) return false;

	Point p5(2*p2.x-p1.x, 2*p2.y-p1.y);
	Point p6(2*p3.x-p4.x, 2*p3.y-p4.y);
	if (p5.dist(p6) >= gap) return false;
	
	p3.x = p2.x + p4.x - p3.x;
	p3.y = p2.y + p4.y - p3.y;
	double a = p1.dist(p2);
	double b = p2.dist(p3);
	double c = p1.dist(p3);
	return c*c*1.2 > a*a + b*b;
}

//-- "Winding Number" algorithm for the inclusion of a point in polygon, adapted from:
// http://geomalgorithms.com/a03-inclusion.html

// tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
int Contour::isLeft(Point P0, Point P1, Point P2)
{
	return ((P1.x - P0.x) * (P2.y - P0.y)
		- (P2.x - P0.x) * (P1.y - P0.y));
}
bool Contour::IsInside(Point P)
{
	if (!is_closed || path.size() < 4) return false;
	int n = int(path.size() - 1);
	std::vector<Point> &V = path;

	int wn = 0;    // the  winding number counter (=0 only when P is outside)

	// loop through all edges of the polygon
	for (int i = 0; i < n; i++) {		// edge from V[i] to V[i+1]
		if (V[i].y <= P.y) {			// start y <= P.y
			if (V[i + 1].y > P.y)		// an upward crossing
				if (isLeft(V[i], V[i + 1], P) > 0)	// P left of  edge
					++wn;				// have  a valid up intersect
		}
		else {							// start y > P.y (no test needed)
			if (V[i + 1].y <= P.y)		// a downward crossing
				if (isLeft(V[i], V[i + 1], P) < 0)	// P right of  edge
					--wn;				// have  a valid down intersect
		}
	}
	return wn != 0;
}
//-- End of "Winding Number" algorithm


//----------------------------- Particle --------------------------------

void Particle::fromContour(Contour &cont)
{
	fill.resize(0);
	bnd = cont.bnd;
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (cont.IsInside(Point(x, y))) {
				if (!is_in) {
					hs.xl = x;
					is_in = true;
				}
				hs.xr = x;
			} else {
				if (is_in) {
					fill.push_back(hs);
					is_in = false;
				}
			}
		}
		if (is_in) fill.push_back(hs);
	}
	if (fill.size() > 0) {
		x0 = fill[0].xl;
		y0 = fill[0].y;
	} else {
		x0 = y0 = 0;
	}
}

int Particle::update_from_fill()
{
	if (fill.size() > 0) {
		x0 = fill[0].xl;
		y0 = fill[0].y;
	} else {
		x0 = y0 = 0;
	}
	bnd = fill_boundary(fill);
	return fill_area(fill);
}

bool Particle::IsInside(int x, int y)
{
	if (!bnd.IsInside(x, y)) return false;
	for (HSeg &hs : fill) {
		if (hs.y != y) continue;
		if (x>=hs.xl && x<=hs.xr) return true;
	}
	return false;
}

Point Particle::center_mass()
{
	double xc = 0., yc = 0.;
	int npt = 0;
	for (HSeg &hs : fill) {
		int y = hs.y;
		for (int x=hs.xl; x<=hs.xr; x++) {
			xc += x;
			yc += y;
			++npt;
		}
	}
	if (npt == 0) return Point(0, 0);
	return Point(int(xc/npt), int(yc/npt));
}

int Particle::overlay_area(std::vector<HSeg> &other_fill)
{
	return fill_overlay_area(fill, other_fill);
}
int Particle::overlay_area(Particle &other)
{
	if (!bnd.intersects(other.bnd)) return 0;
	return fill_overlay_area(fill, other.fill);
}

//----------------------------- Particle3D --------------------------------

long long Particle3D::update_from_fill()
{
	long long vol = 0;
	x0 = y0 = z0 = 0;
	bnd = Boundary3D(0,0,0,0,0,0);
	bool is_first = true;
	for (int z=0; size_t(z)<fills.size(); z++) {
		std::vector<HSeg> &fill = fills[z];
		if (fill.empty()) continue;
		bnd.zmax = z;
		for (HSeg &hs : fill) {
			if (is_first) {
				is_first = false;
				x0 = hs.xl;
				y0 = hs.y;
				z0 = z;
				bnd.xmin = hs.xl;
				bnd.xmax = hs.xr;
				bnd.ymin = bnd.ymax = hs.y;
				bnd.zmin = z;
			} else {
				if (bnd.xmin > hs.xl) bnd.xmin = hs.xl;
				if (bnd.xmax < hs.xr) bnd.xmax = hs.xr;
				if (bnd.ymin > hs.y) bnd.ymin = hs.y;
				if (bnd.ymax < hs.y) bnd.ymax = hs.y;
			}
			vol += (hs.xr - hs.xl + 1);
		}
	}
	return vol;
}

long long Particle3D::overlay_volume(Particle3D &other)
{
	long long vol = 0;
	
	for (int z=0; size_t(z)<fills.size(); z++) {
		if (size_t(z) >= other.fills.size()) break;
		vol += fill_overlay_area(fills[z], other.fills[z]);
	}
	
	return vol;
}

double Particle3D::iou_score(int max_gap)
{
	double sc = 0.;
	int d = int(fills.size());
	for (int z1=0; z1<d-1; z1++) {
		std::vector<HSeg> &fill1 = fills[z1];
		if (fill1.empty()) continue;
		int a1 = fill_area(fill1);
		for (int z2=z1+1; z2<d; z2++) {
			if (z2 - z1 > max_gap) break;
			std::vector<HSeg> &fill2 = fills[z2];
			if (fill2.empty()) continue;
			int a2 = fill_area(fill2);
			int ovl = fill_overlay_area(fill1, fill2);
			if (ovl > 0)
				sc += double(ovl) / (a1 + a2 - ovl);
		}
	}
	return sc;
}

//----------------------------- Nucleus --------------------------------

long long Nucleus::update_from_fill(double afac)
{
	vol = 0;
	x0 = y0 = z0 = 0;
	bnd = Boundary3D(0,0,0,0,0,0);
	bool is_first = true;
	
	std::vector<int> amap(fills.size(), 0);
	int zmax=0, amax=0;
	
	for (int z=0; size_t(z)<fills.size(); z++) {
		std::vector<HSeg> &fill = fills[z];
		if (fill.empty()) continue;
		bnd.zmax = z;
		int a = amap[z] = fill_area(fill);
		vol += a;
		if (a > amax) {
			amax = a;
			zmax = z;
		}
		Boundary b = fill_boundary(fill);
		if (is_first) {
			is_first = false;
			bnd.zmin = z;
			bnd.set2d(b);
			HSeg &hs = fill[0];
			x0 = hs.xl;
			y0 = hs.y;
			z0 = z;
		} else {
			bnd.combo2d(b);
		}
	}
	
	int acut = int(amax * afac);
	if (acut < 1) acut = 1;
	for (z_lo=0; z_lo<zmax; z_lo++) {
		if (amap[z_lo] >= acut) break;
	}
	for (z_hi=int(amap.size()-1); z_hi>zmax; z_hi--) {
		if (amap[z_hi] >= acut) break;
	}
	zdelta = z_hi - z_lo + 1;
	return vol;
}

//----------------------------- Cell --------------------------------

long long Cell::update_from_fill()
{
	vol = 0;
	x0 = y0 = z0 = 0;
	bnd = Boundary3D(0,0,0,0,0,0);
	bool is_first = true;
	
	for (int z=0; size_t(z)<fills.size(); z++) {
		std::vector<HSeg> &fill = fills[z];
		if (fill.empty()) continue;
		bnd.zmax = z;
		vol += fill_area(fill);
		Boundary b = fill_boundary(fill);
		if (is_first) {
			is_first = false;
			bnd.zmin = z;
			bnd.set2d(b);
			HSeg &hs = fill[0];
			x0 = hs.xl;
			y0 = hs.y;
			z0 = z;
		} else {
			bnd.combo2d(b);
		}
	}
	return vol;
}

//----------------------------- Global utility functions --------------------------------

void reverse_path(std::vector<Point> & path)
{
	int iLow=0, iHigh=int(path.size()-1);
	while (iLow < iHigh) {
		Point tmp = path[iLow];
		path[iLow] = path[iHigh];
		path[iHigh] = tmp;
		++iLow;
		--iHigh;
	} 
}

void remove_lead_from_path(std::vector<Point> & path, size_t k)
{
	if (k == 0) return;
	size_t i = 0;
	for (; k<path.size(); k++) {
		path[i] = path[k];
		++i;
	}
	path.resize(i);
}

double path_length(std::vector<Point> & path)
{
	double plen = 0.;
	for(size_t j=1; j<path.size(); j++) {
		plen += path[j-1].dist(path[j]);
	}
	return plen;
}

std::vector<Point> path_from_flat(const std::vector<int>& flat_cont, int iStart, int len)
{
	std::vector<Point> path;
	for (int j=0; j<len; j++) {
		path.push_back(Point(flat_cont[iStart], flat_cont[iStart+1]));
		iStart += 2;
	}
	return path;
}

Boundary boundary_around(int x, int y, int dist, int w, int h)
{
	int sz = dist + dist;
	int xtop = w - 1;
	int xmin = x - dist;
	int xmax = x + dist;
	if (xmin < 0) {
		xmin = 0;
		xmax = sz;
		if (xmax > xtop) xmax = xtop;
	}
	if (xmax > xtop) {
		xmax = xtop;
		xmin = xmax - sz;
		if (xmin < 0) xmin = 0;
	}
	int ytop = h - 1;
	int ymin = y - dist;
	int ymax = y + dist;
	if (ymin < 0) {
		ymin = 0;
		ymax = sz;
		if (ymax > ytop) ymax = ytop;
	}
	if (ymax > ytop) {
		ymax = ytop;
		ymin = ymax - sz;
		if (ymin < 0) ymin = 0;
	}
	return Boundary(xmin, ymin, xmax, ymax);
}

Boundary fill_boundary(std::vector<HSeg> & fill)
{
	Boundary b(0,0,0,0);
	if (fill.size() == 0) return b;
	HSeg hs0 = fill[0];
	b.ymin = b.ymax = hs0.y;
	b.xmin = hs0.xl;
	b.xmax = hs0.xr;
	for (HSeg &s : fill) {
		if (b.ymin > s.y) b.ymin = s.y;
		if (b.ymax < s.y) b.ymax = s.y;
		if (b.xmin > s.xl) b.xmin = s.xl;
		if (b.xmax < s.xr) b.xmax = s.xr;
	}
	return b;
}

double fill_centroid(std::vector<HSeg> & fill, double *px, double *py)
{
	double xc = 0., yc = 0.;
	int npt = 0;
	for (HSeg &hs : fill) {
		int y = hs.y;
		for (int x=hs.xl; x<=hs.xr; x++) {
			xc += x;
			yc += y;
			++npt;
		}
	}
	if (npt == 0) {
		*px = *py = 0.;
		return 0.;
	}
	xc /= npt;
	yc /= npt;
	*px = xc;
	*py = yc;

	double r = 0.;
	for (HSeg &hs : fill) {
		double dy2 = yc - hs.y;
		dy2 *= dy2;
		for (int x=hs.xl; x<=hs.xr; x++) {
			double dx = xc - x;
			r += (dy2 + dx*dx);
		}
	}
	return sqrt(2.*r/npt);
}

double fill_circularity(std::vector<HSeg> & fill)
{
	double xc, yc;
	double rad = fill_centroid(fill, &xc, &yc);
	if (rad < 0.0000001) return 0.;
	double radsq = rad * rad;
	int inside = 0;
	for (HSeg &hs : fill) {
		double dy2 = yc - hs.y;
		dy2 *= dy2;
		for (int x=hs.xl; x<=hs.xr; x++) {
			double dx = xc - x;
			if (dx*dx + dy2 <= radsq)
				++inside;
		}
	}

	return double(inside) / (M_PI * radsq);
}

int fill_overlay_area(std::vector<HSeg> &fill, std::vector<HSeg> &other_fill)
{
	int a = 0;
	for (HSeg &hs0 : fill) {
		for (HSeg &hs1 : other_fill) {
			if (hs0.y != hs1.y) continue;
			int xmin = hs0.xl < hs1.xl ? hs1.xl : hs0.xl;
			int xmax = hs0.xr > hs1.xr ? hs1.xr : hs0.xr;
			if (xmax < xmin) continue;
			a += (xmax-xmin+1);
		}
	}
	return a;
}

void fill_from_contour(std::vector<HSeg> &fill, Contour &cont)
{
	Boundary bnd = cont.bnd;
	HSeg hs;
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		bool is_in = false;
		hs.y = y;
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			if (cont.IsInside(Point(x, y))) {
				if (!is_in) {
					hs.xl = x;
					is_in = true;
				}
				hs.xr = x;
			} else {
				if (is_in) {
					fill.push_back(hs);
					is_in = false;
				}
			}
		}
		if (is_in) fill.push_back(hs);
	}
}

void write_particle_data(std::vector<Slice> particles, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,y,xL,xR\n";
		int id = 1;
		
		for (Slice & ptc : particles) {
			for (HSeg &hs : ptc.fill) {
				fos << id << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
			}
			++id;
		}

		fb.close();
	}
}

void write_cell_data(std::vector<Particle3D> &cells, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,Frame,y,xL,xR\n";
		int id = 1;
		
		for (Particle3D &cell : cells) {
			if (cell.empty()) continue;
			for (int z=0; size_t(z)<cell.fills.size(); z++) {
				for (HSeg &hs : cell.fills[z]) {
					fos << id << "," << z << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
				}
			}
			++id;
		}

		fb.close();
	}
}

void write_cell_data(std::vector<Cell> &cells, const char *outfn)
{
	std::filebuf fb;
	fb.open(outfn, std::ios::out);
	if (fb.is_open())
	{
		std::ostream fos(&fb);
		fos << "ID,Frame,y,xL,xR\n";
		int id = 1;
		
		for (Cell &cell : cells) {
			if (cell.empty()) continue;
			for (int z=0; size_t(z)<cell.fills.size(); z++) {
				// if (cell.fills[z].empty()) continue;
				for (HSeg &hs : cell.fills[z]) {
					fos << id << "," << z << "," << hs.y << "," << hs.xl << "," << hs.xr << "\n";
				}
			}
			++id;
		}

		fb.close();
	}
}

int read_cell_data(const char *infn, std::vector<Particle3D> &cells, int w, int h, int d)
{
	CsvCellDataReader rdr(infn);
	if (rdr.eof || !rdr.is_3d) return -1;

// std::cout << "Reading: " << infn << std::endl;

	size_t j0 = cells.size();
	
	int cur_id = -1;
	Particle3D *pcell = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (z<0 || z>=d || y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			cells.resize(cells.size()+1);
			pcell = & cells[cells.size()-1];
			pcell->fills.resize(d);
			cur_id = id;
		}
		if (!pcell) continue;
		pcell->fills[z].push_back(HSeg(y, xl, xr));
	}
	
	for (size_t j=j0; j<cells.size(); j++) {
		Particle3D& cell = cells[j];
		for (int z=0; z<d; z++) cell.fills[z].shrink_to_fit();
		cell.fills.shrink_to_fit();
		cell.update_from_fill();
	}
	return cur_id;
}

int read_particle_data(const char *infn, std::vector<Particle> &particles, int w, int h)
{
	CsvCellDataReader rdr(infn);
	if (rdr.eof || !rdr.is_2d) return -1;

// std::cout << "Reading: " << infn << std::endl;

	size_t j0 = particles.size();
	
	int cur_id = -1;
	Particle *pptc = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			particles.resize(particles.size()+1);
			pptc = & particles[particles.size()-1];
			cur_id = id;
		}
		if (!pptc) continue;
		pptc->fill.push_back(HSeg(y, xl, xr));
	}
	
	for (size_t j=j0; j<particles.size(); j++) {
		Particle& ptc = particles[j];
		ptc.fill.shrink_to_fit();
		ptc.update_from_fill();
	}
	
	return cur_id;
}

//--------------------- Global data --------------------------------

NbrPoint hood_pts[37] = {
	{0, 0},

	{0, -1},
	{1, 0},
	{0, 1},
	{-1, 0},

	{-1, -1},
	{1, -1},
	{1, 1},
	{-1, 1},

	{-1, -2},
	{0, -2},
	{1, -2},
	{2, -1},
	{2, 0},
	{2, 1},
	{1, 2},
	{0, 2},
	{-1, 2},
	{-2, 1},
	{-2, 0},
	{-2, -1},

	{-2, -2},
	{2, -2},
	{2, 2},
	{-2, 2},
	{-1, -3},
	{0, -3},
	{1, -3},
	{3, -1},
	{3, 0},
	{3, 1},
	{1, 3},
	{0, 3},
	{-1, 3},
	{-3, 1},
	{-3, 0},
	{-3, -1}
};

NbrPoint3D hood3d_pts[27] = {
	{0, 0, 0},

	{1, 0, 0},
	{0, -1, 0},
	{-1, 0, 0},
	{0, 1, 0},
	{0, 0, 1},
	{0, 0, -1},

	{-1, -1, 0},
	{1, -1, 0},
	{1, 1, 0},
	{-1, 1, 0},
	{1, 0, 1},
	{0, -1, 1},
	{-1, 0, 1},
	{0, 1, 1},
	{1, 0, -1},
	{0, -1, -1},
	{-1, 0, -1},
	{0, 1, -1},

	{-1, -1, 1},
	{1, -1, 1},
	{1, 1, 1},
	{-1, 1, 1},
	{-1, -1, -1},
	{1, -1, -1},
	{1, 1, -1},
	{-1, 1, -1}
};

WalkAngle wangs[AUNITMAX] = {
	{ 1,  0.000000000000,  1,  1.000000000000},
	{ 1,  0.024930691738,  1,  0.999689182001},
	{ 1,  0.049845885661,  1,  0.998756921219},
	{ 1,  0.074730093586,  1,  0.997203797181},
	{ 1,  0.099567846596,  1,  0.995030775365},
	{ 1,  0.124343704647,  1,  0.992239206600},
	{ 1,  0.149042266176,  1,  0.988830826225},
	{ 1,  0.173648177667,  1,  0.984807753012},
	{ 1,  0.198146143199,  1,  0.980172487849},
	{ 1,  0.222520933956,  1,  0.974927912182},
	{ 1,  0.246757397690,  1,  0.969077286229},
	{ 1,  0.270840468143,  1,  0.962624246950},
	{ 1,  0.294755174411,  1,  0.955572805786},
	{ 1,  0.318486650252,  1,  0.947927346167},
	{ 1,  0.342020143326,  1,  0.939692620786},
	{ 1,  0.365341024366,  1,  0.930873748644},
	{ 1,  0.388434796275,  1,  0.921476211870},
	{ 1,  0.411287103131,  1,  0.911505852312},
	{ 1,  0.433883739118,  1,  0.900968867902},
	{ 1,  0.456210657353,  1,  0.889871808811},
	{ 1,  0.478253978621,  1,  0.878221573370},
	{ 1,  0.500000000000,  1,  0.866025403784},
	{ 1,  0.521435203379,  1,  0.853290881632},
	{ 1,  0.542546263866,  1,  0.840025923151},
	{ 1,  0.563320058064,  1,  0.826238774316},
	{ 1,  0.583743672235,  1,  0.811938005716},
	{ 1,  0.603804410325,  1,  0.797132507223},
	{ 1,  0.623489801859,  1,  0.781831482468},
	{ 1,  0.642787609687,  1,  0.766044443119},
	{ 1,  0.661685837597,  1,  0.749781202968},
	{ 1,  0.680172737771,  1,  0.733051871830},
	{ 1,  0.698236818086,  1,  0.715866849260},
	{ 1,  0.715866849260,  1,  0.698236818086},
	{ 1,  0.733051871830,  1,  0.680172737771},
	{ 1,  0.749781202968,  1,  0.661685837597},
	{ 1,  0.766044443119,  1,  0.642787609687},
	{ 1,  0.781831482468,  1,  0.623489801859},
	{ 1,  0.797132507223,  1,  0.603804410325},
	{ 1,  0.811938005716,  1,  0.583743672235},
	{ 1,  0.826238774316,  1,  0.563320058064},
	{ 1,  0.840025923151,  1,  0.542546263866},
	{ 1,  0.853290881632,  1,  0.521435203379},
	{ 1,  0.866025403784,  1,  0.500000000000},
	{ 1,  0.878221573370,  1,  0.478253978621},
	{ 1,  0.889871808811,  1,  0.456210657353},
	{ 1,  0.900968867902,  1,  0.433883739118},
	{ 1,  0.911505852312,  1,  0.411287103131},
	{ 1,  0.921476211870,  1,  0.388434796275},
	{ 1,  0.930873748644,  1,  0.365341024366},
	{ 1,  0.939692620786,  1,  0.342020143326},
	{ 1,  0.947927346167,  1,  0.318486650252},
	{ 1,  0.955572805786,  1,  0.294755174411},
	{ 1,  0.962624246950,  1,  0.270840468143},
	{ 1,  0.969077286229,  1,  0.246757397690},
	{ 1,  0.974927912182,  1,  0.222520933956},
	{ 1,  0.980172487849,  1,  0.198146143199},
	{ 1,  0.984807753012,  1,  0.173648177667},
	{ 1,  0.988830826225,  1,  0.149042266176},
	{ 1,  0.992239206600,  1,  0.124343704647},
	{ 1,  0.995030775365,  1,  0.099567846596},
	{ 1,  0.997203797181,  1,  0.074730093586},
	{ 1,  0.998756921219,  1,  0.049845885661},
	{ 1,  0.999689182001,  1,  0.024930691738},
	{ 1,  1.000000000000,  1,  0.000000000000},
	{ 1,  0.999689182001, -1,  0.024930691738},
	{ 1,  0.998756921219, -1,  0.049845885661},
	{ 1,  0.997203797181, -1,  0.074730093586},
	{ 1,  0.995030775365, -1,  0.099567846596},
	{ 1,  0.992239206600, -1,  0.124343704647},
	{ 1,  0.988830826225, -1,  0.149042266176},
	{ 1,  0.984807753012, -1,  0.173648177667},
	{ 1,  0.980172487849, -1,  0.198146143199},
	{ 1,  0.974927912182, -1,  0.222520933956},
	{ 1,  0.969077286229, -1,  0.246757397690},
	{ 1,  0.962624246950, -1,  0.270840468143},
	{ 1,  0.955572805786, -1,  0.294755174411},
	{ 1,  0.947927346167, -1,  0.318486650252},
	{ 1,  0.939692620786, -1,  0.342020143326},
	{ 1,  0.930873748644, -1,  0.365341024366},
	{ 1,  0.921476211870, -1,  0.388434796275},
	{ 1,  0.911505852312, -1,  0.411287103131},
	{ 1,  0.900968867902, -1,  0.433883739118},
	{ 1,  0.889871808811, -1,  0.456210657353},
	{ 1,  0.878221573370, -1,  0.478253978621},
	{ 1,  0.866025403784, -1,  0.500000000000},
	{ 1,  0.853290881632, -1,  0.521435203379},
	{ 1,  0.840025923151, -1,  0.542546263866},
	{ 1,  0.826238774316, -1,  0.563320058064},
	{ 1,  0.811938005716, -1,  0.583743672235},
	{ 1,  0.797132507223, -1,  0.603804410325},
	{ 1,  0.781831482468, -1,  0.623489801859},
	{ 1,  0.766044443119, -1,  0.642787609687},
	{ 1,  0.749781202968, -1,  0.661685837597},
	{ 1,  0.733051871830, -1,  0.680172737771},
	{ 1,  0.715866849260, -1,  0.698236818086},
	{ 1,  0.698236818086, -1,  0.715866849260},
	{ 1,  0.680172737771, -1,  0.733051871830},
	{ 1,  0.661685837597, -1,  0.749781202968},
	{ 1,  0.642787609687, -1,  0.766044443119},
	{ 1,  0.623489801859, -1,  0.781831482468},
	{ 1,  0.603804410325, -1,  0.797132507223},
	{ 1,  0.583743672235, -1,  0.811938005716},
	{ 1,  0.563320058064, -1,  0.826238774316},
	{ 1,  0.542546263866, -1,  0.840025923151},
	{ 1,  0.521435203379, -1,  0.853290881632},
	{ 1,  0.500000000000, -1,  0.866025403784},
	{ 1,  0.478253978621, -1,  0.878221573370},
	{ 1,  0.456210657353, -1,  0.889871808811},
	{ 1,  0.433883739118, -1,  0.900968867902},
	{ 1,  0.411287103131, -1,  0.911505852312},
	{ 1,  0.388434796275, -1,  0.921476211870},
	{ 1,  0.365341024366, -1,  0.930873748644},
	{ 1,  0.342020143326, -1,  0.939692620786},
	{ 1,  0.318486650252, -1,  0.947927346167},
	{ 1,  0.294755174411, -1,  0.955572805786},
	{ 1,  0.270840468143, -1,  0.962624246950},
	{ 1,  0.246757397690, -1,  0.969077286229},
	{ 1,  0.222520933956, -1,  0.974927912182},
	{ 1,  0.198146143199, -1,  0.980172487849},
	{ 1,  0.173648177667, -1,  0.984807753012},
	{ 1,  0.149042266176, -1,  0.988830826225},
	{ 1,  0.124343704647, -1,  0.992239206600},
	{ 1,  0.099567846596, -1,  0.995030775365},
	{ 1,  0.074730093586, -1,  0.997203797181},
	{ 1,  0.049845885661, -1,  0.998756921219},
	{ 1,  0.024930691738, -1,  0.999689182001},
	{ 1,  0.000000000000, -1,  1.000000000000},
	{-1,  0.024930691738, -1,  0.999689182001},
	{-1,  0.049845885661, -1,  0.998756921219},
	{-1,  0.074730093586, -1,  0.997203797181},
	{-1,  0.099567846596, -1,  0.995030775365},
	{-1,  0.124343704647, -1,  0.992239206600},
	{-1,  0.149042266176, -1,  0.988830826225},
	{-1,  0.173648177667, -1,  0.984807753012},
	{-1,  0.198146143199, -1,  0.980172487849},
	{-1,  0.222520933956, -1,  0.974927912182},
	{-1,  0.246757397690, -1,  0.969077286229},
	{-1,  0.270840468143, -1,  0.962624246950},
	{-1,  0.294755174411, -1,  0.955572805786},
	{-1,  0.318486650252, -1,  0.947927346167},
	{-1,  0.342020143326, -1,  0.939692620786},
	{-1,  0.365341024366, -1,  0.930873748644},
	{-1,  0.388434796275, -1,  0.921476211870},
	{-1,  0.411287103131, -1,  0.911505852312},
	{-1,  0.433883739118, -1,  0.900968867902},
	{-1,  0.456210657353, -1,  0.889871808811},
	{-1,  0.478253978621, -1,  0.878221573370},
	{-1,  0.500000000000, -1,  0.866025403784},
	{-1,  0.521435203379, -1,  0.853290881632},
	{-1,  0.542546263866, -1,  0.840025923151},
	{-1,  0.563320058064, -1,  0.826238774316},
	{-1,  0.583743672235, -1,  0.811938005716},
	{-1,  0.603804410325, -1,  0.797132507223},
	{-1,  0.623489801859, -1,  0.781831482468},
	{-1,  0.642787609687, -1,  0.766044443119},
	{-1,  0.661685837597, -1,  0.749781202968},
	{-1,  0.680172737771, -1,  0.733051871830},
	{-1,  0.698236818086, -1,  0.715866849260},
	{-1,  0.715866849260, -1,  0.698236818086},
	{-1,  0.733051871830, -1,  0.680172737771},
	{-1,  0.749781202968, -1,  0.661685837597},
	{-1,  0.766044443119, -1,  0.642787609687},
	{-1,  0.781831482468, -1,  0.623489801859},
	{-1,  0.797132507223, -1,  0.603804410325},
	{-1,  0.811938005716, -1,  0.583743672235},
	{-1,  0.826238774316, -1,  0.563320058064},
	{-1,  0.840025923151, -1,  0.542546263866},
	{-1,  0.853290881632, -1,  0.521435203379},
	{-1,  0.866025403784, -1,  0.500000000000},
	{-1,  0.878221573370, -1,  0.478253978621},
	{-1,  0.889871808811, -1,  0.456210657353},
	{-1,  0.900968867902, -1,  0.433883739118},
	{-1,  0.911505852312, -1,  0.411287103131},
	{-1,  0.921476211870, -1,  0.388434796275},
	{-1,  0.930873748644, -1,  0.365341024366},
	{-1,  0.939692620786, -1,  0.342020143326},
	{-1,  0.947927346167, -1,  0.318486650252},
	{-1,  0.955572805786, -1,  0.294755174411},
	{-1,  0.962624246950, -1,  0.270840468143},
	{-1,  0.969077286229, -1,  0.246757397690},
	{-1,  0.974927912182, -1,  0.222520933956},
	{-1,  0.980172487849, -1,  0.198146143199},
	{-1,  0.984807753012, -1,  0.173648177667},
	{-1,  0.988830826225, -1,  0.149042266176},
	{-1,  0.992239206600, -1,  0.124343704647},
	{-1,  0.995030775365, -1,  0.099567846596},
	{-1,  0.997203797181, -1,  0.074730093586},
	{-1,  0.998756921219, -1,  0.049845885661},
	{-1,  0.999689182001, -1,  0.024930691738},
	{-1,  1.000000000000, -1,  0.000000000000},
	{-1,  0.999689182001,  1,  0.024930691738},
	{-1,  0.998756921219,  1,  0.049845885661},
	{-1,  0.997203797181,  1,  0.074730093586},
	{-1,  0.995030775365,  1,  0.099567846596},
	{-1,  0.992239206600,  1,  0.124343704647},
	{-1,  0.988830826225,  1,  0.149042266176},
	{-1,  0.984807753012,  1,  0.173648177667},
	{-1,  0.980172487849,  1,  0.198146143199},
	{-1,  0.974927912182,  1,  0.222520933956},
	{-1,  0.969077286229,  1,  0.246757397690},
	{-1,  0.962624246950,  1,  0.270840468143},
	{-1,  0.955572805786,  1,  0.294755174411},
	{-1,  0.947927346167,  1,  0.318486650252},
	{-1,  0.939692620786,  1,  0.342020143326},
	{-1,  0.930873748644,  1,  0.365341024366},
	{-1,  0.921476211870,  1,  0.388434796275},
	{-1,  0.911505852312,  1,  0.411287103131},
	{-1,  0.900968867902,  1,  0.433883739118},
	{-1,  0.889871808811,  1,  0.456210657353},
	{-1,  0.878221573370,  1,  0.478253978621},
	{-1,  0.866025403784,  1,  0.500000000000},
	{-1,  0.853290881632,  1,  0.521435203379},
	{-1,  0.840025923151,  1,  0.542546263866},
	{-1,  0.826238774316,  1,  0.563320058064},
	{-1,  0.811938005716,  1,  0.583743672235},
	{-1,  0.797132507223,  1,  0.603804410325},
	{-1,  0.781831482468,  1,  0.623489801859},
	{-1,  0.766044443119,  1,  0.642787609687},
	{-1,  0.749781202968,  1,  0.661685837597},
	{-1,  0.733051871830,  1,  0.680172737771},
	{-1,  0.715866849260,  1,  0.698236818086},
	{-1,  0.698236818086,  1,  0.715866849260},
	{-1,  0.680172737771,  1,  0.733051871830},
	{-1,  0.661685837597,  1,  0.749781202968},
	{-1,  0.642787609687,  1,  0.766044443119},
	{-1,  0.623489801859,  1,  0.781831482468},
	{-1,  0.603804410325,  1,  0.797132507223},
	{-1,  0.583743672235,  1,  0.811938005716},
	{-1,  0.563320058064,  1,  0.826238774316},
	{-1,  0.542546263866,  1,  0.840025923151},
	{-1,  0.521435203379,  1,  0.853290881632},
	{-1,  0.500000000000,  1,  0.866025403784},
	{-1,  0.478253978621,  1,  0.878221573370},
	{-1,  0.456210657353,  1,  0.889871808811},
	{-1,  0.433883739118,  1,  0.900968867902},
	{-1,  0.411287103131,  1,  0.911505852312},
	{-1,  0.388434796275,  1,  0.921476211870},
	{-1,  0.365341024366,  1,  0.930873748644},
	{-1,  0.342020143326,  1,  0.939692620786},
	{-1,  0.318486650252,  1,  0.947927346167},
	{-1,  0.294755174411,  1,  0.955572805786},
	{-1,  0.270840468143,  1,  0.962624246950},
	{-1,  0.246757397690,  1,  0.969077286229},
	{-1,  0.222520933956,  1,  0.974927912182},
	{-1,  0.198146143199,  1,  0.980172487849},
	{-1,  0.173648177667,  1,  0.984807753012},
	{-1,  0.149042266176,  1,  0.988830826225},
	{-1,  0.124343704647,  1,  0.992239206600},
	{-1,  0.099567846596,  1,  0.995030775365},
	{-1,  0.074730093586,  1,  0.997203797181},
	{-1,  0.049845885661,  1,  0.998756921219},
	{-1,  0.024930691738,  1,  0.999689182001}
};

