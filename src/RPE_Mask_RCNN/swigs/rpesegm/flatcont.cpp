
#include "flatcont.h"

Boundary boundary_from_path(std::vector<Point> & path)
{
	Boundary bnd;
	if (path.size() > 0) {
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
	return bnd;
}

double out_of_boundary(Boundary & bnd, const std::vector<Point> & path, double *p_len)
{
	if (path.size() < 2) {
		if (p_len) *p_len = 0.;
		return 0.;
	}
	double t_dist = 0., o_dist = 0.;
	for (size_t j=1; j<path.size()-1; j++) {
		Point p1 = path[j-1];
		Point p2 = path[j];
		double d = p1.dist(p2);
		t_dist += d;
		int inside = 0;
		if (bnd.IsInside(p1)) ++inside;
		if (bnd.IsInside(p2)) ++inside;
		if (inside == 1) {
			o_dist += d * 0.5;
		} else if (inside == 0) {
			o_dist += d;
		}
	}
	if (p_len) *p_len = t_dist;
	return o_dist;
}

/*
// Simplify path by removing pixels in the middle of straight lines
static void simplify_path_0(std::vector<Point> &path)
{
	std::vector<Point> spath;
	spath.push_back(path[0]);
	for (size_t j=1; j<path.size()-1; j++) {
		Point p = path[j];
		Point pl = path[j-1];
		Point pr = path[j+1];
		if (p.x-pl.x != pr.x-p.x || p.y-pl.y != pr.y-p.y)
			spath.push_back(p);
	}
	spath.push_back(path[path.size()-1]);
	path.resize(spath.size());
	for (size_t j=0; j<spath.size(); j++)
		path[j] = spath[j];
}
*/

// 0 - good to go
// 1 - bad, but there is hope
// 2 - too bad, no hope
static int same_segment(std::vector<Point> &path, size_t i0, size_t i1)
{
	if (i0+2 > i1) return 1;
	Point p0 = path[i0];
	Point p1 = path[i1];
	int dx = p0.x - p1.x;
	int dy = p0.y - p1.y;
	if (dx == 0 && dy == 0) return 2;
	if (dx < 0) dx = -dx;
	if (dy < 0) dy = -dy;
	int nbad = 0;
	if (dx >= dy) {
		if (p1.x < p0.x) p0.swap(p1);
		dy = p1.y - p0.y;
		for (size_t j=i0+1; j<i1; j++) {
			Point p2 = path[j];
			int y2 = p0.y + ((p2.x - p0.x) * dy) / dx;
			if (y2 != p2.y) ++nbad;
		}
	} else {
		if (p1.y < p0.y) p0.swap(p1);
		dx = p1.x - p0.x;
		for (size_t j=i0+1; j<i1; j++) {
			Point p2 = path[j];
			int x2 = p0.x + ((p2.y - p0.y) * dx) / dy;
			if (x2 != p2.x) ++nbad;
		}
	}
	if (!nbad) return 0;
	double df = p0.dist(p1);
	double dl = 0.;
	for (size_t i=i0+1; i<=i1; i++) {
		dl += path[i-1].dist(path[i]);
	}
	return (dl-df > 2.) ? 2 : 1;
}

void simplify_path(std::vector<Point> &path)
{
	std::vector<Point> spath;
	size_t i0 = 0;
	while (i0 < path.size()) {
		spath.push_back(path[i0]);
		if (i0 + 2 >= path.size()) {
			++i0;
			continue;
		}
		size_t best_i = 0;
		for (size_t i1=i0+2; i1<path.size(); i1++) {
			int rc = same_segment(path, i0, i1);
			if (rc == 2) break;
			if (!rc) best_i = i1;
		}
		if (best_i > i0) {
			i0 = best_i;
		} else {
			++i0;
		}
	}
	path.resize(spath.size());
	for (size_t j=0; j<spath.size(); j++)
		path[j] = spath[j];
}

static bool detect_contour(Raster8 & msk, Contour& cont, int x0, int y0, unsigned char bordc, bool simplify=true)
{
	unsigned char fgd = msk.value(x0, y0);
	unsigned char tmpbord = fgd + 1;
	std::vector<Point> path;
	std::vector<Point> savepts;
	
	while (true) {
		// std::cout << x0 << "," << y0 << " path=" << path.size() << std::endl;
		if (path.empty()) {
			// Find first border pixel around (x0, y0)
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				int x = x0 + hood_pts[j].dx;
				int y = y0 + hood_pts[j].dy;
				if (msk.value(x, y) == bordc) {
					path.push_back(Point(x,y));
					savepts.push_back(Point(x,y));
					msk.setValue(x, y, tmpbord);
					break;
				}
			}
		} else if (path.size() == 1) {
			// Find second border pixel next to first
			Point pcur = path[path.size()-1];
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				int x = pcur.x + hood_pts[j].dx;
				int y = pcur.y + hood_pts[j].dy;
				if (msk.value(x, y) == bordc && msk.touches(x, y, fgd)) {
					path.push_back(Point(x,y));
					savepts.push_back(Point(x,y));
					msk.setValue(x, y, tmpbord);
					break;
				}
			}
			if (path.size() == 1) {
				path.clear();
			}
		} else {
			// Find next border pixel close to the current (-1) but as far as possible from the previous (-2)
			Point pcur = path[path.size()-1];
			Point prev = path[path.size()-1];
			int xn = -1, yn = -1;
			double cdist = 0.;
			for (int j=1; j<HOOD_SIZE_MOORE; j++) {
				int x = pcur.x + hood_pts[j].dx;
				int y = pcur.y + hood_pts[j].dy;
				if (msk.value(x, y) == bordc && msk.touches(x, y, fgd)) {
					Point pnext(x, y);
					double d = pnext.dist(prev);
					if (d > cdist) {
						cdist = d;
						xn = x;
						yn = y;
					}
				}
			}
			if (xn >= 0 && yn >= 0) {
				path.push_back(Point(xn, yn));
				savepts.push_back(Point(xn, yn));
				msk.setValue(xn, yn, tmpbord);
				// Check if close enough to path start (0)
				if (path.size() > 10 && path[0].dist(xn, yn) < 2) break;
			} else {
				// No matches -- go to the previous point and try again
				path.resize(path.size() - 1);
			}
		}
		
		if (path.empty()) break;
	}
	
	for (Point pt : savepts) {
		msk.setValue(pt.x, pt.y, bordc);
	}
	Boundary bnd = boundary_from_path(path);
	if (path.empty()) {
		bnd.xmin = bnd.xmax = x0;
		bnd.ymin = bnd.ymax = y0;
	}
	bnd.expand(3);
	msk.clip(bnd, 1);
	msk.replaceColor(bnd, tmpbord, bordc);
	if (path.size() < 5) return false;
	
	if (simplify)
		simplify_path(path);
	cont.set_path(path);
	return true;
}

bool detect_particle_contour(Raster8 & msk, Contour& cont, Particle &ptc, unsigned char bordc, bool simplify)
{
	int peri = 0;
	for (HSeg & hs : ptc.fill) {
		for (int x=hs.xl; x<=hs.xr; x++) {
			if (msk.touches(x, hs.y, bordc)) ++peri;
		}
	}
	if (peri < 10) return false;
	for (HSeg & hs : ptc.fill) {
		for (int x=hs.xl; x<=hs.xr; x++) {
			if (!msk.touches(x, hs.y, bordc, HOOD_SIZE_MOORE))
				continue;
			if (!detect_contour(msk, cont, x, hs.y, bordc, simplify))
				continue;
			if (cont.plen > peri/2) return true;
		}
	}
	return false;
}

std::vector<int> particles_to_contours(Raster8 &msk, std::vector<Particle> & particles)
{
	std::vector<int> res;

	msk.fill(0x80);
	for (Particle& ptc : particles)
		msk.paintParticle(ptc, 0);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0);
	msk.fillBorder(0x10, 1);
	
	std::vector<Contour> contours;

	Boundary bnd = msk.getBoundary();
	bnd.expand(-BORDER_SIZE);

	int tot = 0;
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *b = msk.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (b[x0] != 0) continue;
			Particle ptc;
			// std::cout << "particle at " << x0 << "," << y0 << std::endl;
			msk.detectParticle(ptc, x0, y0, 0x50);
			contours.resize(contours.size()+1);
			Contour& cont = contours[contours.size()-1];
			double pct = 100.;
			if (detect_contour(msk, cont, x0, y0, 0xFF, false)) {
				++tot;
				double peri;
				double oob = out_of_boundary(bnd, cont.path, &peri);
				if (peri > 0.001) pct = 100.*oob/peri;
			}
			msk.paintParticle(ptc, 0x40);
			if (pct >= AT_BORDER_IMPORT_PCT) {
				contours.resize(contours.size()-1);
			} else {
				simplify_path(cont.path);
			}
		}
	}
	
	for (Contour& cont : contours) {
		res.push_back(int(cont.path.size()));
		for (Point &pt : cont.path) {
			res.push_back(pt.x);
			res.push_back(pt.y);
		}
	}
	
	// std::cout << "particles_to_contours(): " << contours.size() << " out of " << tot << std::endl;
	
	return res;
}

std::vector<CellContour> unflatten_cell_contours(const std::vector<int>& flat_cells)
{
	std::vector<CellContour> res;
	size_t iStart = 0;
	int cur_z = -1;
	int cur_id = -1;
	int next_idx = 0;
	while (iStart < flat_cells.size()) {
		int z = flat_cells[iStart++];
		int id = flat_cells[iStart++];
		int len = flat_cells[iStart++];
		Contour cont(flat_cells, int(iStart), len);
		iStart += (len + len);
		if (cur_z != z || cur_id != id) {
			cur_z = z;
			cur_id = id;
			next_idx = 0;
		}
		res.resize(res.size()+1);
		CellContour& cc = res[res.size()-1];
		cc.cell_id = id;
		cc.z = z;
		cc.cont_idx = next_idx++;
		cc.cont.set_path(cont.path);
		cc.valid = cc.cont.path.size() >= 5;
	}
	return res;
}

