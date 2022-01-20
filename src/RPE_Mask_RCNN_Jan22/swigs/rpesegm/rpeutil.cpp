
#include "raster.h"
#include "flatcont.h"
#include <set>
#include "csv.h"
#include "rpeutil.h"

const double DNA_FILL = 0.6;

// Output is converted to Python list of integers formatted as:
// <n1>(<xi><yi>.. repeated n1 times)<n2> (<xi><yi>..repeated n2 times) ...
std::vector<int> import_contours(int w, int h, const char *in_csv)
{
	std::vector<Particle> particles;
	int rc = read_particle_data(in_csv, particles, w, h);
	// std::cout << "import_contours(): " << rc << " particles read from " << in_csv << std::endl;
	
	Raster8 msk(w, h, NULL);
	return particles_to_contours(msk, particles);
}

std::vector<int> import_reshape_contours(unsigned char *mask, int hm, int wm, const char *rs_csv)
{
	Raster8 msk(wm, hm, mask);
	Boundary mbnd = msk.getBoundary();
	std::vector<Particle> particles;
	
	CsvReader rdr(rs_csv);
	if (rdr.next()) {
		CsvRow headers = rdr.GetRow();
		rdr.setHeaders(headers);

		int xs_idx = rdr.header_index("XStart");
		int ys_idx = rdr.header_index("YStart");
		while (rdr.next()) {
			int x0 = rdr.GetInt(xs_idx, -1);
			int y0 = rdr.GetInt(ys_idx, -1);
			if (!mbnd.IsInside(x0, y0)) continue;
			particles.resize(particles.size()+1);
			Particle & ptc = particles[particles.size()-1];
			msk.detectParticle(ptc, x0, y0, 0x80);
		}
	}
	// std::cout << "import_reshape_contours(): " << particles.size() << " particles read from " << rs_csv << std::endl;
	return particles_to_contours(msk, particles);
}

struct ZContour
{
	int z;
	Contour cont;
};

// Output is converted to Python list of integers formatted as:
// <z1><id1><n1>(<xi><yi>..repeated n1 times) <z2><id2><n2>(<xi><yi>..repeated n2 times) ...
// id-s are 0-based; result is sorted by (z,id)
// Contours with the same id belong to the same cell (there can be more than 1 per z-frame).
std::vector<int> import_3d_contours(int w, int h, int d, const char *in_csv)
{
	std::vector<int> res;

	std::vector<Particle3D> cells;
	int rc = read_cell_data(in_csv, cells, w, h, d);
	// std::cout << "import_3d_contours(): " << rc << " cells read from " << in_csv << std::endl;

	Raster8 msk(w, h, NULL);
	std::vector<std::vector<ZContour>> cell_conts(cells.size());

	for (int z=0; z<d; z++) {
		msk.fill(0x80);
		for (Particle3D& cell : cells) {
			if (!cell.fills[z].empty())
				msk.paintParticleFill(cell.fills[z], 0);
		}
		msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0);
		msk.fillBorder(0x10, 1);

		for (int idx=0; size_t(idx)<cells.size(); idx++) {
			std::vector<HSeg> & fill = cells[idx].fills[z];
			if (fill.empty()) continue;
			std::vector<ZContour> & conts = cell_conts[idx];
			
			for (HSeg& hs : fill) {
				int y0 = hs.y;
				unsigned char *p = msk.scanLine(y0);
				for (int x0=hs.xl; x0<=hs.xr; x0++) {
					if (p[x0] != 0) continue;
					
					Particle ptc;
					msk.detectParticle(ptc, x0, y0, 0x50);
					
					conts.resize(conts.size()+1);
					ZContour & zc = conts[conts.size()-1];
					zc.z = z;
					if (!detect_particle_contour(msk, zc.cont, ptc, 0xFF, false)) {
						conts.resize(conts.size()-1);
					}
					msk.paintParticle(ptc, 0x40);
				}
			}
		}
	}
	
	Boundary bnd = msk.getBoundary();
	bnd.expand(-BORDER_SIZE);
	int next_id = 0;
	for (int idx=0; size_t(idx)<cell_conts.size(); idx++) {
		std::vector<ZContour> & conts = cell_conts[idx];
		if (conts.empty()) continue;
		double sum_peri = 0;
		double out_peri = 0.;
		for (ZContour &zc : conts) {
			double peri;
			double oob = out_of_boundary(bnd, zc.cont.path, &peri);
			sum_peri += peri;
			out_peri += oob;
		}
		if (sum_peri <= 0.001) continue;
		if (100.*out_peri/sum_peri >= AT_BORDER_IMPORT_PCT) continue;
		for (ZContour &zc : conts) {
			simplify_path(zc.cont.path);
			res.push_back(zc.z);
			res.push_back(next_id);
			res.push_back(int(zc.cont.path.size()));
			for (Point &pt : zc.cont.path) {
				res.push_back(pt.x);
				res.push_back(pt.y);
			}
		}
		++next_id;
	}
	
	// std::cout << "Imported " << next_id << " cells out of " << cells.size() << std::endl;
	
	return res;
}

static std::vector<int> find_split_slices(std::vector<CellContour> & cell_conts,
		Slice & ptc, int cell_id, int z0, int d)
{
	std::vector<int> res;
	std::vector<int> to_split;
	std::vector<int> to_keep;
	
	for (int j=0; size_t(j)<cell_conts.size(); j++) {
		CellContour& cc = cell_conts[j];
		if (!cc.valid || cc.cell_id != cell_id || cc.z < z0-2 || cc.z > z0+2) continue;
		int ovl = ptc.overlay_area(cc.ptc);
		double ratio = double(ovl + ovl) / double(ptc.area + cc.ptc.area);
		cc.valid = false;
		if (ratio > 0.85)
			to_split.push_back(j);
		else if (ratio < 0.25)
			to_keep.push_back(j);
		else cc.valid = true;
	}
	
	if (to_split.size() < 2 || to_keep.size() < 1)
		return res;
	
	for (int dz=0; dz<d; dz++) {
		if (z0-dz < 0 && z0+dz >= d) break;
		for (int j=0; size_t(j)<cell_conts.size(); j++) {
			CellContour& cc = cell_conts[j];
			if (!cc.valid || cc.cell_id != cell_id || (cc.z != z0-dz && cc.z != z0+dz)) continue;
			cc.valid = false;
			
			int ovl = 0;
			double best_split = 0.;
			double split_pct = 0.;
			int zdif = z0 - cc.z;
			if (zdif >= -5 && zdif <= 5) {
				ovl = ptc.overlay_area(cc.ptc);
				best_split = double(ovl + ovl) / double(ptc.area + cc.ptc.area);
				split_pct = double(ovl) / double(cc.ptc.area + 1);
			}
			
			double best_keep = 0.;
			int split_ovl = 0;
			int keep_ovl = 0;
			for (int js=0; size_t(js)<to_split.size(); js++) {
				CellContour& ccs = cell_conts[to_split[js]];
				int zdif = ccs.z - cc.z;
				if (zdif < -5 || zdif > 5) continue;
				ovl = ccs.ptc.overlay_area(cc.ptc);
				if (ovl == 0) continue;
				double ratio = double(ovl + ovl) / double(ccs.ptc.area + cc.ptc.area);
				if (ratio > best_split) {
					split_pct = double(ovl) / double(cc.ptc.area + 1);
					best_split = ratio;
					split_ovl = ovl;
				}
			}
			for (int jk=0; size_t(jk)<to_keep.size(); jk++) {
				CellContour& cck = cell_conts[to_keep[jk]];
				int zdif = cck.z - cc.z;
				if (zdif < -5 || zdif > 5) continue;
				int ovl = cck.ptc.overlay_area(cc.ptc);
				if (ovl == 0) continue;
				double ratio = double(ovl + ovl) / double(cck.ptc.area + cc.ptc.area + 1);
				if (ratio > best_keep) {
					best_keep = ratio;
					keep_ovl = ovl;
				}
			}
			if (best_split > 0.9 || best_keep < 0.25 || best_split > best_keep + 0.1) {
				to_split.push_back(j);
			}
			else
				to_keep.push_back(j);
		}
	}
	
	for (int j : to_split) {
		CellContour& cc = cell_conts[j];
		res.push_back(cc.z);
		res.push_back(cc.cont_idx);
	}
	
	return res;
}

static std::vector<int> find_join_slices(std::vector<CellContour> & cell_conts,
		Slice & ptc, int cell_id, int z0, int d)
{
	std::vector<int> res;
	
	std::vector<int> to_join;
	std::set<int> join_ids;
	std::set<int> keep_ids;
	
	join_ids.insert(cell_id);
	
	for (int j=0; size_t(j)<cell_conts.size(); j++) {
		CellContour& cc = cell_conts[j];
		if (!cc.valid || cc.z < z0-2 || cc.z > z0+2) continue;
		if (cc.cell_id == cell_id) {
			to_join.push_back(j);
			cc.valid = false;
		}
		int ovl = ptc.overlay_area(cc.ptc);
		double ratio = double(ovl + ovl) / double(ptc.area + cc.ptc.area);
		if (ratio > 0.85) {
			join_ids.insert(cc.cell_id);
		} else if (ratio < 0.25) {
			keep_ids.insert(cc.cell_id);
		}
	}
	
	// std::cout << "cell ID: " << cell_id << std::endl;
	for (int idk : keep_ids) {
		if (idk == cell_id) {
			join_ids.clear();
			break;
		}
		join_ids.erase(idk);
	}
	if (join_ids.size() < 2) return res;
	
	for (int dz=0; dz<d; dz++) {
		if (z0-dz < 0 && z0+dz >= d) break;
		
		for (int j=0; size_t(j)<cell_conts.size(); j++) {
			CellContour& cc = cell_conts[j];
			if (!cc.valid || (cc.z != z0-dz && cc.z != z0+dz)) continue;
			cc.valid = false;
			if (join_ids.find(cc.cell_id) == join_ids.end()) continue;
			
			int ovl = 0;
			double best_join = 0.;
			double join_pct = 0.;
			int zdif = z0 - cc.z;
			if (zdif >= -5 && zdif <= 5) {
				ovl = ptc.overlay_area(cc.ptc);
				best_join = double(ovl + ovl) / double(ptc.area + cc.ptc.area);
				join_pct = double(ovl) / double(cc.ptc.area + 1);
			}
			int join_ovl = 0;
			for (int jj=0; size_t(jj)<to_join.size(); jj++) {
				CellContour& ccj = cell_conts[to_join[jj]];
				int zdif = ccj.z - cc.z;
				if (zdif < -5 || zdif > 5) continue;
				ovl = ccj.ptc.overlay_area(cc.ptc);
				double ratio = double(ovl + ovl) / double(ccj.ptc.area + cc.ptc.area);
				if (ratio > best_join) {
					join_pct = double(ovl) / double(cc.ptc.area + 1);
					best_join = ratio;
					join_ovl = ovl;
				}
			}
			
			if (best_join > 0.6667 && join_pct > 0.75) {
				to_join.push_back(j);
			} else if (best_join < 0.3333) {
				join_ids.erase(cc.cell_id);
				if (join_ids.size() < 2) return res;
			}
		}
	}

	for (int idj : join_ids) {
		if (keep_ids.find(idj) != keep_ids.end() || idj == cell_id) continue;
		// std::cout << "OK to join ID " << idj << " -> " << cell_id << std::endl;
		res.push_back(cell_id);
		res.push_back(idj);
	}

	return res;
}

std::vector<int> cell_for_contour(int w, int h, int d, int z0,
		const std::vector<int>& flat_cells, const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	
	std::vector<CellContour> cell_conts = unflatten_cell_contours(flat_cells);

	Raster8 msk(w, h, NULL);
	msk.fill(0x80);
	
	for (CellContour& cc : cell_conts) {
		if (!cc.valid) continue;
		cc.cont.bnd.expand(1);
		msk.clip(cc.cont.bnd, 2);
		msk.fillContour(cc.cont, 0xA0, 0xA0);
		cc.ptc.bnd = cc.cont.bnd;
		cc.ptc.area = msk.rescanParticle(cc.ptc, 0xA0);
		msk.paintParticle(cc.ptc, 0x80);
	}
	
	Contour cont(flat_cont);
	simplify_path(cont.path);
	cont._update_from_path();

	Slice ptc;
	msk.clip(ptc.bnd, 2);

	msk.fillContour(cont, 0x50, 0x50);
	
	for (CellContour& cc : cell_conts) {
		if (!cc.valid || cc.z != z0) continue;
		msk.paintParticle(cc.ptc, 0xC0);
	}
	
	msk.expandBorders(0xC0, 0xC1, HOOD_SIZE_MOORE, 0x50);
	msk.replaceColor(0xC0, 0x80);
	msk.replaceColor(0xC1, 0x80);
	msk.fillBorder(0x80, 2);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0x50);
	
	ptc.bnd = cont.bnd;
	ptc.area = msk.rescanParticle(ptc, 0x50);
	// std::cout << "cell_for_contour(): ptc.area = " << ptc.area << std::endl;
	
	if (ptc.area < 20 || !detect_particle_contour(msk, cont, ptc, 0xFF) ||
			cont.path.size() < 5) {
		// std::cout << "detect_particle_contour() failed!" << std::endl;
		res.push_back(-1);
		res.push_back(0);
		return res;
	}
	
	int min_ovl_area = int(0.666667 * ptc.area); // at least 2/3 area overlay
	
	int zs = z0 - 2; if (zs < 0) zs = 0;
	int ze = z0 + 2; if (ze >= d) ze = d - 1;
	
	int best_id = -1;
	double best_ratio = 0.;
	int best_z = -1;
	
	for (CellContour& cc : cell_conts) {
		if (!cc.valid || cc.z < zs || cc.z > ze) continue;
		int ovl = ptc.overlay_area(cc.ptc);
		if (ovl < min_ovl_area) continue;
		double ratio = double(ovl + ovl) / double(ptc.area + cc.ptc.area);
		if (ratio > best_ratio) {
			best_ratio = ratio;
			best_id = cc.cell_id;
			best_z = cc.z;
		}
	}
	
	if (cont.path.size() >= 3) {
		if (cont.path[0].equals(cont.path[cont.path.size()-1]))
			cont.path.resize(cont.path.size()-1);
	}

	res.push_back(best_id);
	res.push_back(int(cont.path.size()));
	for (Point &pt : cont.path) {
		res.push_back(pt.x);
		res.push_back(pt.y);
	}
	
	// std::cout << "best_id=" << best_id << " best_ratio=" << best_ratio << std::endl;

	if (best_ratio > 0.85) {
		// Slices to split from best_id
		std::vector<int> to_split = find_split_slices(cell_conts, ptc, best_id, z0, d);
		res.push_back(int(to_split.size()/2));
		for (int v : to_split)
			res.push_back(v);
		if (to_split.empty()) {
			// Cell IDs to join best_id
			std::vector<int> to_join = find_join_slices(cell_conts, ptc, best_id, z0, d);
			res.push_back(int(to_join.size()/2));
			for (int v : to_join)
				res.push_back(v);
		} else {
			res.push_back(0);
		}
	}
	else {
		res.push_back(0);
		res.push_back(0);
	}

	return res;
}

std::vector<int> validate_contour(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	Contour cont(flat_cont);
	cont.optimize(1.5, 2.);
	simplify_path(cont.path);
	cont._update_from_path();

	Raster8 msk(w, h, NULL);
	msk.fill(0x80);
	
	Particle ptc;
	msk.clip(ptc.bnd, 2);
	
	msk.fillContour(cont, 0x50, 0x50);
	ptc.bnd = cont.bnd;
	int o_area = msk.rescanParticle(ptc, 0x50);

	size_t i = 0;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour _cont(flat_slices, int(i), len);
		i += (len + len);
		if (len > 3 && cont.bnd.intersects(_cont.bnd)) {
			msk.clip(_cont.bnd, 2);
			msk.fillContour(_cont, 0xC0, 0xC0);
		}
	}

	msk.expandBorders(0xC0, 0xC1, HOOD_SIZE_MOORE, 0x50);
	msk.replaceColor(0xC0, 0x80);
	msk.replaceColor(0xC1, 0x80);
	msk.fillBorder(0x80, 2);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0x50);
	
	ptc.bnd = cont.bnd;
	int c_area = msk.rescanParticle(ptc, 0x50);
	
	if (c_area >= 20) {
		if (o_area > c_area + 10) {
			if (!detect_particle_contour(msk, cont, ptc, 0xFF) || cont.path.size() < 5)
				c_area = 0;
		}
	} else c_area = 0;
	if (c_area <= 0) {
		res.push_back(0);
		return res;
	}
	
	if (cont.path[0].equals(cont.path[cont.path.size()-1]))
		cont.path.resize(cont.path.size()-1);

	res.push_back(int(cont.path.size()));
	for (Point &pt : cont.path) {
		res.push_back(pt.x);
		res.push_back(pt.y);
	}

	return res;
}

std::vector<int> simplify_contour(const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	Contour cont(flat_cont);
	simplify_path(cont.path);

	res.push_back(int(cont.path.size()));
	for (Point &pt : cont.path) {
		res.push_back(pt.x);
		res.push_back(pt.y);
	}

	return res;
}

std::vector<int> intersecting_slices(const std::vector<int>& flat_slices, const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	
	Contour cont(flat_cont);
	if (cont.path.size() < 5) return res;
	
	size_t i = 0;
	int idx = -1;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour _cont(flat_slices, int(i), len);
		i += (len + len);
		++idx;
		//if (len < 3) continue;
		bool found = false;
		for (Point pt : _cont.path) {
			if (cont.IsInside(pt)) {
				found = true;
				res.push_back(idx);
				break;
			}
		}
		if (found) continue;
		for (Point pt : cont.path) {
			if (_cont.IsInside(pt)) {
				found = true;
				res.push_back(idx);
				break;
			}
		}
	}

	return res;
}

// Output:
// indexes to delete: <nidx><idx0, idx1, ... (repeated nidx times)>,
// followed by flat contours to add: <n1><(x0,y0), (x1,y1)...n1 times> <n2><(x0,y0), (x1,y1)...n2 times> ...
std::vector<int> cut_slices(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	Contour cut(flat_cont);
	if (cut.is_closed) {
		// Un-"close" contour if the ends are far away
		size_t last = cut.path.size() - 1;
		if (cut.path[last-1].dist(cut.path[last]) > cut.plen*0.1) {
			cut.path.resize(cut.path.size()-1);
			cut.is_closed = false;
			cut._update_from_path();
		}
	}
	Raster8 msk(w, h, NULL);
	msk.fill(0x80);
	// msk.paintContourBorder(cut, 0xC0);
	
	res.push_back(0);

	int cont_idx = -1;
	std::vector<Slice> slices;
	size_t i = 0;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour cont(flat_slices, int(i), len);
		i += (len + len);
		++cont_idx;
		if (!cont.is_closed) continue;
		cont.bnd.expand(1);
		msk.clip(cont.bnd, 2);
		msk.fillContour(cont, 0x50, 0x50);
		msk.paintContourBorder(cut, 0xC0);
		
		size_t o_size = slices.size();
		for (int y0=cont.bnd.ymin; y0<=cont.bnd.ymax; y0++) {
			unsigned char *p = msk.scanLine(y0);
			for (int x0=cont.bnd.xmin; x0<=cont.bnd.xmax; x0++) {
				if (p[x0] != 0x50) continue;
				slices.resize(slices.size() + 1);
				Slice& ptc = slices[slices.size()-1];
				ptc.area = msk.detectParticle(ptc, x0, y0, 0);
			}
		}
		if (slices.size() == o_size+1) {
			// No split
			slices.resize(slices.size()-1);
		}
		if (slices.size() > o_size) {
			res.push_back(cont_idx);
			++res[0];
		}
	}
	
	if (slices.empty()) {
		res.clear();
		res.push_back(0);
		return res;
	}
	
	msk.fill(0x80);
	for (Slice& ptc : slices)
		msk.paintParticle(ptc, 0);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0);
	msk.fillBorder(0x10, 1);
	
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *b = msk.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (b[x0] != 0) continue;
			Particle ptc;
			int ar = msk.detectParticle(ptc, x0, y0, 0x50);
			Contour cont;
			if (detect_particle_contour(msk, cont, ptc, 0xFF)) {
				res.push_back(int(cont.path.size()));
				for (Point &pt : cont.path) {
					res.push_back(pt.x);
					res.push_back(pt.y);
				}
			}
			msk.paintParticle(ptc, 0x40);
		}
	}

	return res;
}

// Output:
// indexes to delete: <nidx><idx0, idx1, ... (repeated nidx times)>,
// followed by flat contours to add: <n1><(x0,y0), (x1,y1)...n1 times> <n2><(x0,y0), (x1,y1)...n2 times> ...
std::vector<int> join_slices(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont)
{
	std::vector<int> res;
	Contour jcont(flat_cont);
	Raster8 msk(w, h, NULL);
	msk.fill(0x80);
	
	std::vector<Slice> particles;
	
	size_t i = 0;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour cont(flat_slices, int(i), len);
		i += (len + len);
		cont.bnd.expand(1);
		msk.clip(cont.bnd, 2);
		msk.fillContour(cont, 0x50, 0x50);
		particles.resize(particles.size() + 1);
		Slice &ptc = particles[particles.size() - 1];
		ptc.bnd = cont.bnd;
		ptc.area = msk.rescanParticle(ptc, 0x50);
		msk.fill(ptc.bnd, 0x80);
	}

	jcont.bnd.expand(1);
	msk.clip(jcont.bnd, 2);
	msk.fillContour(jcont, 0x50, 0x50);
	
	res.push_back(0);
	Boundary bnd = jcont.bnd;
	for (int idx=0; size_t(idx)<particles.size(); idx++) {
		Slice & ptc = particles[idx];
		if (msk.countColor(ptc, 0x50) == 0) continue;
		msk.paintParticle(ptc, 0x51);
		bnd.combo(ptc.bnd);
		res.push_back(idx);
		++res[0];
	}
	msk.replaceColor(bnd, 0x51, 0x50);

	msk.expandBorders(0x50, 0x51, HOOD_SIZE_NEUMANN, 0x80);
	msk.expandBorders(0x80, 0x52, HOOD_SIZE_MOORE, 0x51);
	msk.replaceColor(0x51, 0x50);
	msk.replaceColor(0x52, 0x80);
	
	Slice jptc;
	jptc.bnd = bnd;
	msk.rescanParticle(jptc, 0x50);
	
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_NEUMANN, 0x50);
	msk.fillBorder(0x10, 1);
	
	if (!detect_particle_contour(msk, jcont, jptc, 0xFF)) {
		res.clear();
		res.push_back(0);
	} else {
		res.push_back(int(jcont.path.size()));
		for (Point &pt : jcont.path) {
			res.push_back(pt.x);
			res.push_back(pt.y);
		}
	}
	
	return res;
}

// 'flat_cells' parameter is expected to be a Python list of integers formatted as:
// <z1><id1><n1>(<xi><yi>..repeated n1 times) <z2><id2><n2>(<xi><yi>..repeated n2 times) ...
// id-s are 0-based; result is sorted by (z,id)
// Contours with the same id belong to the same cell (there can be more than 1 per z-frame).
// [This is the same format as output from import_3d_contours()]
void export_3d_contours(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const std::vector<int> & flat_cells, const char *out_csv, bool validate, bool separate)
{
	std::vector<CellContour> cell_conts = unflatten_cell_contours(flat_cells);
	// std::cout << "export_3d_contours(): Loaded " << cell_conts.size() << " contours" << std::endl;
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	
	Raster8 msk = mstack.getPlane(0);
	msk.fill(0x80);
	
	std::vector<Cell> cells;
	int cur_cell_id = -1;
	int cur_z = -1;
	Cell *pcell = NULL;
	Boundary cur_bnd;
	for (CellContour& cc : cell_conts) {
		if (!cc.valid) continue;
// std::cout << "cc cell_id=" << cc.cell_id << " z=" << cc.z << " cont_idx=" << cc.cont_idx << std::endl;
		if ((cur_z != cc.z || cur_cell_id != cc.cell_id) && pcell != NULL && cur_z >= 0) {
			msk.rescanParticleFill(cur_bnd, pcell->fills[cur_z], 0x50);
			msk.fill(cur_bnd, 0x80);
		}
		if (cc.cell_id != cur_cell_id) {
			int idx = int(cells.size());
			cells.resize(cells.size()+1);
			pcell = & cells[idx];
			pcell->fills.resize(mstack.d);
			pcell->idx = cc.cell_id;
			cur_cell_id = idx;
			cur_z = -1;
		}
		cc.cont.bnd.expand(1);
		msk.clip(cc.cont.bnd, 2);
		msk.fillContour(cc.cont, 0x50, 0x50);
		if (cur_z != cc.z) {
			cur_bnd = cc.cont.bnd;
			cur_z = cc.z;
		} else {
			cur_bnd.combo(cc.cont.bnd);
		}
	}
	if (pcell) {
		msk.rescanParticleFill(cur_bnd, pcell->fills[cur_z], 0x50);
	}
	
	// std::cout << "export_3d_contours(): decoded " << cells.size() << " cells" << std::endl;
	
	// Make sure no slice overlaps, paint cells of 0 on 0x80 background
	mstack.fill(0x80);
	for (int z=0; z<mstack.d; z++) {
		Raster8 _msk = mstack.getPlane(z);
		for (Cell &cell : cells) {
			if (cell.fills[z].empty()) continue;
			Boundary bnd = fill_boundary(cell.fills[z]);
			_msk.paintParticleFillInto(cell.fills[z], 0x50, 0x80);
			_msk.rescanParticleFill(bnd, cell.fills[z], 0x50);
			_msk.paintParticleFill(cell.fills[z], 0);
		}
	}
	for (Cell &cell : cells) {
		cell.update_from_fill();
	}
	
	// See if some cells need to be merged
	if (validate)
	for (size_t i=1; i<cells.size(); i++) {
		Cell & other = cells[i];
		Boundary o_bnd = other.bnd.boundary2d();
		for (size_t j=0; j<i; j++) {
			Cell & cell = cells[j];
			if (cell.idx < 0) continue;
			if (!cell.bnd.intersects2d(o_bnd)) continue;
			int z, o_z;
			if (cell.bnd.zmax+1 == other.bnd.zmin) {
				z = cell.bnd.zmax;
				o_z = other.bnd.zmin;
			} else if (other.bnd.zmax+1 == cell.bnd.zmin) {
				z = cell.bnd.zmin;
				o_z = other.bnd.zmax;
			} else continue;
			int ar = fill_area(cell.fills[z]);
			int o_ar = fill_area(other.fills[o_z]);
			int ovl = fill_overlay_area(cell.fills[z], other.fills[o_z]);
			double pct = double(ovl) / double(ar + 1);
			double o_pct = double(ovl) / double(o_ar + 1);
			if ((pct > 0.85 && o_pct > 0.6) || (pct > 0.6 && o_pct > 0.85)) {
				for (z=other.bnd.zmin; z<=other.bnd.zmax; z++) {
					cell.fills[z] = std::move(other.fills[z]);
				}
				other.clear();
				other.idx = -1;
				cell.update_from_fill();
				break;
			}
		}
	}

	// See if some cells need to be split
	if (validate)
	for (size_t j=0; j<cells.size(); j++) {
		Cell * pcell = & cells[j];
		if (pcell->idx < 0) continue;
		int ztop = pcell->bnd.zmax;
		int zbot = pcell->bnd.zmin;
		int z = ztop;
		while (z > zbot) {
			if (pcell->fills[z-1].empty() || fill_overlay_area(pcell->fills[z], pcell->fills[z-1]) < 5) {
//std::cout << "Split cell ID=" << pcell->idx << " at z=" << z << ".." << ztop << " x0=" << pcell->bnd.xmin << " y0=" << pcell->bnd.ymin << std::endl;
				
				int idx = int(cells.size());
				cells.resize(cells.size()+1);
				Cell * ncell = & cells[idx];
				ncell->fills.resize(mstack.d);
				ncell->idx = idx;
				pcell = & cells[j];		// Needs to be updated after cells.resize()

				for (int zz=z; zz<=ztop; zz++)
					ncell->fills[zz] = std::move(pcell->fills[zz]);
				ncell->update_from_fill();
				pcell->update_from_fill();
//std::cout << "New cell ID=" << ncell->idx << " z=" << ncell->bnd.zmin << ".." << ncell->bnd.zmax << std::endl;
//std::cout << "Remaining z=" << pcell->bnd.zmin << ".." << pcell->bnd.zmax << std::endl;

				for (ztop=z-1; ztop>=zbot; ztop--)
					if (!pcell->fills[ztop].empty()) break;
				z = ztop;
			} else --z;
		}
	}

	if (separate) {
		// std::cout << "force separation" << std::endl;
		for (Cell &cell : cells) {
			if (cell.idx < 0) continue;
			mstack.forceSeparation(cell, 0, 0x80, 0x50);
		}
	}

	mstack.expandBorders(0x80, 0, 0xFF);
	
	std::cout << "Write: " << out_csv << std::endl;
	write_cell_data(cells, out_csv);
}

void export_2d_contours(unsigned char *mask, int hm, int wm,
		const std::vector<int>& flat_slices, const char *out_csv)
{
	Raster8 msk(wm, hm, mask);
	msk.fill(0x80);
	
	std::vector<Slice> particles;
	
	size_t i = 0;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour cont(flat_slices, int(i), len);
		i += (len + len);
		if (len < 5) continue;
		cont.bnd.expand(1);
		msk.clip(cont.bnd, 2);
		msk.fillContour(cont, 0x50, 0x50);
		particles.resize(particles.size() + 1);
		Slice &ptc = particles[particles.size() - 1];
		ptc.bnd = cont.bnd;
		ptc.area = msk.rescanParticle(ptc, 0x50);
		msk.fill(ptc.bnd, 0x80);
	}

	for (Slice & ptc : particles) {
		msk.paintParticleInto(ptc, 0x50, 0x80);
		msk.rescanParticle(ptc, 0x50);
		msk.paintParticle(ptc, 0);
	}
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
	
	std::cout << "Write: " << out_csv << std::endl;
	write_particle_data(particles, out_csv);
}

std::vector<int> cells_at_border(int w, int h, const std::vector<int> & flat_cells)
{
	std::vector<int> res;
	Boundary bnd(BORDER_SIZE, BORDER_SIZE, w-1-BORDER_SIZE, h-1-BORDER_SIZE);

	std::vector<CellContour> cell_conts = unflatten_cell_contours(flat_cells);
	
	int cur_cell_id = -1;
	double sum_peri = 0.;
	double out_peri = 0.;
	for (CellContour& cc : cell_conts) {
		if (cc.cell_id != cur_cell_id) {
			if (sum_peri > 0.001) {
				double pct = 100.*out_peri/sum_peri;
				if (pct >= AT_BORDER_PCT) res.push_back(cur_cell_id);
			}
			sum_peri = out_peri = 0.;
		}
		cur_cell_id = cc.cell_id;
		if (!cc.valid) continue;
		double peri;
		double oob = out_of_boundary(bnd, cc.cont.path, &peri);
		sum_peri += peri;
		out_peri += oob;
	}
	if (sum_peri > 0.001) {
		double pct = 100.*out_peri/sum_peri;
		if (pct >= AT_BORDER_PCT) res.push_back(cur_cell_id);
	}
	return res;
}

std::vector<int> contours_at_border(int w, int h, const std::vector<int> & flat_slices)
{
	std::vector<int> res;
	Boundary bnd(BORDER_SIZE, BORDER_SIZE, w-1-BORDER_SIZE, h-1-BORDER_SIZE);

	int idx = 0;
	size_t i = 0;
	while (i < flat_slices.size()) {
		int len = flat_slices[i++];
		Contour cont(flat_slices, int(i), len);
		i += (len + len);

		double peri;
		double oob = out_of_boundary(bnd, cont.path, &peri);
		if (peri > 0.001) {
			double pct = 100.*oob/peri;
			if (pct >= AT_BORDER_PCT) res.push_back(idx);
		}
		++idx;
	}
	
	return res;
}


static void filterParticlesIn(Raster8& msk, Boundary& bnd, unsigned char oldc, unsigned char newc, int minsz, unsigned char bk)
{
	for (int y0=bnd.ymin; y0<=bnd.ymax; y0++) {
		for (int x0=bnd.xmin; x0<=bnd.xmax; x0++) {
			if (msk.value(x0, y0) != oldc) continue;
			std::vector<HSeg> fill;
			msk.findParticleFill(fill, x0, y0, newc);
			if (fill_area(fill) < minsz)
				msk.paintParticleFill(fill, bk);
		}
	}
}


// roi = [ymin xmin ymax+1 xmax+1]
std::vector<int> mask_rcnn_to_particles(unsigned char *ptmask, int npts, int hptm, int wptm,
		int *rois, int nrois, int roisz, int x_orig, int y_orig)
{
	std::vector<int> res;
	
	size_t fr_size = size_t(hptm) * size_t(wptm);
	size_t roi_size = size_t(roisz);
	
	for (int ipt=0; ipt<npts; ipt++) {
		Raster8 msk(wptm, hptm, ptmask + (ipt * fr_size));
		int *roi = rois + (ipt * roi_size);
		Slice ptc;
		ptc.bnd = Boundary(roi[1], roi[0], roi[3], roi[2]);
		ptc.bnd.expand(1);
		msk.clip(ptc.bnd, 1);
		
		filterParticlesIn(msk, ptc.bnd, 0xFF, 0xA0, 100, 0x80);
		msk.bordersAround(ptc.bnd, 0xA0, 0, 0x80, HOOD_SIZE_NEUMANN);
		msk.bordersAround(ptc.bnd, 0, 0x80, 0x40, HOOD_SIZE_MOORE);
		msk.replaceColor(ptc.bnd, 0x40, 0);
		msk.replaceColor(ptc.bnd, 0x80, 0xA0);
		
		ptc.area = msk.rescanParticle(ptc, 0xA0);
		
		res.push_back(int(ptc.fill.size()));
		for (HSeg & hs : ptc.fill) {
			res.push_back(hs.y+y_orig);
			res.push_back(hs.xl+x_orig);
			res.push_back(hs.xr+x_orig);
		}
	}
	
	return res;
}

static void fix_borders_dna(Raster8& msk, std::vector<Slice>& particles)
{
	msk.fillBorder(0x10, 1);
	for (Slice& ptc : particles) {
		if (ptc.area <= 0) continue;
		// Filter out "phantom slices"
		double ff = double(msk.countColor(ptc, 0xFF)) / double(ptc.area);
		if (ff < DNA_FILL) {
			ptc.area = -1;
		} else {
			msk.paintParticle(ptc, 0xC0);
		}
	}
	for (int pass=0; pass<2; pass++) {
		int nbsz = (pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
		for (Slice& ptc : particles) {
			if (ptc.area <= 0) continue;
			Boundary bnd = ptc.bnd;
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.paintParticle(ptc, 0xD0);
			msk.bordersAround(bnd, 0xD0, 0xFF, 0xA0, nbsz);
			
			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				unsigned char *p = msk.scanLine(y);
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					if (p[x] != 0xA0) continue;
					unsigned char xc = 0xD0;
					for (int j=1; j<HOOD_SIZE_MOORE; j++) {
						unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
						if (c == 0xC0) {
							xc = 0;
							break;
						}
					}
					p[x] = xc;
				}
			}

			ptc.bnd = bnd;
			ptc.area = msk.rescanParticle(ptc, 0xD0);
			msk.paintParticle(ptc, 0xC0);
		}
	}
}

static void fix_borders_actin(Raster8& msk, std::vector<Slice>& particles)
{
	msk.fill(0);
	for (Slice& ptc : particles) {
		if (ptc.area <= 0) continue;
		msk.paintParticle(ptc, 0xA0);
		ptc.bnd.expand(1);
		msk.clip(ptc.bnd, 1);
		
		filterParticlesIn(msk, ptc.bnd, 0xA0, 0xC0, 100, 0x80);
		msk.replaceColor(ptc.bnd, 0x80, 0);
		
		ptc.area = msk.rescanParticle(ptc, 0xC0);
	}
	msk.fillBorder(0x10, 1);
	
//	msk.filterParticles(0xA0, 0xC0, 200, 0x80);
//	msk.replaceColor(0x80, 0);
//	for (Slice& ptc : particles) {
//		if (ptc.area <= 0) continue;
//		msk.paintParticleInto(ptc, 0xA0, 0xC0);
//		ptc.area = msk.rescanParticle(ptc, 0xA0);
//		msk.paintParticle(ptc, 0xC0);
//	}
	
	msk.expandBorders(0xC0, 0x5, HOOD_SIZE_MOORE, 0);
	msk.expandBorders(0, 0x15, HOOD_SIZE_MOORE, 0x5);
	msk.replaceColor(0x15, 0);
	msk.replaceColor(0x5, 0x40);
	
	msk.expandBorders(0xC0, 0x40, HOOD_SIZE_FATCROSS, 0);
	msk.expandBorders(0x40, 0x30, HOOD_SIZE_FATCROSS, 0);
	msk.replaceColor(0x30, 0x40);
	msk.filterParticles(0, 0x5, 300, 0x80);
	msk.replaceColor(0x5, 0);
	msk.replaceColor(0x80, 0x40);
	msk.expandBorders(0, 0x60, HOOD_SIZE_FATCROSS, 0x40);
	msk.replaceColor(0x60, 0);
	msk.expandBorders(0, 0x60, HOOD_SIZE_NEUMANN, 0x40);
	msk.replaceColor(0x60, 0);

	for (int pass=0; pass<8; pass++) {
		int nbsz = (pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
		for (Slice& ptc : particles) {
			if (ptc.area <= 0) continue;
			Boundary bnd = ptc.bnd;
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.paintParticle(ptc, 0xD0);
			msk.bordersAround(bnd, 0xD0, 0x40, 0xA0, nbsz);
			
			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				unsigned char *p = msk.scanLine(y);
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					if (p[x] != 0xA0) continue;
					unsigned char xc = 0xD0;
					for (int j=1; j<HOOD_SIZE_MOORE; j++) {
						unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
						if (c == 0xC0) {
							xc = 0;
							break;
						}
					}
					p[x] = xc;
				}
			}

			ptc.bnd = bnd;
			ptc.area = msk.rescanParticle(ptc, 0xD0);
			msk.paintParticle(ptc, 0xC0);
		}
	}
}

std::vector<int> recombine_flat_particles(unsigned char *mask, int hm, int wm,
		std::vector<int> flat_particles, std::vector<double> scores, int fix_borders)
{
	std::vector<Slice> particles;
	size_t iStart = 0;
	while (iStart < flat_particles.size()) {
		int hlen = flat_particles[iStart++];
		particles.resize(particles.size()+1);
		Slice& ptc = particles[particles.size()-1];
		for (int i=0; i<hlen; i++) {
			ptc.fill.push_back(HSeg(flat_particles[iStart], flat_particles[iStart+1], flat_particles[iStart+2]));
			iStart += 3;
		}
		ptc.area = ptc.update_from_fill();
		if (ptc.area < 50) {
			ptc.area = 0;
		}
	}
	
	// std::cout << "Decoded: " << particles.size() << " objects; scores: " << scores.size() << std::endl;
	
	// De-dupe
	for (size_t j=0; j<particles.size(); j++) {
		if (j >= scores.size()) break;
		Slice& ptc0 = particles[j];
		if (ptc0.area <= 0) continue;
		int maxovl0 = ptc0.area / 2;
		for (size_t i=0; i<j; i++) {
			Slice& ptc = particles[i];
			if (ptc.area <= 0) continue;
			int ovl = ptc.overlay_area(ptc0);
			if (ovl <= maxovl0 && ovl <= ptc.area / 2) continue;
			if (scores[j] <= scores[i]) {
				ptc0.area = -1;
				break;
			} else {
				ptc.area = -1;
			}
		}
	}
	
	Raster8 msk(wm, hm, mask);
	
	if (fix_borders == FIX_BORDERS_DNA) {
		fix_borders_dna(msk, particles);
	} else if (fix_borders == FIX_BORDERS_ACTIN) {
		fix_borders_actin(msk, particles);
	}
	
	msk.fill(0);
	msk.fillBorder(0x10, 1);
	
	for (size_t j=0; j<particles.size(); j++) {
		Slice& ptc = particles[j];
		if (ptc.area <= 0) continue;
		Boundary bnd = ptc.bnd;
		bnd.expand(1);
		msk.clip(bnd, 1);
		
		msk.paintParticleInto(ptc, 0xA0, 0);
		
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != 0xA0) continue;
				bool has_nbr = false;
				for (int j=1; j<HOOD_SIZE_NEUMANN; j++) {
					unsigned char c = msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy);
					if (c == 0xA0) {
						has_nbr = true;
						break;
					}
				}
				if (!has_nbr) p[x] = 0;
			}
		}
		
		msk.bordersAround(bnd, 0xA0, 0, 0x40, HOOD_SIZE_MOORE);
		
		ptc.bnd = bnd;
		ptc.area = msk.rescanParticle(ptc, 0xA0);
		if (ptc.area <= 0) continue;
		msk.paintParticle(ptc, 0xFF);
		
	}

	msk.fillBorder(0, 1);
	msk.replaceColor(0x40, 0);
	
	int n_scores = 0;
	std::vector<int> res;
	res.push_back(n_scores);
	for (size_t j=0; j<particles.size(); j++) {
		Slice& ptc = particles[j];
		if (ptc.area <= 0) continue;

		if (j < scores.size()) {
			double cx, cy;
			fill_centroid(ptc.fill, &cx, &cy);
			res.push_back(int(cx));
			res.push_back(int(cy));
			res.push_back(int(scores[j] * 1000.));
			++n_scores;
		}
	}
	
	res[0] = n_scores;
	return res;
}

void mask2to3(unsigned char *mask, int hm, int wm)
{
	Raster8 msk(wm, hm, mask);
	msk.replaceColor(0, 0x80);
	msk.replaceColor(0xFF, 0);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
}

void test_fix_borders_actin(unsigned char *mask, int hm, int wm)
{
	Raster8 msk(wm, hm, mask);
	msk.fillBorder(0x10, 1);
	std::vector<Slice> particles;
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (p[x0] != 0xFF) continue;
			particles.resize(particles.size()+1);
			Slice& ptc = particles[particles.size()-1];
			ptc.area = msk.detectParticle(ptc, x0, y0, 0xA0);
		}
	}
	std::cout << "Detected " << particles.size() << " cells." << std::endl;
	fix_borders_actin(msk, particles);
}

