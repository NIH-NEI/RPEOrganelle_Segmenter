
#include "raster.h"
#include "assembly3d.h"
#include "ml3d.h"
#include "csv.h"
#include "rpesegm.h"

const double DNA_FILL = 0.6;

const double BASIC_OVERLAY = 0.8;
const double BASIC_FULL_OVERLAY = 0.9;
const double BASIC_OVL_RATIO = 0.75;

static void assemble_basic_cell(Particle3D& cell,
		std::vector<std::vector<AsmSlice>> &all_slices, int z0, int orig_idx)
{
	int d = int(all_slices.size());
	for (int z=z0+2; z<d; z++) {
		int z1 = z - 1;
		if (cell.fills[z1].empty()) --z1;
		if (cell.fills[z1].empty()) break;
		std::vector<HSeg>& my_fill = cell.fills[z1];
		int my_area = fill_area(my_fill);
		Boundary my_bnd = fill_boundary(my_fill);
		for (AsmSlice& ptc : all_slices[z]) {
			if (ptc.taken()) continue;
			if (!my_bnd.intersects(ptc.bnd)) continue;
			int ovl = ptc.overlay_area(my_fill);
			if (ovl < 1) continue;
			if (ovl < int(ptc.area * BASIC_FULL_OVERLAY)) {
				double ovl_ratio = (double(ovl) + ovl) / (double(ptc.area) + my_area);
				if (ovl_ratio < BASIC_OVL_RATIO) continue;
			}
			cell.add_slice(ptc, z);
			ptc.orig_idx = orig_idx;
		}
	}
	for (int z=z0-1; z>=0; z--) {
		int z1 = z + 1;
		if (cell.fills[z1].empty()) ++z1;
		if (cell.fills[z1].empty()) break;
		std::vector<HSeg>& my_fill = cell.fills[z1];
		int my_area = fill_area(my_fill);
		Boundary my_bnd = fill_boundary(my_fill);
		for (AsmSlice& ptc : all_slices[z]) {
			if (ptc.taken()) continue;
			if (!my_bnd.intersects(ptc.bnd)) continue;
			int ovl = ptc.overlay_area(my_fill);
			if (ovl < 1) continue;
			if (ovl < int(ptc.area * BASIC_FULL_OVERLAY)) {
				double ovl_ratio = (double(ovl) + ovl) / (double(ptc.area) + my_area);
				if (ovl_ratio < BASIC_OVL_RATIO) continue;
			}
			cell.add_slice(ptc, z);
			ptc.orig_idx = orig_idx;
		}
	}
}

static void missing_slice(std::vector<Particle3D>& cells, Raster8& msk, size_t idx, int z)
{
	Particle3D& cell = cells[idx];
	msk.fill(0);
	int zdn, zup;
	for (zdn=z-1; zdn>=cell.bnd.zmin; zdn--)
		if (!cell.fills[zdn].empty()) break;
	for (zup=z+1; zup<=cell.bnd.zmax; zup++)
		if (!cell.fills[zup].empty()) break;
	msk.paintParticleFill(cell.fills[zdn], 1);
	msk.paintParticleFill(cell.fills[zup], 1);
	msk.paintParticleFillADD(cell.fills[zdn], 1);
	msk.paintParticleFillADD(cell.fills[zup], 1);
	msk.fillBorder(0, 1);
	Boundary bnd = cell.bnd.boundary2d();
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = msk.scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			unsigned char c = p[x];
			if (c < 2) c = 0;
			else if (c == 2) {
				int cnt = 0;
				for (int j=1; j<HOOD_SIZE_MOORE; j++) {
					cnt += int(msk.value(x+hood_pts[j].dx, y+hood_pts[j].dy));
				}
				c = (cnt >= 13) ? 3 : 0;
			} else c = 3;
			p[x] = c;
		}
	}
	msk.rescanParticleFill(bnd, cell.fills[z], 3);
	
	msk.fill(bnd, 0);
	
	Boundary bnda = bnd;
	std::vector<size_t> crop_cells;
	for (size_t i=0; i<cells.size(); i++) {
		if (i == idx) continue;
		Particle3D& pt = cells[i];
		if (pt.fills[z].empty()) continue;
		Boundary b = pt.bnd.boundary2d();
		if (!b.intersects(bnd)) continue;
		crop_cells.push_back(i);
		msk.paintParticleFill(pt.fills[z], 0x80);
		bnda.combo(b);
	}
	bnda.expand(1);
	msk.clip(bnda, 2);
	msk.expandBorders(0, 0x40, HOOD_SIZE_FATCROSS, 0x80);
	
	msk.paintParticleFillInto(cell.fills[z], 0xC0, 0);
	msk.paintParticleFillInto(cell.fills[z], 0xC0, 0x40);
	msk.replaceColor(bnda, 0x40, 0x80);
	msk.expandBorders(0x80, 0x20, HOOD_SIZE_FATCROSS, 0xC0);
	msk.rescanParticleFill(bnd, cell.fills[z], 0xC0);
	
	for (size_t i : crop_cells) {
		Particle3D& pt = cells[i];
		Boundary b = pt.bnd.boundary2d();
		msk.paintParticleFillInto(pt.fills[z], 0x60, 0x80);
		msk.rescanParticleFill(b, pt.fills[z], 0x60);
		msk.paintParticleFill(pt.fills[z], 0x80);
	}
}

static void missing_slice_2(std::vector<Particle3D>& cells, Raster8& msk, size_t idx, int z)
{
	Particle3D& cell = cells[idx];
	msk.fill(0);
	int zdn, zup;
	for (zdn=z-1; zdn>=cell.bnd.zmin; zdn--)
		if (!cell.fills[zdn].empty()) break;
	for (zup=z+1; zup<=cell.bnd.zmax; zup++)
		if (!cell.fills[zup].empty()) break;
	msk.paintParticleFillADD(cell.fills[zdn], 1);
	msk.paintParticleFillADD(cell.fills[zup], 1);
	msk.fillBorder(0, 1);
	Boundary bnd = cell.bnd.boundary2d();
	msk.expandBordersInto(bnd, 2, 1, 4);
	msk.replaceColor(bnd, 1, 0);

	for (size_t i=0; i<cells.size(); i++) {
		if (i == idx) continue;
		Particle3D& pt = cells[i];
		if (pt.fills[z].empty()) continue;
		Boundary b = pt.bnd.boundary2d();
		if (!b.intersects(bnd)) continue;
		msk.paintParticleFill(pt.fills[z], 0x20);
	}
	bnd.expand(1);
	msk.clip(bnd, 1);
	msk.expandBordersInto(bnd, 0x20, 2, 4);
	
	msk.rescanParticleFill(bnd, cell.fills[z], 2);
}

void assemble_basic_body(Raster3D &mstack,
		std::vector<std::vector<AsmSlice>>& all_slices,
		const char *csvfile, int postproc)
{
	std::vector<Particle3D> cells;

	for (int z0=mstack.d-2; z0>=0; z0--) {
		int z = z0 + 1;
		std::vector<AsmSlice>& slices0 = all_slices[z0];
		std::vector<AsmSlice>& slices = all_slices[z];
		for (AsmSlice& ptc0 : slices0) {
			if (ptc0.taken()) continue;
			int min_ovl0 = int(BASIC_OVERLAY * ptc0.area);
			for (AsmSlice& ptc : slices) {
				if (ptc.taken()) continue;
				int ovl = ptc0.overlay_area(ptc);
				if (ovl < min_ovl0 || ovl < int(BASIC_OVERLAY * ptc.area))
					continue;
				int orig_idx = int(cells.size());
				cells.resize(cells.size()+1);
				Particle3D& cell = cells[orig_idx];
				cell.fills.resize(mstack.d);
				cell.bnd.zmin = z0;
				cell.bnd.zmax = z;
				cell.bnd.set2d(ptc0.bnd);
				cell.add_slice(ptc, z);
				cell.add_slice(ptc0, z0);
				ptc0.orig_idx = ptc.orig_idx = orig_idx;
				assemble_basic_cell(cell, all_slices, z0, orig_idx);
				break;
			}
		}
	}
	
	//std::cout << "Assembled " << cells.size() << " cells." << std::endl;
	
	if (postproc >= 1) {
		//mstack.fill(0);
		for (Particle3D& cell : cells) {
			cell.update_from_fill();
			//mstack.paintParticle(cell, 0x80);
		}
		
		Raster8 msk = Raster8(mstack.w, mstack.h, NULL);
		
		int nc = 0;
		int ns = 0;
		for (size_t idx=0; idx<cells.size(); idx++) {
			Particle3D& cell = cells[idx];
			int nempty = 0;
			for (int z=cell.bnd.zmin; z<=cell.bnd.zmax; z++) {
				if (cell.fills[z].empty()) {
					++nempty;
					missing_slice(cells, msk, idx, z);
					//Raster8 m = mstack.getPlane(z);
					//m.paintParticleFill(cell.fills[z], 0xFF);
				}
			}
			if (nempty > 0) {
				++nc;
				ns += nempty;
			}
		}
		//std::cout << "Empty slices: " << ns << "; Cells w. ES: " << nc << std::endl;
		//return;
	}
	
	mstack.fill(0x80);
	for (Particle3D& pt : cells) {
		pt.update_from_fill();
		mstack.forceBorder(pt, 0, 0xFF, 0x20, 0xE0);
	}

	std::cout << "Write object data to " << csvfile << std::endl;
	write_cell_data(cells, csvfile);
}

void assemble_basic(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, int postproc)
{
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);

	std::vector<std::vector<AsmSlice>> all_slices(mstack.d);
	
	size_t n_slices = 0;
	for (int z=0; z<mstack.d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		AsmSlice::detect_all(all_slices[z], msk, 0xFF, 0xA0, z);
		n_slices += all_slices[z].size();
	}
	
	assemble_basic_body(mstack, all_slices, csvfile, postproc);
}

void assemble_ml(std::vector<std::vector<int>> all_flat_particles,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile, int postproc)
{
	AssemblerML aml(mask3d, zm3d, hm3d, wm3d);
	aml.GOOD_IOU = 0.60;
	aml.ACCEPTABLE_IOU = 0.60;
	
	aml.load_flat_slices(all_flat_particles);
	
	aml.find_seeds();
	aml.consolidate_seeds();
	aml.cells_from_seeds();
	
	if (postproc == POSTPROC_DNA) {
		// std::cout << "fix_borders_dna" << std::endl;
		aml.fix_borders_dna();
	} else if (postproc == POSTPROC_ACTIN) {
		// std::cout << "fix_borders_actin" << std::endl;
		aml.fix_borders_actin();
	}
	aml.fill_gaps();
	
	aml.paint_final();
	
	std::cout << "Write " << aml.cells.size() << " objects to " << csvfile << std::endl;
	write_cell_data(aml.cells, csvfile);
}

void segment_dna(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	int otsu1, int otsu2)
{
	Raster16 img(w, h, data);
	Raster8 ang(w, h, angle);
	Raster8 msk(w, h, mask);
	
	Raster8 out(w, h, NULL);
	int bsz = 3;

	out.fill(0);
	
	// Detect particle contours based on strong gradient criteria
	int zcount = HOOD_SIZE_RAD3 * 33 / 100;
	DnaGradientWalker wk(img, ang, out, otsu1);
	wk.set_zero_mask(&msk, zcount);
	wk.setStep(4.);
	wk.bend = -7;
	wk.awiggle = 21;
	wk.maxgap = 10.;
	
	std::vector<Contour> closed_contours;
	std::vector<Contour> open_contours;

	for (int y0=3; y0<h-3; y0++) {
		unsigned char *mp = out.scanLine(y0);
		for (int x0=3; x0<w-3; x0++) {
			int v = img.value(x0, y0);
			if (v > otsu1 &&
					img.local_rank(x0, y0, HOOD_SIZE_FATCROSS) < 2 &&
					msk.countColors(x0, y0, HOOD_SIZE_RAD3, 0) > zcount) {
				std::vector<Point> path;
				int l = wk.walk(path, x0, y0);
				if (wk.in_loop) {
					path.push_back(path[0]);
					closed_contours.push_back(Contour(path));
					out.paintPath(path, 0xFF, 0x80);
				} else if (l > 20) {
					open_contours.push_back(Contour(path));
					out.paintPath(path, 0xFF, 0x80);
				}
			}
		
		}
	}
	
	out.fill(0);

	for (Contour &cont : closed_contours) {
		cont.optimize();
		out.fillContour(cont, 0xE0, 0xE0);
	}
	
	for (Contour &cont : open_contours) {
		Contour mcont(cont.path);
		if (mcont.can_close(mcont.plen*0.2) && wk.close_path(mcont.path)) {
			mcont.is_closed = true;
			mcont.optimize();
			closed_contours.push_back(Contour(mcont.path));
			out.fillContour(closed_contours[closed_contours.size()-1], 0xC0, 0xE0);
			cont.path.resize(0);
		}
	}
	
	// Try to close up some open contours using weaker gradient criteria
	wk.zcount = HOOD_SIZE_RAD3 * 2 / 100;
	wk.setStep(4.);
	wk.bend = -7;
	wk.awiggle = 21;
	wk.maxgap = 10.;
	wk.touch_angle = 94;

	for (Contour &cont : open_contours) {
		if (cont.path.size() == 0) continue;
		if (size_t(out.countColorsNb(cont.path, 0, HOOD_SIZE_MOORE) * 120/100) < cont.path.size()) {
			cont.path.resize(0);
			continue;
		}
		Point p0 = cont.path[0];
		wk.in_loop = false;
		wk.at_border = 0;
		Contour mcont(cont.path);
		
		wk.walk_one(mcont.path, AUNIT90, wk.bend);
		reverse_path(mcont.path);
		wk.walk_one(mcont.path, AUNIT270, -wk.bend);
		reverse_path(mcont.path);
		if (wk.in_loop) {
			mcont.path.push_back(mcont.path[0]);
			mcont.is_closed = true;
			mcont.optimize();
			closed_contours.push_back(Contour(mcont.path));
			out.fillContour(closed_contours[closed_contours.size()-1], 0xA0, 0xC0);
			cont.path.resize(0);
		} else if (mcont.can_close(mcont.plen*0.2) && wk.close_path(mcont.path)) {
			mcont.is_closed = true;
			mcont.optimize();
			closed_contours.push_back(Contour(mcont.path));
			out.fillContour(closed_contours[closed_contours.size()-1], 0x60, 0xA0);
			cont.path.resize(0);
		}

	}
	
	// Convert closed contours to particles
	out.fill(0);
	std::vector<Particle> particles;
	for (size_t j=0; j<closed_contours.size(); j++) {
		particles.resize(particles.size()+1);
		Particle & ptc = particles[particles.size()-1];
		
		ptc.fromContour(closed_contours[j]);
		double ff = double(msk.countColor(ptc, 0xFF)) / double(fill_area(ptc.fill));
		if (ff < DNA_FILL) {
			particles.resize(particles.size()-1);
		} else {
			out.paintParticle(ptc, 0xFF);
		}
	}
	
	// Create a 2-pix border inside old particles
	for (Particle &ptc : particles) {
		out.shrinkParticle(ptc, 0xFF, 0xF0,	0x80, HOOD_SIZE_FATCROSS);
	}
	
	for (int y=0; y<h; y++) {
		unsigned char *ob = out.scanLine(y);
		unsigned char *mb = msk.scanLine(y);
		for (int x=0; x<w; x++) {
			if (ob[x] != 0) continue;
			if (mb[x] == 0xFF) ob[x] = 0xA0;
		}
	}
	
	// Filter out black holes <300px inside new particles
	out.fillBorder(0x10, 2);
	out.filterParticles(0, 0x20, 300, 0xA0);
	out.replaceColor(0x20, 0);
	// Create a 2-pix border inside new particles
	out.expandBorders(0, 0x40, HOOD_SIZE_FATCROSS, 0xA0);
	out.expandBorders(0x80, 0x40, HOOD_SIZE_FATCROSS, 0xA0);
	out.replaceColor(0x40, 0x80);

	// Dilate new particle borders
	unsigned char curc = 0x80;
	unsigned char nextc = 0xB0;
	for (int pass=0; pass<8; pass++) {
		out.filterParticles(0xA0, 0x20, 500, 0xF0);
		out.replaceColor(0x20, 0xA0);
		out.expandBorders(curc, nextc, HOOD_SIZE_FATCROSS, 0xA0);
		curc = nextc;
		nextc += 4;
	}
	
	out.replaceColor(0xF0, 0xA0);	

	size_t old_np = particles.size();

	// Now detect new particles
	size_t np = old_np;
	particles.resize(np + 1);
	for (int y=0; y<h; y++) {
		unsigned char *b = out.scanLine(y);
		for (int x=0; x<w; x++) {
			if (b[x] != 0xA0) continue;
			Particle &ptc = particles[np];
			int a = out.detectParticle(ptc, x, y, 0xE0);
			if (a < 10) {
				out.paintParticle(ptc, 0x80);
			} else {
				out.paintParticle(ptc, 0xAA);
				++np;
				particles.resize(np + 1);
			}
		}
	}
	particles.resize(np);
	
	// Erode new particle borders
	// Alternate Von Neumann and Moore neighborhoods to reduce direction bias
	nextc = curc - 4;
	for (int pass=0; pass<8; pass++) {
		for (int k=0; k<2; k++) {
			for (size_t j=old_np; j<np; j++) {
				Particle &ptc = particles[j];
				out.expandParticle(ptc, 0xAA, 0xE0,	curc, (k & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE);
			}
		}
		out.replaceColor(curc, nextc);
		curc = nextc;
		nextc -= 4;
	}
	// Remove the 2pix border around all particles.
	out.replaceColor(curc, 0x80);
	out.replaceColor(0xAA, 0xFF);
	for (int k=0; k<4; k++) {
		for (Particle &ptc : particles) {
			out.expandParticle(ptc, 0xFF, 0xF0,	0x80, (k & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE);
		}
	}
	
	// Paint final picture (background=0, particle=0xFF)
	msk.fill(0);
	for (Particle &ptc : particles)
		msk.paintParticle(ptc, 0xFF);
}

void assemble_dna(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, bool validate)
{
	Assembler3d asmb(wm3d, hm3d, zm3d, mask3d);
	int n_particles = asmb.detect_slices();
	int n_origs = asmb.initial_linkage();

	int n_chips = asmb.absorb_chips();
	int n_conflicts = asmb.find_conflicts();
	int n_nonconflicts = asmb.resolve_nonconflicts();
	if (validate) {
		int n_resolved = asmb.resolve_conflicts();
		int n_crumbs = asmb.bread_crumbs();
	}

	asmb.mstack.fill(0x80);
	std::vector<Particle3D> cells;
	for (AsmOrigin &orig : asmb.origins) {
		if (!orig.valid) continue;
		cells.resize(cells.size()+1);
		Particle3D &pt = cells[cells.size()-1];
		asmb.toParticle3D(orig, pt);
		asmb.mstack.paintParticle(pt, 0);
	}
	
	asmb.all_slices.clear();
	asmb.origins.clear();

	for (Particle3D &pt : cells) {
		asmb.mstack.forceSeparation(pt, 0, 0x80, 0x50);
	}
	
	// Remove small filaments and protrusions
	asmb.mstack.expandBorders(0x80, 0, 0xC0, HOOD3D_26);
	asmb.mstack.expandBorders(0, 0xC0, 0x40, HOOD3D_6);
	asmb.mstack.replaceColor(0x40, 0);
	asmb.mstack.expandBorders(0, 0xC0, 0x40, HOOD3D_6);
	asmb.mstack.replaceColor(0x40, 0);
	asmb.mstack.replaceColor(0xC0, 0x80);
	
	asmb.mstack.rescanShrunkParticles(cells, 0, 0x50);
	asmb.mstack.expandBorders(0x80, 0, 0xFF);
	
	std::cout << "Write dna data to " << csvfile << std::endl;
	write_cell_data(cells, csvfile);
}

void segment_actin(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	int otsu1, int otsu2)
{
	Raster16 img(w, h, data);
	Raster8 ang(w, h, angle);
	Raster8 msk(w, h, mask);
	
	Raster8 out(w, h, NULL);

	for (long long i=0; i<msk.len; i++) {
		if (msk.buf[i] != 0) msk.buf[i] = 0xC0;
	}
	msk.fillBorder(0x10, 2);
	msk.expandBorders(0, 0x20, HOOD_SIZE_FATCROSS, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.expandBorders(0, 0x20, HOOD_SIZE_FATCROSS, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.expandBorders(0, 0x20, HOOD_SIZE_MOORE, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.fillBorder(0, 2);
	memcpy(out.buf, msk.buf, msk.len);
	
	std::vector<std::vector<Point>> open_paths;
	
	int bsz = 3;

	GradientWalker wk(img, ang, out, otsu1);
	wk.setStep(5.);
	wk.bend = 7;
	wk.awiggle = 28;
	wk.maxgap = 15.;
	wk.touch_angle = 95;

	for (int y0=bsz; y0<h-bsz; y0++) {
		unsigned short *dsrc = img.scanLine(y0);
		unsigned char *pout = out.scanLine(y0);
		for (int x0=bsz; x0<w-bsz; x0++) {
			if (pout[x0] > 0) continue;
			if (dsrc[x0] < otsu2 || !img.is_local_max(x0, y0)) continue;
			// out.setValue(x0, y0, 0x80);
			open_paths.resize(open_paths.size()+1);
			std::vector<Point> &path = open_paths[open_paths.size()-1];
			int l = wk.walk(path, x0, y0);
			// std::cout << "l=" << l << " touch=" << wk.touch << std::endl;
			if (l < 25 && wk.touch < 2) {
				open_paths.resize(open_paths.size()-1);
			} else {
				if (wk.in_loop)
					path.push_back(path[wk.loop_start]);
				out.paintPath(path, 0xFF, 0x80);
			}
		}
	}

	out.expandBorders(0x80, 0x40, HOOD_SIZE_FATCROSS, 0);
	out.fillBorder(0x10, 1);
	out.filterParticles(0, 0x30, 100, 0x60);
	out.replaceColor(0x30, 0);
	out.replaceColor(0x60, 0x40);

	out.filterParticles(0, 0x30, 2000, 0x41);
	out.replaceColor(0x30, 0);
	out.replaceColor(0x80, 0xFE);

	wk.setStep(4.);
	wk.bend = 7;
	wk.awiggle = 21;
	wk.maxgap = 10.;
	wk.touch_angle = 95;
	
	for (int y0=bsz; y0<h-bsz; y0++) {
		unsigned short *dsrc = img.scanLine(y0);
		unsigned char *pout = out.scanLine(y0);
		for (int x0=bsz; x0<w-bsz; x0++) {
			if (pout[x0] > 0) continue;
			if (dsrc[x0] < otsu1 || img.local_rank(x0, y0, HOOD_SIZE_FATCROSS) > 2) continue;
			open_paths.resize(open_paths.size()+1);
			std::vector<Point> &path = open_paths[open_paths.size()-1];
			int l = wk.walk(path, x0, y0);
			if (l < 25 && wk.touch < 2) {
				open_paths.resize(open_paths.size()-1);
			} else {
				if (wk.in_loop)
					path.push_back(path[wk.loop_start]);
				out.paintPath(path, 0xA0, 0x60);
			}
		}
	}
	
	out.fill(0);
	for (std::vector<Point> &path : open_paths) {
		for (size_t j=1; j<path.size(); j++)
			out.paintSegment(path[j-1], path[j], 0xFF);
	}
	for (long long i=0; i<out.len; i++) {
		if (msk.buf[i] != 0 && out.buf[i] == 0) out.buf[i] = 0xC0;
	}
	
	out.fillBorder(0x10, 2);
	out.expandBorders(0xFF, 0x80, HOOD_SIZE_FATCROSS, 0);
	out.expandBorders(0x80, 0x78, HOOD_SIZE_FATCROSS, 0);
	out.expandBorders(0x78, 0x70, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 20, 0x30);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x30, 0x80);

	// out.filterParticles(0, 0x20, 300, 0x40);
	// out.replaceColor(0x20, 0);	
	// out.expandBorders(0x78, 0x70, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x70, 0x68, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x68, 0x60, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x60, 0x58, HOOD_SIZE_FATCROSS, 0);

	out.replaceColor(0x40, 0);	
	
	out.replaceColor(0xC0, 0);
	out.fillBorder(0, 3);
	out.fillBorder(0x10, 1);
	out.fillParticle(1, 1, 0x10);
	out.filterParticles(0, 0x80, 50000, 0x20);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x10, 0x80);
	out.fillBorder(0x10, 1);

	out.filterParticles(0, 0x20, 50, 0x80);
	out.replaceColor(0x20, 0);
	
	// Try to rejoin small chips to bigger ones nearby
	out.filterParticles(0, 0x20, 200, 0xA0);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x58, 0x60);
	for (int pass=0; pass<5; pass++) {
		out.expandBorders(0xA0, 0x30, HOOD_SIZE_FATCROSS, 0x60);
		out.replaceColor(0x30, 0xA0);
	}
	out.replaceColor(0xA0, 0);
	
	std::vector<Slice> particles;
	for (int y0=1; y0<h-1; y0++) {
		unsigned char *p = out.scanLine(y0);
		for (int x0=1; x0<w-1; x0++) {
			if (p[x0] != 0) continue;
			particles.resize(particles.size()+1);
			Slice &ptc = particles[particles.size()-1];
			ptc.area = out.detectParticle(ptc, x0, y0, 0x30);
		}
	}
	for (Slice &ptc : particles)
		out.paintParticle(ptc, 0);
	int npass = 8;
	for (int lc=0x60; lc<=0x80; lc+=0x8) {
		unsigned char curc = (unsigned char)(lc);
		unsigned char nextc = curc + 0x8;
//		if (curc == 0x80) {
//			out.replaceColor(0xFF, 0x80);
//			npass = 6;
//		}
		for (int pass=0; pass<npass; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
			for (Slice &ptc : particles)
				out.expandParticle(ptc, 0, 0x30, curc, nbsz);
		}
		if (curc == 0x80) break;
		out.replaceColor(curc, nextc);
		npass = 4;
	}
	
	// Perform a couple more erosion iterations to compensate for GradientWalker bias
	out.replaceColor(0xFF, 0x80);
	for (Slice &ptc : particles)
		out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_NEUMANN);
	for (Slice &ptc : particles)
		out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_MOORE);

	out.replaceColor(out.forcedc, 0x80);
	//out.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
	out.fillBorder(0x80, 1);
	out.replaceColor(0, 0xFF);
	out.replaceColor(0x80, 0);

/*
    // Overlay DNA for debugging	
	for (long long i=0; i<msk.len; i++) {
		if (msk.buf[i] != 0 && out.buf[i] == 0) out.buf[i] = 0xA0;
	}
*/

	memcpy(msk.buf, out.buf, msk.len);
	return;
}

void segment_actin_z01(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	unsigned char *mask2,
	int otsu1, int otsu2)
{

	Raster16 img(w, h, data);
	Raster8 ang(w, h, angle);
	Raster8 msk(w, h, mask);
	Raster8 z01(w, h, mask2);
	
	Raster8 out(w, h, NULL);

	for (long long i=0; i<z01.len; i++) {
		if (z01.buf[i] != 0xFF) z01.buf[i] = 0;
	}
	z01.expandBorders(0xFF, 0x80, HOOD_SIZE_FATCROSS, 0);
	for (int mpass=0; mpass<2; mpass++) {
		z01.expandBorders(0x80, 0x60, HOOD_SIZE_FATCROSS, 0);
		z01.replaceColor(0x60, 0x80);
	}

	for (long long i=0; i<msk.len; i++) {
		if (msk.buf[i] != 0) msk.buf[i] = 0xC0;
	}
	msk.fillBorder(0x10, 2);
	msk.expandBorders(0, 0x20, HOOD_SIZE_FATCROSS, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.expandBorders(0, 0x20, HOOD_SIZE_FATCROSS, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.expandBorders(0, 0x20, HOOD_SIZE_MOORE, 0xC0);
	msk.replaceColor(0x20, 0);
	msk.fillBorder(0, 2);
	memcpy(out.buf, msk.buf, msk.len);

	std::vector<std::vector<Point>> open_paths;
	
	int bsz = 3;

	GradientWalker wk(img, ang, out, otsu1);
	wk.setStep(5.);
	wk.bend = 7;
	wk.awiggle = 28;
	wk.maxgap = 15.;
	wk.touch_angle = 95;

	for (int y0=bsz; y0<h-bsz; y0++) {
		unsigned short *dsrc = img.scanLine(y0);
		unsigned char *pout = out.scanLine(y0);
		unsigned char *zout = z01.scanLine(y0);
		for (int x0=bsz; x0<w-bsz; x0++) {
			if (pout[x0] > 0) continue;
			if (zout[x0] == 0) continue;
			if (dsrc[x0] < otsu2 || !img.is_local_max(x0, y0)) continue;
			// out.setValue(x0, y0, 0x80);
			open_paths.resize(open_paths.size()+1);
			std::vector<Point> &path = open_paths[open_paths.size()-1];
			int l = wk.walk(path, x0, y0);
			// std::cout << "l=" << l << " touch=" << wk.touch << std::endl;
			if (l < 25 && wk.touch < 2) {
				open_paths.resize(open_paths.size()-1);
			} else {
				if (wk.in_loop)
					path.push_back(path[wk.loop_start]);
				out.paintPath(path, 0xFF, 0x80);
			}
		}
	}

	out.expandBorders(0x80, 0x40, HOOD_SIZE_FATCROSS, 0);
	out.fillBorder(0x10, 1);
	out.filterParticles(0, 0x30, 100, 0x60);
	out.replaceColor(0x30, 0);
	out.replaceColor(0x60, 0x40);

	out.filterParticles(0, 0x30, 2000, 0x41);
	out.replaceColor(0x30, 0);
	out.replaceColor(0x80, 0xFE);

	for (int mpass=0; mpass<1; mpass++) {
		z01.expandBorders(0x80, 0x60, HOOD_SIZE_FATCROSS, 0);
		z01.replaceColor(0x60, 0x80);
	}

	wk.setStep(4.);
	wk.bend = 7;
	wk.awiggle = 21;
	wk.maxgap = 10.;
	wk.touch_angle = 95;
	
	for (int y0=bsz; y0<h-bsz; y0++) {
		unsigned short *dsrc = img.scanLine(y0);
		unsigned char *pout = out.scanLine(y0);
		unsigned char *zout = z01.scanLine(y0);
		for (int x0=bsz; x0<w-bsz; x0++) {
			if (pout[x0] > 0) continue;
			if (zout[x0] == 0) continue;
			if (dsrc[x0] < otsu1 || img.local_rank(x0, y0, HOOD_SIZE_FATCROSS) > 2) continue;
			open_paths.resize(open_paths.size()+1);
			std::vector<Point> &path = open_paths[open_paths.size()-1];
			int l = wk.walk(path, x0, y0);
			if (l < 25 && wk.touch < 2) {
				open_paths.resize(open_paths.size()-1);
			} else {
				if (wk.in_loop)
					path.push_back(path[wk.loop_start]);
				out.paintPath(path, 0xA0, 0x60);
			}
		}
	}
	
	out.fill(0);
	for (std::vector<Point> &path : open_paths) {
		for (size_t j=1; j<path.size(); j++)
			out.paintSegment(path[j-1], path[j], 0xFF);
	}
	for (long long i=0; i<out.len; i++) {
		if (msk.buf[i] != 0 && out.buf[i] == 0) out.buf[i] = 0xC0;
	}

	out.fillBorder(0x10, 2);
	out.expandBorders(0xFF, 0x80, HOOD_SIZE_FATCROSS, 0);
	out.expandBorders(0x80, 0x78, HOOD_SIZE_FATCROSS, 0);
	out.expandBorders(0x78, 0x70, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 20, 0x30);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x30, 0x80);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x70, 0x68, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x68, 0x60, HOOD_SIZE_FATCROSS, 0);

	out.filterParticles(0, 0x20, 300, 0x40);
	out.replaceColor(0x20, 0);
	out.expandBorders(0x60, 0x58, HOOD_SIZE_FATCROSS, 0);

	out.replaceColor(0x40, 0);	
	
	out.replaceColor(0xC0, 0);
	out.fillBorder(0, 3);
	out.fillBorder(0x10, 1);
	out.fillParticle(1, 1, 0x10);
	out.filterParticles(0, 0x80, 25000, 0x20);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x10, 0x80);
	out.fillBorder(0x10, 1);

	out.filterParticles(0, 0x20, 50, 0x80);
	out.replaceColor(0x20, 0);
	
	// Try to rejoin small chips to bigger ones nearby
	out.filterParticles(0, 0x20, 200, 0xA0);
	out.replaceColor(0x20, 0);
	out.replaceColor(0x58, 0x60);
	for (int pass=0; pass<5; pass++) {
		out.expandBorders(0xA0, 0x30, HOOD_SIZE_FATCROSS, 0x60);
		out.replaceColor(0x30, 0xA0);
	}
	out.replaceColor(0xA0, 0);
	
	std::vector<Slice> particles;
	for (int y0=1; y0<h-1; y0++) {
		unsigned char *p = out.scanLine(y0);
		for (int x0=1; x0<w-1; x0++) {
			if (p[x0] != 0) continue;
			particles.resize(particles.size()+1);
			Slice &ptc = particles[particles.size()-1];
			ptc.area = out.detectParticle(ptc, x0, y0, 0x30);
		}
	}
	for (Slice &ptc : particles)
		out.paintParticle(ptc, 0);
	int npass = 8;
	for (int lc=0x60; lc<=0x80; lc+=0x8) {
		unsigned char curc = (unsigned char)(lc);
		unsigned char nextc = curc + 0x8;
//		if (curc == 0x80) {
//			out.replaceColor(0xFF, 0x80);
//			npass = 6;
//		}
		for (int pass=0; pass<npass; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
			for (Slice &ptc : particles)
				out.expandParticle(ptc, 0, 0x30, curc, nbsz);
		}
		if (curc == 0x80) break;
		out.replaceColor(curc, nextc);
		npass = 4;
	}

	// Perform a couple more erosion iterations to compensate for GradientWalker bias
	out.replaceColor(0xFF, 0x80);
	for (Slice &ptc : particles)
		out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_NEUMANN);
	for (Slice &ptc : particles)
		out.expandParticle(ptc, 0, 0x30, 0x80, HOOD_SIZE_MOORE);

	out.replaceColor(out.forcedc, 0x80);
	out.fillBorder(0x80, 1);
	out.replaceColor(0, 0xFF);
	out.replaceColor(0x80, 0);

/*
    // Overlay DNA for debugging	
	for (long long i=0; i<msk.len; i++) {
		if (msk.buf[i] != 0 && out.buf[i] == 0) out.buf[i] = 0xA0;
	}
*/
/*
	// Overlay Z01 for debugging	
	for (long long i=0; i<msk.len; i++) {
		if (z01.buf[i] == 0xFF) out.buf[i] = 0xA0;
		else if (z01.buf[i] == 0x80 && out.buf[i] == 0) out.buf[i] = 0x40;
	}
*/

	memcpy(msk.buf, out.buf, msk.len);
}

void assemble_actin(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, const char *dna_csvfile, const char *z01_csvfile, bool validate)
{
	ActinAssembler3d asmb(wm3d, hm3d, zm3d, mask3d);
	int n_slices = asmb.detect_slices();
	
	int rc = asmb.read_z01_data(z01_csvfile);
	int n_zsl = asmb.z01_slice_linkage();
	
	rc = asmb.read_dna_data(dna_csvfile);
	int n_nucsl = asmb.nuc_slice_linkage();
	
	int n_origs = asmb.initial_linkage();
	int n_chips = asmb.absorb_chips();
	if (validate) {
		int n_links = asmb.nuc_cell_linkage();
		int n_nuc_conflicts = asmb.resolve_nuclear_conflicts();
		int n_conflicts = asmb.resolve_conflicts();
		int n_restored = asmb.fill_missing_slices();
	}
	
	asmb.mstack.fill(0x80);
	std::vector<Particle3D> cells;
	for (AsmOrigin &orig : asmb.origins) {
		if (!orig.valid) continue;
		cells.resize(cells.size()+1);
		Particle3D &pt = cells[cells.size()-1];
		asmb.toParticle3D(orig, pt);
		// asmb.mstack.paintParticle(pt, 0);
	}
	
	asmb.all_slices.clear();
	asmb.origins.clear();
	
	for (int z=0; z<asmb.mstack.d; z++) {
		Raster8 msk = asmb.mstack.getPlane(z);
		for (Particle3D& pt : cells) {
			if (pt.fills[z].empty()) continue;
			Boundary bnd = pt.bnd.boundary2d();
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.paintParticleFillInto(pt.fills[z], 0xC0, 0x80);
			msk.rescanParticleFill(bnd, pt.fills[z], 0xC0);
			msk.bordersAround(bnd, 0xC0, 0x80, 0x40, HOOD_SIZE_MOORE);
			msk.replaceColor(bnd, 0xC0, 0);
		}
		msk.replaceColor(0x40, 0x80);
	}
	
	asmb.mstack.expandBorders(0x80, 0, 0xFF);
	
	std::cout << "Write Actin data to " << csvfile << std::endl;
	write_cell_data(cells, csvfile);
}

void detect_top_bottom(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
	const char *actin_csvfile, const char *dna_csvfile)
{
	TopBottomDetector tbd(wm3d, hm3d, zm3d, mask3d);
	int dna_id = tbd.read_dna_data(dna_csvfile);
	int cell_id = tbd.read_actin_data(actin_csvfile);
	
// std::cout << tbd.nuclei.size() << " nuclei read from " << dna_csvfile << std::endl;
// std::cout << tbd.cells.size() << " cells read from " << actin_csvfile << std::endl;
	
	tbd.map_cells_to_nuclei();
	tbd.actin_cell_limits();
	tbd.detect_top_bottom();

	tbd.mstack.fill(0x80);
	for (Cell &cell : tbd.cells) {
		tbd.mstack.paintParticle(cell, 0);
	}
	for (Cell &cell : tbd.cells) {
		tbd.mstack.forceSeparation(cell, 0, 0x80, 0x50);
	}

	// Remove small filaments and protrusions
	tbd.mstack.expandBorders(0x80, 0, 0xA0, HOOD3D_26);
	tbd.mstack.expandBorders(0xA0, 0, 0xC0, HOOD3D_26);
	tbd.mstack.replaceColor(0xA0, 0xC0);
	for (int pass=0; pass<4; pass++) {
		tbd.mstack.expandBorders(0, 0xC0, 0x40, HOOD3D_6);
		tbd.mstack.replaceColor(0x40, 0);
	}
	tbd.mstack.replaceColor(0xC0, 0x80);
	
	tbd.mstack.replaceColor(0xA0, 0x80);
	tbd.mstack.rescanShrunkCells(tbd.cells, 0, 0x50);
	tbd.mstack.expandBorders(0x80, 0, 0xFF);
	
/* DEBUG -- overlay DNA data
	for (int z=0; z<tbd.mstack.d; z++) {
		Raster8 msk = tbd.mstack.getPlane(z);
		for (Nucleus &nuc : tbd.nuclei) {
			if (z < nuc.bnd.zmin || z > nuc.bnd.zmax) continue;
			msk.paintParticleFillInto(nuc.fills[z], 0x40, 0);
			msk.paintParticleFillInto(nuc.fills[z], 0xA0, 0x80);
		}
	}
*/
	std::cout << "Write Actin data to " << actin_csvfile << std::endl;
	write_cell_data(tbd.cells, actin_csvfile);
}

// Used by extract_z01()
static void smooth_z_mask(Raster8 &xmsk)
{
	Raster8 tmsk(xmsk.w, xmsk.h, NULL);
	memcpy(tmsk.buf, xmsk.buf, xmsk.len);
	for (int y=2; y<xmsk.h-2; y++) {
		unsigned char *p = tmsk.scanLine(y);
		for (int x=2; x<xmsk.w-2; x++) {
			int csum = 0;
			for (int j=0; j<HOOD_SIZE_FATCROSS; j++) {
				csum += int(xmsk.value(x+hood_pts[j].dx, y+hood_pts[j].dy));
			}
			p[x] = (unsigned char)(csum / HOOD_SIZE_FATCROSS);
		}
	}
	memcpy(xmsk.buf, tmsk.buf, xmsk.len);
}

// Input: 3D array of 16-bit (n_frames, height, width) containing Z01 channel raw data;
// Output:	Frame 0 replaced with aggregated 2D Z01 data
//			Frame 1 - Lower Z (ranging 0 through n_frames-1)
//			Frame 2 - Upper Z (ranging 0 through n_frames-1)
void extract_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d)
{
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);

	// Detect frame range loz..hiz suitable for computing Otsu
	// Try to exclude first 3 frames often containing heavy static
	double *fr_dev = new double[dstack.d];
	double hidev = 0.;
	int loz = 0, hiz = 0;
	for (int z=dstack.d-1; z>=0; z--) {
		Raster16 dat = dstack.getPlane(z);
		double dev;
		dat.mean_std(&dev);
		fr_dev[z] = dev;
		if (hidev < dev && z > 2) {
			hidev = dev;
			hiz = z;
		}
	}
	double middev = hidev * 0.75;
	double upperdev = hidev * 1.02;
	for (loz=hiz-1; loz>=0; loz--) {
		if (fr_dev[loz] < middev || fr_dev[loz] >= upperdev) break;
	}
	++loz;
	for (; hiz<dstack.d; hiz++) {
		if (fr_dev[hiz] < middev) break;
	}
	--hiz;
	delete [] fr_dev;
	
	// std::cout << "loz=" << loz << " hiz=" << hiz << std::endl;
	if (loz < 0) loz = 0;
	if (hiz >= dstack.d) hiz = dstack.d - 1;
	if (loz > hiz) loz=hiz;
	
	Histogram hist(0x1000);
	for (int z=loz; z<=hiz; z++) {
		Raster16 dat = dstack.getPlane(z);
		for (int y=0; y<dat.h; y++) {
			hist.add_row16(dat.scanLine(y), dat.w);
		}
	}
	
	unsigned short otsu = hist.otsu16();
	// Half Otsu seems to work best for any week W1-W6
	otsu /= 2;

	Raster3D mstack(wd3d, hd3d, zd3d, NULL);
	for (long long i=0; i<mstack.len; i++) {
		mstack.buf[i] = (dstack.buf[i] > otsu) ? 0xFF : 0;
	}
	
	Raster8 tmsk(wd3d, hd3d, NULL);
	Raster8 xmsk(wd3d, hd3d, NULL);
	xmsk.fill(0xFF);
	
	// Start from the top, follow the pattern.
	// Exclude anything that conflicts with pattern found 8 (or more) frames above the current frame.
	for (int z=mstack.d-1; z>=0; z--) {
		Raster8 msk = mstack.getPlane(z);
		msk.expandBorders(0xFF, 0x80, HOOD_SIZE_MOORE, 0);
		if (z < mstack.d-1) {
			Raster8 top = mstack.getPlane(z+1);
			for (long long i=0; i<msk.len; i++)
				if (top.buf[i] == 0xFF && msk.buf[i] != 0xFF) msk.buf[i] = 0x80;
		}
		if (z > 0) {
			Raster8 bot = mstack.getPlane(z-1);
			for (long long i=0; i<msk.len; i++)
				if (bot.buf[i] == 0xFF && msk.buf[i] != 0xFF) msk.buf[i] = 0x80;
		}
		msk.expandBorders(0x80, 0x40, HOOD_SIZE_MOORE, 0);
		msk.replaceColor(0x40, 0xFF);
		msk.replaceColor(0x80, 0xFF);
		msk.fillBorder(0x10, 1);
		msk.filterParticles(0xFF, 0x80, 1000, 0);
		
		if (z+6 < mstack.d-1)
		{
			Raster8 top = mstack.getPlane(z+1);
			Raster8 top2 = mstack.getPlane(z+2);
			for (int y0=1; y0<msk.h-1; y0++) {
				unsigned char *p = msk.scanLine(y0);
				for (int x0=1; x0<msk.w-1; x0++) {
					if (p[x0] != 0x80) continue;
					Slice ptc;
					ptc.area = msk.detectParticle(ptc, x0, y0, 0xA0);
					
					int ncgood = 0;
					for (HSeg &hs : ptc.fill) {
						unsigned char *p1 = top.scanLine(hs.y);
						unsigned char *p2 = top2.scanLine(hs.y);
						for (int x=hs.xl; x<=hs.xr; x++) {
							if (p1[x] >= 0x80 || p2[x] >= 0x80)
								++ncgood;
						}
					}
					
					int ncbad = ptc.area - xmsk.countColor(ptc, 0xFF);
					
					if (ncbad > int(ptc.area * 0.4) ||
							(ncbad > ncgood && ncbad > int(ptc.area * 0.05))) {
						msk.paintParticle(ptc, 0x40);
					}
				}
			}
			msk.replaceColor(0xA0, 0x80);
			msk.replaceColor(0x40, 0);
		}
		
		
		for (int pass=0; pass<3; pass++) {
			msk.expandBorders(0x80, 0x40, HOOD_SIZE_RAD3, 0);
			msk.replaceColor(0x40, 0x80);
		}
		msk.filterParticles(0, 0x20, 10000, 0x40);
		msk.replaceColor(0x40, 0x80);
		msk.replaceColor(0x20, 0);
		
		if (z+8 < mstack.d-1) {
			Raster8 top = mstack.getPlane(z+8);
			for (long long i=0; i<xmsk.len; i++)
				if (top.buf[i] >= 0x80 && xmsk.buf[i] == 0xFF)
					xmsk.buf[i] = (unsigned char)(z);
		}
		
		Raster16 dat = dstack.getPlane(z);
		for (long long i=0; i<msk.len; i++) {
			if (msk.buf[i] == 0x80 && dat.buf[i] > otsu)
				msk.buf[i] = 0xFF;
		}
	}
	
	// Find upper Z for every pixel in the detected pattern
	xmsk.fill(0xFF);
	for (int z=mstack.d-1; z>=0; z--) {
		int zz = z + 1;
		if (zz >= mstack.d) zz = mstack.d - 1;
		else if (zz < 3) zz = 3;
		Raster8 msk = mstack.getPlane(z);
		for (long long i=0; i<xmsk.len; i++)
			if (msk.buf[i] >= 0x80 && xmsk.buf[i] == 0xFF)
				xmsk.buf[i] = (unsigned char)(zz);
	}
	
	// Try to fill the gaps (averaging Z-values at the perimeter)
	xmsk.fillBorder(0xFE, 1);
	for (int y0=1; y0<xmsk.h-1; y0++) {
		unsigned char *p = xmsk.scanLine(y0);
		for (int x0=1; x0<xmsk.w-1; x0++) {
			if (p[x0] != 0xFF) continue;
			Slice ptc;
			ptc.area = xmsk.detectParticle(ptc, x0, y0, 0xF0);
			
			Boundary bnd = ptc.bnd;
			bnd.expand(3);
			xmsk.clip(bnd);
			bnd.expand(-2);
			
			int npts = 1;
			int csum = dstack.d - 1;

			for (int y=bnd.ymin; y<=bnd.ymax; y++) {
				for (int x=bnd.xmin; x<=bnd.xmax; x++) {
					unsigned char c = xmsk.value(x,y);
					if (c >= 0xF0) continue;
					for (int j=1; j<HOOD_SIZE_FATCROSS; j++) {
						if (xmsk.value(x+hood_pts[j].dx, y+hood_pts[j].dy) >= 0xF0) {
							++npts;
							csum += int(c);
							break;
						}
					}
				}
			}
			if (npts > 0) {
				csum = (csum / npts) + 3;
				if (csum >= dstack.d) csum = dstack.d - 1;
			}

			xmsk.paintParticle(ptc, (unsigned char)csum);
		}
	}
	
	// Fill border values by duplicating the values next to border
	for (int y0=1; y0<xmsk.h-1; y0++) {
		unsigned char *p = xmsk.scanLine(y0);
		p[0] = p[1];
		int x = xmsk.w-1;
		p[x] = p[x-1];
	}
	memcpy(xmsk.scanLine(0), xmsk.scanLine(1), xmsk.w * sizeof(unsigned char));
	memcpy(xmsk.scanLine(xmsk.h-1), xmsk.scanLine(xmsk.h-2), xmsk.w * sizeof(unsigned char));
	
	// Smooth it a little
	smooth_z_mask(xmsk);
	
	// Now find lower Z 
	tmsk.fill(0xFF);
	for (int z=0; z<mstack.d-1; z++) {
		int zz = z - 1;
		if (zz < 0) zz = 0;
		else if (zz+3 >= mstack.d-1) zz = mstack.d-4;
		Raster8 msk = mstack.getPlane(z);
		for (long long i=0; i<tmsk.len; i++)
			if (msk.buf[i] >= 0x80 && tmsk.buf[i] == 0xFF)
				tmsk.buf[i] = (unsigned char)(zz);
	}
	// Fill gaps with values UpperZ - 8
	for (long long i=0; i<xmsk.len; i++) {
		if (tmsk.buf[i] == 0xFF) {
			int c = int(xmsk.buf[i]) - 8;
			if (c < 0) c = 0;
			tmsk.buf[i] = (unsigned char)(c);
		}
	}
	
	// Smooth lower Z
	smooth_z_mask(tmsk);

	// Aggregate pixel values by taking an average of 3 brightest pixels in the range LowerZ..UpperZ
	// Store the result in frame 0
	for (int y=0; y<dstack.h; y++) {
		unsigned short* tgt = dstack.scanLine(y, 0);
		unsigned char *plo = tmsk.scanLine(y);
		unsigned char *phi = xmsk.scanLine(y);
		for (int x=0; x<dstack.w; x++) {
			int z0 = int(plo[x]);
			int z1 = int(phi[x]);
			if (z1-z0 < 4) {
				z1 = z0 + 4;
				if (z1 >= dstack.d) z1 = dstack.d - 1;
			}
			if (z1-z0 < 2) {
				z0 = z1 - 2;
				if (z0 < 0) z0 = 0;
			}
			unsigned short v1=0, v2=0, v3=0;
			for (int z=z0; z<=z1; z++) {
				unsigned short v = dstack.value(x, y, z);
				if (v > v1) {
					v3 = v2;
					v2 = v1;
					v1 = v;
				} else if (v > v2) {
					v3 = v2;
					v2 = v;
				} else if (v > v3) {
					v3 = v;
				}
			}
			tgt[x] = (unsigned short)((int(v1) + int(v2) + int(v3)) / 3);
		}
	}

	// Store LowerZ values in frame 1, and UpperZ, in frame 2
	unsigned short *fr1 = dstack.scanLine(0, 1);
	unsigned short *fr2 = dstack.scanLine(0, 2);
	for (long long i=0; i<xmsk.len; i++) {
		fr1[i] = (unsigned short)(tmsk.buf[i]); // << 3;
		fr2[i] = (unsigned short)(xmsk.buf[i]); // << 3;
	}
}

void combo_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d)
{
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	
	Histogram hist(0x4000);
	for (int z=0; z<dstack.d; z++) {
		for (int y=0; y<dstack.h; y++)
			hist.add_row16(dstack.scanLine(y, z), dstack.w);
	}
	unsigned short otsu = hist.otsu16();
	// std::cout << "OTSU: " << otsu << std::endl;
	
	std::vector<int> highs(dstack.d, 0);
	std::vector<int> matches(dstack.d, 0);
	int maxhigh = 1;
	for (int z=1; z<dstack.d; z++) {
		int ntot = 1;
		int nmatch = 0;
		// int nmism = 0;
		for (int y=0; y<dstack.h; y++) {
			unsigned short *p0 = dstack.scanLine(y, z-1);
			unsigned short *p1 = dstack.scanLine(y, z);
			for (int x=0; x<dstack.w; x++) {
				bool m0 = p0[x] > otsu;
				bool m1 = p1[x] > otsu;
				if (m0 || m1) {
					++ntot;
					if (m0 && m1) ++nmatch;
				}
			}
		}
		if (maxhigh < ntot) maxhigh = ntot;
		highs[z] = ntot + 1;
		matches[z] = nmatch;
	}
	
	std::vector<double> ratings(dstack.d, 0.);
	double maxrat = 0.;
	int bestz = 0;
	for (int z=1; z<dstack.d; z++) {
		int high = highs[z];
		int match = matches[z];
		double rating = (double(match) * high) / (double(high - match) * maxhigh);
		ratings[z] = rating;
		if (maxrat < rating) {
			maxrat = rating;
			bestz = z;
		}
	}
	
	ratings[0] = ratings[1];
	
	std::vector<int> locmin(dstack.d, 0);
	for (int z=0; z<dstack.d; z++) {
		int z0 = z - 3; if (z0 < 0) z0 = 0;
		int z1 = z + 3; if (z1 >= dstack.d) z1 = dstack.d - 1;
		double xrat = maxrat + 1.;
		int xz = -1;
		for (int zz=z0; zz<=z1; zz++) {
			if (ratings[zz] < xrat) {
				xrat = ratings[zz];
				xz = zz;
			}
		}
		if (xz == z) locmin[z] = 1;
	}
	
	double rthresh = maxrat * 0.07;
	int zmin, zmax;
	for (zmin=bestz; zmin>0; zmin--) {
		if (ratings[zmin] < rthresh || locmin[zmin] != 0) {
			--zmin;
			break;
		}
	}
	for (zmax=bestz; zmax<dstack.d; zmax++) {
		if (locmin[zmax] != 0) break;
		if (ratings[zmax] < rthresh)
			break;
	}

	std::vector<double> vmax1(dstack.w);
	std::vector<double> vmax2(dstack.w);
	std::vector<double> vmax3(dstack.w);
	
	for (int y=0; y<dstack.h; y++) {
		for (size_t i=0; i<vmax1.size(); i++) {
			vmax1[i] = vmax2[i] = vmax3[i] = 0.;
		}
		for (int z=zmin; z<=zmax; z++) {
			unsigned short *s = dstack.scanLine(y, z);
			for (int x=0; x<dstack.w; x++) {
				double v = double(s[x]);
				if (vmax1[x] < v) {
					vmax3[x] = vmax2[x];
					vmax2[x] = vmax1[x];
					vmax1[x] = v;
				}
				else if (vmax2[x] < v) {
					vmax3[x] = vmax2[x];
					vmax2[x] = v;
				}
				else if (vmax3[x] < v) {
					vmax3[x] = v;
				}
			}
		}
		unsigned short *t = dstack.scanLine(y, 0);
		for (int x=0; x<dstack.w; x++) {
			t[x] = (unsigned short)((vmax1[x] + vmax2[x] + vmax3[x]) * 0.33333);
		}
	}
	
}

void normalize_frame(unsigned short *sdata, int shd, int swd, int kernel_size)
{
	int hf_size = kernel_size / 2;
	Raster16 dat(swd, shd, sdata);
	Raster16 davg(swd, shd, NULL);

	SummedArea<double> sumar(swd, shd);
	for (int y=0; y<dat.h; y++)
		sumar.add_line(y, dat.scanLine(y));
	
	double amax = 0.;
	double amin = 1e20;
	for (int y=0; y<dat.h; y++) {
		unsigned short *t = davg.scanLine(y);
		for (int x=0; x<dat.w; x++) {
			double dv = sumar.avg_around(x, y, hf_size);
			if (amax < dv) amax = dv;
			if (amin > dv) amin = dv;
			t[x] = (unsigned short)(dv);
		}
	}
	
	double adelta = amax / 20.;
	if (adelta < 1.) adelta = 1.;
	amax = (amax + adelta) / 10.;
	for (int y=0; y<dat.h; y++) {
		unsigned short *p = davg.scanLine(y);
		unsigned short *t = dat.scanLine(y);
		for (int x=0; x<dat.w; x++) {
			double sc = amax / (p[x] + adelta);
			if (sc > 1.) sc = 1.;
			// p[x] = (unsigned short)(sc * 250.);
			t[x] = (unsigned short)(t[x] * sc);
		}
	}
}

int adjusted_threshold(unsigned short *data, int hd, int wd, int thresh, double pct)
{
	if (thresh < 1) thresh = 1;
	int athresh = thresh;
	
	size_t *hist = new size_t[thresh+1];
	memset(hist, 0, (thresh+1)*sizeof(size_t));

	size_t len = size_t(hd) * size_t(wd);
	for (size_t i=0; i<len; i++) {
		int v = int(data[i]);
		if (v > thresh) v = thresh;
		++hist[v];
	}
	
	size_t acc = 0;
	for (int lvl=0; lvl<thresh; lvl++) {
		acc += hist[lvl];
		double p = double(acc) / len;
		if (p >= pct) {
			athresh = lvl;
			break;
		}
	}
	
	delete [] hist;
	
	return athresh;
}

// mask -- segmentation data frame (0=background, 0xFF=foreground);
// csvfile -- file name to write 2d particle data (ID,y,xL,xR) to;
// On return 'mask' replaced with: 0x80=background, 0=foreground, 0xFF=foreground/border
void export_2d_segmentation(unsigned char *mask, int hm, int wm, const char *csvfile)
{
	Raster8 msk(wm, hm, mask);
	msk.fillBorder(0x10, 1);
	std::vector<Particle> particles;
	for (int y0=1; y0<hm-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=1; x0<wm-1; x0++) {
			if (p[x0] != 0xFF) continue;
			particles.resize(particles.size() + 1);
			Particle& ptc = particles[particles.size() - 1];
			msk.detectParticle(ptc, x0, y0, 0x80);
		}
	}
	
	msk.fill(0x80);
	for (Particle & ptc : particles)
		msk.paintParticle(ptc, 0);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
	
	CsvWriter wr(csvfile);
	if (!wr.is_open()) return;
	
	wr.append("ID");
	wr.append("y");
	wr.append("xL");
	wr.append("xR");
	wr.next();
	int id = 1;
	for (Particle & ptc : particles) {
		for (HSeg& hs : ptc.fill) {
			wr.append(id);
			wr.append(hs.y);
			wr.append(hs.xl);
			wr.append(hs.xr);
			wr.next();
		}
		++id;
	}

	wr.close();
}

// csvfile -- file name to read 2d particle data (ID,y,xL,xR) from;
// On return 'mask' replaced with: 0x80=background, 0=foreground, 0xFF=foreground/border
// Return value: # of particles read from csv (0 if csv does not exist or has wrong format)
int import_2d_segmentation(unsigned char *mask, int hm, int wm, const char *csvfile)
{
	std::vector<Particle> particles;
	if (read_particle_data(csvfile, particles, wm, hm) < 0)
		return 0;
	if (particles.empty())
		return 0;
	Raster8 msk(wm, hm, mask);
	msk.fill(0x80);
	for (Particle & ptc : particles)
		msk.paintParticle(ptc, 0);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
	return int(particles.size());
}

// csvfile -- file name to read 3d particle data (ID,Frame,y,xL,xR) from;
// border -- a boolean value, if true:
// 		on return 'mask' replaced with: 0=background, 0x80=foreground, 0xFF=2D border
// if false:
// 		on return 'mask' replaced with: 0=background, 0xFF=foreground/border
// Return value: # of cells read from csv (0 if csv does not exist or has wrong format)
int import_3d_segmentation(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, bool border)
{
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	mstack.fill(0);

	std::vector<Particle3D> cells;
	if (read_cell_data(csvfile, cells, wm3d, hm3d, zm3d) < 0)
		return 0;
	if (cells.empty())
		return 0;
	
	if (border) {
		for (Particle3D& cell : cells)
			mstack.paintParticle(cell, 0x80);
		for (int z=0; z<mstack.d; z++) {
			Raster8 msk = mstack.getPlane(z);
			msk.expandBorders(0, 0xFF, HOOD_SIZE_MOORE, 0x80);
		}
	} else {
		for (Particle3D& cell : cells)
			mstack.paintParticle(cell, 0xFF);
	}

	return int(cells.size());
}

//--- Output gray mask color coding ---
// 0x30 - RPE cell borders(outer)
// 0x40 - RPE unmatched
// 0x50 - REShAPE cell borders (outer)
// 0x60 - REShAPE unmatched
// 0x80 - RPE/REShAPE bad match
// 0xD0 - RPE/REShAPE half match
// 0xFF - RPE/REShAPE good match
void compare_with_reshape(int w, int h,
	unsigned char *mask,
	unsigned char *rs_mask,
	const char *rs_csvfile,
	const char *out_csvfile)
{
	SegmentationComparator sc(w, h, mask, rs_mask);
	sc.load_cells();
	// std::cout << "Loaded " << sc.cells.size() << " RPE_Map cells." << std::endl;
	sc.load_rs_cells(rs_csvfile);
	// std::cout << "Loaded " << sc.rs_cells.size() << " REShAPE cells from " << rs_csvfile << std::endl;
	sc.match_cells();

	// Paint fills
	sc.msk.fill(0);
	for (int idx=0; size_t(idx)<sc.cells.size(); idx++) {
		if (sc.has_rs_match(idx)) continue;
		sc.msk.paintParticle(sc.cells[idx], 0x40);
	}
	
	CsvWriter wr(out_csvfile);
	wr.append("ID");
	wr.append("XStart");
	wr.append("YStart");
	wr.append("Area");
	// wr.append("MatchArea");
	wr.append("Match");
	if (wr.next()) {
		std::cout << "Write: " << out_csvfile << std::endl;
	}
	
	std::vector<Slice *> matches;
	for (int rs_idx=0; size_t(rs_idx)<sc.rs_cells.size(); rs_idx++) {
		Slice &rcell = sc.rs_cells[rs_idx];
		int minovl = int(sc.movl * rcell.area);
		int ovl = sc.find_matches(matches, rs_idx);
		
		wr.append(rs_idx+1);
		wr.append(rcell.x0);
		wr.append(rcell.y0);
		wr.append(rcell.area);
		// wr.append(ovl);
		
		sc.msk.paintParticleInto(rcell, 0x20, 0);
		sc.msk.paintParticleInto(rcell, 0x80, 0x40);
		for (Slice *pcell : matches) {
			sc.msk.paintParticleInto(*pcell, 0x40, 0);
			sc.msk.paintParticleInto(*pcell, 0x80, 0x60);
		}
		unsigned char c = 0x80;
		if (ovl >= minovl) {
			if (matches.size() == 1) {
				c = 0xFF;
				wr.append("Good");
			} else {
				c = 0xD0;
				wr.append("Half");
			}
		} else {
			wr.append((matches.size() == 0) ? "No" : "Bad");
		}
		wr.next();
		for (Slice *pcell : matches)
			sc.msk.paintParticleInto(*pcell, c, 0x20);
		sc.msk.paintParticleInto(rcell, 0x60, 0x20);
	}
	
	// Paint borders
	sc.rmsk.fill(0);
	for (Slice &cell : sc.rs_cells)
		sc.rmsk.paintParticle(cell, 0x20);
	sc.rmsk.expandBorders(0x20, 0x50, HOOD_SIZE_MOORE, 0);
	sc.rmsk.replaceColor(0x20, 0);
	for (Slice &cell : sc.cells)
		sc.rmsk.paintParticleInto(cell, 0x20, 0);
	sc.rmsk.expandBorders(0x20, 0x30, HOOD_SIZE_MOORE, 0);
	sc.rmsk.replaceColor(0x20, 0);
	for (long long i=0; i<sc.msk.len; i++) {
		if (sc.rmsk.buf[i] != 0) sc.msk.buf[i] = sc.rmsk.buf[i];
	}
	
	// Export RPE Map segmentation	
	sc.rmsk.fill(0);
	for (Slice &cell : sc.cells)
		sc.rmsk.paintParticle(cell, 0x80);
	sc.rmsk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
}

void colorize_reshape_comparison(int w, int h,
	unsigned char *mask3d,
	unsigned short *data,
	unsigned char *mask2)
{
	double gamma = 0.75;
	long long len = (long long)(h) * (long long)(w);
	
	int hi_val;

	int hlen = 0x2000;
	long long *hist = new long long[hlen];
	memset(hist, 0, size_t(hlen)*sizeof(long long));
	for (long long i=0; i<len; i++) {
		++ hist[size_t(data[i])>>2];
	}
	long long hi_cnt = (long long)(0.99 * len);
	long long cnt = len;
	for (hi_val=hlen-1; hi_val>=0; hi_val--) {
		cnt -= hist[hi_val];
		if (cnt < hi_cnt) break;
	}
	hi_val <<= 2;
	if (hi_val < 1) hi_val = 1;
	
	unsigned char *rgb = mask3d;
	unsigned short *dp = data;
	unsigned char *mp = mask2;
	for (long long i=0; i<len; i++) {
		unsigned short dc = *dp++;
		unsigned char cc = *mp;
		
		double scv = double(dc) / hi_val;
		if (scv > 1.) scv = 1.;
		int c = (int)(pow(scv, gamma) * 0xE0);
		int r=c, g=c, b=c;
		
		*mp++ = (unsigned char)(c);
		
		if (cc == 0x30 || cc == 0x50) {
			r >>= 1; g >>= 1; b >>= 1;
		} else if (cc != 0) {
			r = r*3/4; g = g*3/4; b = b*3/4;
		}

		switch (cc) {
			case 0x30:			// RPE cell borders(outer)
				b += 0xA0;
				break;
			case 0x50:			// REShAPE cell borders (outer)
				r += 0xA0;
				break;
			case 0x40:			// RPE unmatched
				b += 0x40;
				break;
			case 0x60:			// REShAPE unmatched
				r += 0x40;
				break;
			case 0x80:			// REShAPE/RPE bad match
				r += 0x38;
				b += 0x30;
				break;
			case 0xD0:			// REShAPE/RPE half match
				g += 0x38;
				r += 0x40;
				break;
			case 0xFF:			// REShAPE/RPE good match
				g += 0x40;
				break;
		}

		if (r > 0xFF) r = 0xFF;
		if (g > 0xFF) g = 0xFF;
		if (b > 0xFF) b = 0xFF;
		*rgb++ = (unsigned char)(r);
		*rgb++ = (unsigned char)(g);
		*rgb++ = (unsigned char)(b);
	} 
}

void read_reshape(unsigned char *mask, int hm, int wm, const char *csvfile)
{
	Raster8 msk(wm, hm, mask);
	Boundary mbnd = msk.getBoundary();
	std::vector<Slice> cells;
	
	CsvReader rdr(csvfile);
	if (rdr.next()) {
		CsvRow headers = rdr.GetRow();
		rdr.setHeaders(headers);

		int xs_idx = rdr.header_index("XStart");
		int ys_idx = rdr.header_index("YStart");
		while (rdr.next()) {
			int x0 = rdr.GetInt(xs_idx, -1);
			int y0 = rdr.GetInt(ys_idx, -1);
			if (!mbnd.IsInside(x0, y0)) continue;
			cells.resize(cells.size()+1);
			Slice &cell = cells[cells.size()-1];
			cell.area = msk.detectParticle(cell, x0, y0, 0x80);
		}
	}
	msk.fill(0);
	if (cells.empty()) {
		std::cout << "Failed to read anything from " << csvfile << std::endl;
		return;
	} else {
		std::cout << "Loaded " << cells.size() << " cells from " << csvfile << std::endl;
	}
	for (Slice &cell : cells)
		msk.paintParticle(cell, 0x80);
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);
}



static bool match_particle(Slice &mptc, std::vector<Slice> &zparticles)
{
	int m80 = int(mptc.area * 0.8);
	for (Slice &zptc : zparticles) {
		int ovl = mptc.overlay_area(zptc);
		if (ovl < m80) continue;
		if (ovl < int(zptc.area * 0.8)) continue;
		return true;
	}
	return false;
}

void match_z01(unsigned char *mask, int hm, int wm, unsigned char *mask2, int hm2, int wm2)
{
	Raster8 msk(wm, hm, mask);
	Raster8 z01(wm, hm, mask2);
	
	msk.fillBorder(0x10, 1);
	z01.fillBorder(0x10, 1);
	
	std::vector<Slice> mparticles;
	std::vector<Slice> zparticles;
	
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *mp = msk.scanLine(y0);
		unsigned char *zp = z01.scanLine(y0);
		for (int x0=1; x0<msk.w-1; x0++) {
			if (mp[x0] == 0xFF) {
				mparticles.resize(mparticles.size() + 1);
				Slice &ptc = mparticles[mparticles.size() - 1];
				ptc.area = msk.detectParticle(ptc, x0, y0, 0xC0);
			}
			if (zp[x0] == 0x80) {
				zparticles.resize(zparticles.size() + 1);
				Slice &ptc = zparticles[zparticles.size() - 1];
				ptc.area = z01.detectParticle(ptc, x0, y0, 0xC0);
			}
		}
	}
	
	msk.fill(0x80);
	for (Slice &mptc : mparticles) {
		if (match_particle(mptc, zparticles))
			msk.paintParticle(mptc, 0);
		else
			msk.paintParticle(mptc, 0x40);
	}
	msk.expandBorders(0x80, 0xFF, HOOD_SIZE_MOORE, 0);

}

//--- export_mask_id_3d()
static bool compare_cell_sort(SortScoreML& ss1, SortScoreML& ss2)
{
	return ss1.sc > ss2.sc;
}

static void dilate_cells(std::vector<Particle3D>& cells, int w, int h, int d, int nd)
{
	if (cells.empty()) return;
	Raster8 msk(w, h, NULL);
	
	std::vector<SortScoreML> sort_cells(cells.size());
	for (int idx=0; size_t(idx)<cells.size(); idx++) {
		sort_cells[idx].ptid = idx;
		sort_cells[idx].sc = cells[idx].iou_score(3);
	}
	std::sort(sort_cells.begin(), sort_cells.end(), compare_cell_sort);
	
	for (int pass=0; pass<nd; pass++) {
		int nbsz = (pass&1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
		for (int z=0; z<d; z++) {
			msk.fill(0);
			msk.fillBorder(0x10, 1);
			for (int ii=0; size_t(ii)<sort_cells.size(); ii++) {
				Particle3D& cell = cells[sort_cells[ii].ptid];
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				msk.paintParticleFillInto(fill, 0x80, 0);
				Boundary bnd = fill_boundary(fill);
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.expandBordersInto(bnd, 0x80, 0, 0x50,nbsz, false);
				msk.rescanParticleFill(bnd, fill, 0x80);
				msk.paintParticleFill(fill, 0xC0);
			}
		}
	}
}

int export_mask_id_3d(unsigned short *data3d, int zd3d, int hd3d, int wd3d, const char *csvfile, int num_dilations)
{
	std::vector<Particle3D> cells;
	if (read_cell_data(csvfile, cells, wd3d, hd3d, zd3d) < 0)
		return 0;
	if (cells.empty())
		return 0;
	
	// std::cout << "Read " << cells.size() << " cells from " << csvfile << std::endl;
	
	if (num_dilations > 0) {
		dilate_cells(cells, wd3d, hd3d, zd3d, num_dilations);
	}
	
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	memset(dstack.buf, 0, dstack.len * sizeof (unsigned short));
	
	int id = 0;
	for (Particle3D& cell : cells) {
		int zmin = -1;
		int zmax = -1;
		for (int z=0; z<zd3d; z++) {
			std::vector<HSeg> &fill = cell.fills[z];
			if (fill.empty()) continue;
			zmax = z;
			if (zmin < 0) zmin = z;
		}
		if (zmax - zmin < 3) continue;
		++id;
		
		for (int z=zmin; z<=zmax; z++) {
			Raster16 dat = dstack.getPlane(z);
			dat.paintParticleFill(cell.fills[z], (unsigned short)id);
		}
	}
	
	// std::cout << "Validated: " << id << " cells." << std::endl;
	
	return id;
}

int import_mask_id_3d(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile)
{
	int rc = 0;
	
	Raster16_3D dstack(wd3d, hd3d, zd3d, data3d);
	Raster3D mstack(wm3d, hm3d, zm3d, mask3d);
	
	std::vector<Particle3D> cells;
	std::vector<std::vector<AsmSlice>> all_slices(dstack.d);
	
	mstack.fill(0);
	for (int z=0; z<dstack.d; z++) {
		std::vector<AsmSlice>& slices = all_slices[z];
		Raster16 dat = dstack.getPlane(z);
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10);
		
		int max_id = 0;
		for (long long i=0; i<dat.len; i++) {
			int id = int(dat.buf[i]);
			if (max_id < id) max_id = id;
		}
		if (max_id == 0) continue;
		slices.resize(max_id);
		
		for (int idx=0; size_t(idx)<slices.size(); idx++) {
			AsmSlice& pt = slices[idx];
			pt.z = z;
			Boundary bnd(-1, -1, 0, 0);
			unsigned short id = (unsigned short)(idx + 1);
			for (int y=0; y<dat.h; y++) {
				unsigned short *src = dat.scanLine(y);
				unsigned char *tgt = msk.scanLine(y);
				for (int x=0; x<dat.w; x++) {
					if (src[x] != id) continue;
					if (bnd.xmin < 0 || bnd.xmin > x) bnd.xmin = x;
					if (bnd.xmax < x) bnd.xmax = x;
					if (bnd.ymin < 0 || bnd.ymin > y) bnd.ymin = y;
					if (bnd.ymax < y) bnd.ymax = y;
					if (tgt[x] == 0) tgt[x] = 0x80;
				}
			}
			pt.bnd = bnd;
			pt.area = msk.rescanParticle(pt, 0x80);
			msk.paintParticle(pt, 0xFF);
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.expandBordersInto(bnd, 0xFF, 0, 0x40);
		}
		
	}

	assemble_basic_body(mstack, all_slices, csvfile, 0);
	
	return rc;
}



