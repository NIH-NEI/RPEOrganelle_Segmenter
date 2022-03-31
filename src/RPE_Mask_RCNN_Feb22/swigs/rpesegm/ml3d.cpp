#include "ml3d.h"

void AssemblerML::load_flat_slices(std::vector<std::vector<int>> & all_flat_particles)
{
	for (int z=0; z<d; z++) {
		if (size_t(z) >= all_flat_particles.size()) break;
		std::vector<int> & flat_particles = all_flat_particles[z];
		std::vector<SliceML> & particles = all_slices[z];
		
		size_t iStart = 0;
		while (iStart < flat_particles.size()) {
			int hlen = flat_particles[iStart++];
			particles.resize(particles.size()+1);
			SliceML& ptc = particles[particles.size()-1];
			for (int i=0; i<hlen; i++) {
				ptc.fill.push_back(HSeg(flat_particles[iStart], flat_particles[iStart+1], flat_particles[iStart+2]));
				iStart += 3;
			}
			ptc.area = ptc.update_from_fill();
			if (ptc.area < 50) {
				particles.resize(particles.size()-1);
			} else {
				ptc.z = z;
				ptc.idx = int(particles.size()-1);
			}
		}
		
		// std::cout << "z=" << z << "; #slices=" << particles.size() << std::endl;
	}

}

// De-duping based on 'nearby chain score'
void AssemblerML::de_dupe()
{
	for (int z=0; z<d; z++) {
		std::vector<SliceML>& slices = all_slices[z];
		for (SliceML& pt : slices) {
			pt.sc = nearby_chain_score(pt);
		}
		
		for (SliceML& ptc : slices) {
			if (!ptc.valid()) continue;
			std::vector<int> dups = find_dups(ptc);
			if (dups.empty()) continue;
			bool clear_all = true;
			for (int idx : dups) {
				SliceML& pt = slices[idx];
				if (ptc.sc < pt.sc) {
					ptc.area = 0;
					clear_all = false;
					break;
				}
			}
			if (clear_all) {
				for (int idx : dups) {
					slices[idx].area = 0;
				}
			}
		}
	}
}

static bool compare_cell_sort(SortScoreML& ss1, SortScoreML& ss2)
{
	return ss1.sc > ss2.sc;
}

void AssemblerML::initial_linkage()
{
	for (int z2=d-1; z2>1; z2--) {
		std::vector<SliceML> & slices = all_slices[z2];
		
		for (int idx2=0; size_t(idx2)<slices.size(); idx2++) {
			SliceML& ptc = slices[idx2];
			if (!ptc.valid() || ptc.taken()) continue;
			
			std::vector<LinkML> chain;
			SliceML* ppt = &ptc;
			chain.push_back(LinkML(ppt->z, ppt->idx));
			for (int z1=z2-1; z1>=0; z1--) {
				if (ppt->z - z1 > MAX_GAP) break;
				int idx1 = find_below(*ppt, z1, ACCEPTABLE_IOU);
				if (idx1 < 0) continue;
				LinkML ln(z1, idx1);
				SliceML& pt = lnsl(ln);
				if (pt.taken()) break;
				ppt = &pt;
				chain.push_back(ln);
			}
			if (int(chain.size()) < MIN_SEED) continue;
			std::reverse(chain.begin(), chain.end());
			double sc = chain_score(chain);

			int ptid = int(cells.size());
			cells.resize(cells.size()+1);
			sorted_cells.push_back(SortScoreML(ptid, sc));
			Particle3D& cell = cells[ptid];
			cell.fills.resize(d);
			
			cell.bnd.set2d(ptc.bnd);
			cell.bnd.zmin = chain[0].z;
			cell.bnd.zmax = z2;
			for (LinkML& ln : chain) {
				SliceML& pt = lnsl(ln);
				pt.ptid = ptid;
				cell.add_slice(pt, ln.z);
			}
		}
		
	}
	
	for (Particle3D& cell : cells) {
		cell.update_from_fill();
	}
	std::sort(sorted_cells.begin(), sorted_cells.end(), compare_cell_sort);
}

static void frame_seq(std::vector<int> &res, int zlo, int zhi)
{
	int zmid = (zlo + zhi) / 2;
	if (zmid <= zlo || zmid >= zhi) return;
	res.push_back(zmid);
	if (zmid - zlo > 1) {
		frame_seq(res, zlo, zmid);
	}
	if (zhi - zmid > 1) {
		frame_seq(res, zmid, zhi);
	}
}

void AssemblerML::fill_gaps()
{
	mstack.fill(0);
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x20, 1);
	}
	for (SortScoreML& ss : sorted_cells) {
		Particle3D& cell = cells[ss.ptid];
		
		int zlo = -1, zhi = -1;
		for (int z=cell.bnd.zmin; z<=cell.bnd.zmax; z++) {
			std::vector<HSeg>& fill = cell.fills[z];
			if (fill.empty()) continue;
			Raster8 msk = mstack.getPlane(z);
			msk.paintParticleFillInto(fill, 0x60, 0);
			Boundary bnd = fill_boundary(fill);
			msk.rescanParticleFill(bnd, fill, 0x60);
			bnd.expand(1);
			msk.clip(bnd, 1);
			msk.expandBordersInto(bnd, 0x60, 0, 0x50, HOOD_SIZE_MOORE);
			msk.replaceColor(bnd, 0x60, 0x40);
			msk.paintParticleFill(fill, 0x80);
			
			if (zlo < 0) zlo = z;
			zhi = z;
			if (zhi - zlo > 1) {
				// std::cout << "ptid=" << ss.ptid << " sc=" << ss.sc << " zlo=" << zlo << " zgap=" << (zhi-zlo-1) << std::endl;
				std::vector<int> mf;
				frame_seq(mf, zlo, zhi);
				for (int zmiss : mf) {
					interpolate_frame(cell, zmiss);
				}
			}
			zlo = zhi;
		}
		
		//mstack.paintParticle(cell, 0x80);
	}
}

void AssemblerML::smooth_cells()
{
	mstack.fill(0x80);
	for (Particle3D& cell : cells) {
		mstack.paintParticle(cell, 0);
	}
	
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		msk.filterParticles(0x80, 0x40, 50, 0);
		msk.replaceColor(0x40, 0x80);
		msk.filterParticles(0, 0x40, 50, 0x80);
		msk.replaceColor(0x40, 0);
		msk.fillBorder(0x80, 1);
	}

	// Remove small filaments and protrusions
	mstack.expandBorders(0x80, 0, 0xC0, HOOD3D_26);
	mstack.expandBorders(0, 0xC0, 0x40, HOOD3D_6);
	mstack.replaceColor(0x40, 0);
	mstack.expandBorders(0, 0xC0, 0x40, HOOD3D_6);
	mstack.replaceColor(0x40, 0);
	mstack.replaceColor(0xC0, 0x80);
	
	mstack.rescanShrunkParticles(cells, 0, 0x50);
	mstack.expandBorders(0x80, 0, 0xFF);
}

void AssemblerML::fix_borders_actin()
{
	mstack.fill(0);

	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		int npt = 0;
		for (Particle3D& cell : cells) {
			if (cell.fills[z].empty()) continue;
			++npt;
			msk.paintParticleFill(cell.fills[z], 0xC0);
		}
		if (!npt) continue;
		
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
			for (SortScoreML& ss : sorted_cells) {
				Particle3D& cell = cells[ss.ptid];
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				Boundary bnd = fill_boundary(fill);
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.paintParticleFill(fill, 0xD0);
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

				msk.rescanParticleFill(bnd, fill, 0xD0);
				msk.paintParticleFill(fill, 0xC0);
			}
		}
	}
}

void AssemblerML::fix_borders_dna()
{
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(0x10, 1);
		
		int npt = 0;
		// Filter out "phantom slices"
		for (Particle3D& cell : cells) {
			std::vector<HSeg>& fill = cell.fills[z];
			if (fill.empty()) continue;
			double ff = double(msk.countColor(fill, 0xFF)) / fill_area(fill);
			if (ff < DNA_PHANTOM_FILL) {
				fill.clear();
			} else {
				++npt;
				msk.paintParticleFill(fill, 0xC0);
			}
		}
		if (!npt) continue;
		
		for (int pass=0; pass<2; pass++) {
			int nbsz = (pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE;
			for (SortScoreML& ss : sorted_cells) {
				Particle3D& cell = cells[ss.ptid];
				std::vector<HSeg>& fill = cell.fills[z];
				if (fill.empty()) continue;
				Boundary bnd = fill_boundary(fill);
				bnd.expand(1);
				msk.clip(bnd, 1);
				msk.paintParticleFill(fill, 0xD0);
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

				msk.rescanParticleFill(bnd, fill, 0xD0);
				msk.paintParticleFill(fill, 0xC0);
			}
		}
	}

	for (Particle3D& cell : cells)
		cell.update_from_fill();
}

int AssemblerML::filter_cells(int min_height, double min_sc)
{
	int res = 0;
	sorted_cells.resize(cells.size());
	for (int j=0; size_t(j)<cells.size(); j++) {
		cells[j].update_from_fill();
		sorted_cells[j].ptid = j;
		sorted_cells[j].sc = cells[j].iou_score(MAX_GAP);
	}
	std::sort(sorted_cells.begin(), sorted_cells.end(), compare_cell_sort);
	for (SortScoreML& ss : sorted_cells) {
		Particle3D& cell = cells[ss.ptid];
		int height = cell.bnd.zmax - cell.bnd.zmin + 1;
		if (height < min_height || ss.sc < min_sc) {
			// std::cout << "clr: h=" << height << " sc=" << ss.sc << std::endl;
			cell.clear();
		} else {
			++res;
		}
	}
	//std::cout << "Validated #cells=" << res << std::endl;
	
	mstack.fill(0x80);
	for (Particle3D& cell : cells) {
		if (!cell.empty())
			mstack.paintParticle(cell, 0);
	}
	mstack.expandBorders(0x80, 0, 0xFF);

	return res;
}


// Protected

std::vector<int> AssemblerML::find_dups(SliceML& ptc)
{
	std::vector<int> res;
	std::vector<SliceML> & slices = all_slices[ptc.z];
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		if (idx == ptc.idx) continue;
		SliceML& pt = slices[idx];
		if (!pt.valid() || !ptc.bnd.intersects(pt.bnd)) continue;
		double o = double(ptc.overlay_area(pt));
		if (o/ptc.area >= MIN_DUP_OVL || o/pt.area >= MIN_DUP_OVL) {
			res.push_back(idx);
		}
	}

	return res;
}
int AssemblerML::find_below(SliceML& ptc, int z, double good_iou)
{
	if (good_iou <= 0.) good_iou = GOOD_IOU;
	int res = -1;
	double bsf = 0.;
	std::vector<SliceML> & slices = all_slices[z];
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		SliceML& pt = slices[idx];
		if (!pt.valid() || !ptc.bnd.intersects(pt.bnd)) continue;
		int o = ptc.overlay_area(pt);
		double iou = double(o) / (ptc.area + pt.area - o);
		if (iou > bsf) {
			bsf = iou;
			res = idx;
		}
	}
	
	return bsf >= good_iou ? res : -1;
}
std::vector<int> AssemblerML::find_above(SliceML& ptc, int z)
{
	std::vector<int> res;
	std::vector<SliceML> & slices = all_slices[z];
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		SliceML& pt = slices[idx];
		if (!pt.valid() || !ptc.bnd.intersects(pt.bnd)) continue;
		int o = ptc.overlay_area(pt);
		double iou = double(o) / (ptc.area + pt.area - o);
		if (iou >= GOOD_IOU) {
			res.push_back(idx);
		}
	}
	
	return res;
}

double AssemblerML::chain_score(std::vector<LinkML>& chain)
{
	double res = 0.;
	
	for (int lni0=0; lni0<int(chain.size())-1; lni0++) {
		SliceML& pt0 = lnsl(chain[lni0]);
		for (int lni1=lni0+1; size_t(lni1)<chain.size(); lni1++) {
			LinkML& ln = chain[lni1];
			if (ln.z - pt0.z > MAX_GAP) continue;
			SliceML& pt1 = lnsl(ln);
			int o = pt0.overlay_area(pt1);
			double iou = double(o) / (pt0.area + pt1.area - o);
			res += iou;
		}
	}
	
	return res;
}

double AssemblerML::fwd_chain_score(SliceML& ptc, std::vector<LinkML>& chain, int z1)
{
	if (z1 - ptc.z > MAX_GAP || z1 >= d-1) return chain_score(chain);
	int z2 = z1 + 1;
	size_t c_sz = chain.size();

	double res = fwd_chain_score(ptc, chain, z2);
	LinkML& ln1 = chain[chain.size() - 1];
	SliceML& pt1 = lnsl(ln1);
	std::vector<int> zlist = find_above(pt1, z2);
	if (!zlist.empty()) {
		chain.resize(c_sz + 1);
		LinkML& ln = chain[chain.size() - 1];
		ln.z = z2;
		for (int idx : zlist) {
			ln.idx = idx;
			double sc = fwd_chain_score(ptc, chain, z2);
			if (sc > res) res = sc;
		}
		chain.resize(c_sz);
	}

	return res;
}

double AssemblerML::nearby_chain_score(SliceML& ptc)
{
	std::vector<LinkML> chain;
	
	SliceML *ppt = &ptc;
	chain.push_back(LinkML(ppt->z, ppt->idx));
	for (int dz=1; dz<=MAX_GAP; dz++) {
		int z1 = ptc.z - dz;
		if (z1 < 0) break;
		int idx1 = find_below(*ppt, z1);
		if (idx1 >= 0) {
			ppt = &(all_slices[z1][idx1]);
			chain.push_back(LinkML(ppt->z, ppt->idx));
		}
	}
	std::reverse(chain.begin(), chain.end());
	
	return fwd_chain_score(ptc, chain, ptc.z+1);
}

void AssemblerML::interpolate_frame(Particle3D& cell, int zmiss)
{
	Raster8 msk = mstack.getPlane(zmiss);
	int zlo, zhi;
	for (zlo = zmiss-1; zlo>=0; zlo--) {
		if (!cell.fills[zlo].empty()) break;
	}
	for (zhi = zmiss+1; zhi<d; zhi++) {
		if (!cell.fills[zhi].empty()) break;
	}
	if (zlo < 0 || zhi >= d) return;
	
	std::vector<HSeg>& fill_lo = cell.fills[zlo];
	std::vector<HSeg>& fill = cell.fills[zmiss];
	std::vector<HSeg>& fill_hi = cell.fills[zhi];
	Boundary bnd = fill_boundary(fill_lo);
	Boundary bnd_hi = fill_boundary(fill_hi);
	bnd.combo(bnd_hi);
	
	msk.paintParticleFillInto(fill_lo, 0x44, 0);
	msk.paintParticleFillInto(fill_hi, 0x60, 0x44);
	msk.paintParticleFillInto(fill_hi, 0x44, 0);
	
	for (int y0=bnd.ymin+1; y0<=bnd.ymax-1; y0++) {
		unsigned char *b = msk.scanLine(y0);
		for (int x0=bnd.xmin+1; x0<=bnd.xmax-1; x0++) {
			if (b[x0] != 0x44) continue;
			int nc = 0, nb = 0;
			for (int j=1; j<HOOD_SIZE_FATCROSS; j++) {
				unsigned char c = msk.value(x0+hood_pts[j].dx, y0+hood_pts[j].dy);
				if (c == 0x60) {
					++nc;
					break;
				}
				if (c != 0x44 && c != 0x48 && c != 0x30) {
					++nb;
				}
			}
			b[x0] = (nc > 0 || nb == 0) ? 0x48 : 0x30;
		}
	}
	msk.replaceColor(0x30, 0);
	msk.replaceColor(0x48, 0x60);
	
	msk.rescanParticleFill(bnd, fill, 0x60);
	bnd.expand(1);
	msk.clip(bnd, 1);
	msk.expandBordersInto(bnd, 0x60, 0, 0x50, HOOD_SIZE_MOORE);
	msk.replaceColor(bnd, 0x60, 0x40);
	msk.paintParticleFill(fill, 0xFF);

}





