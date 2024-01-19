#include "assembly3d.h"

//----------- Assembler3d -----------

int Assembler3d::detect_slices()
{
	int np = 0;
	for (int z=0; z<d; z++) {
		std::vector<AsmSlice> & slices = all_slices[z];
		Raster8 msk = mstack.getPlane(z);
		msk.fillBorder(bkg, 1);
		AsmSlice::detect_all(slices, msk, fgd, nvc, z);
		np += int(slices.size());
// std::cout << "plane " << z << " : detected " << slices.size() << " particles" << std::endl;
	}
	return np;
}

int Assembler3d::initial_linkage()
{
	for (int z0=d-2; z0>=0; z0--) {
		std::vector<AsmSlice> & slices = all_slices[z0];
		for (int idx=0; size_t(idx)<slices.size(); idx++) {
			AsmSlice &ptc = slices[idx];
			if (!ptc.valid()) continue;
			int other_idx = -1;
			int z;
			for (int j=1; j<=3; j++) {
				z = z0 + j;
				if (z >= d) break;
				other_idx = find_full_match(ptc, z);
				if (other_idx >= 0) break;
			}
			if (other_idx < 0) continue;
			AsmSlice &other = all_slices[z][other_idx];
			if (other.orig_idx < 0) {
				int oidx = int(origins.size());
				origins.resize(origins.size()+1);
				AsmOrigin & orig = origins[oidx];
				orig.init(d, oidx);
				orig.zmap[z] = other_idx;
				orig.zmap[z0] = idx;
				ptc.orig_idx = other.orig_idx = orig.idx;
			} else {
				AsmOrigin & orig = origins[other.orig_idx];
				orig.zmap[z0] = idx;
				ptc.orig_idx = orig.idx;
			}
		}
	}
	int valid_cnt = 0;
	for (AsmOrigin & orig : origins) {
		if (!orig.valid) continue;
		if (orig.count() >= 3)
			++valid_cnt;
		else
			invalidate_orig(orig);
	}
	
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		std::vector<AsmSlice> & slices = all_slices[z];
		for (AsmOrigin & orig : origins) {
			if (!orig.valid) continue;
			int idx = orig.zmap[z];
			if (idx >= 0)
				msk.paintParticle(slices[idx], vc);
		}
	}
	
	return valid_cnt;
}

int Assembler3d::absorb_chips()
{
	int cnt = 0;
	std::vector<int> ichips;
	unsigned char mainc = tmpc + 1;
	unsigned char chipc = tmpc + 2;
	for (AsmOrigin &orig : origins) {
		if (!orig.valid) continue;
		if (calc_average(orig) == 0) continue;

		for (int z=orig.zmin; z<=orig.zmax; z++) {
			int idx = orig.zmap[z];
			if (idx < 0) continue;
			int nc = find_chips(ichips, orig, z);
			if (nc == 0) continue;
			cnt += nc;
			Raster8 msk = mstack.getPlane(z);
			std::vector<AsmSlice> & slices = all_slices[z];
			AsmSlice &mptc = slices[idx];
			msk.paintParticle(mptc, mainc);
			Boundary bnd = mptc.bnd;
			for (int ic : ichips) {
				AsmSlice &cptc = slices[ic];
				msk.paintParticle(cptc, chipc);
				bnd.combo(cptc.bnd);
				cptc.invalidate();
			}
			
			for (int y0=bnd.ymin+1; y0<=bnd.ymax-1; y0++) {
				unsigned char *p = msk.scanLine(y0);
				for (int x0=bnd.xmin+1; x0<=bnd.xmax-1; x0++) {
					if (p[x0] != bkg) continue;
					bool touch_main = false;
					bool touch_chip = false;
					bool touch_other = false;
					for (int j=1; j<HOOD_SIZE_MOORE; j++) {
						unsigned char c = msk.value(x0+hood_pts[j].dx, y0+hood_pts[j].dy);
						if (c == tmpc) continue;
						if (c == mainc) touch_main = true;
						else if (c == chipc) touch_chip = true;
						else if (c != bkg) {
							touch_other = true;
							break;
						}
					}
					if (touch_other) continue;
					if (touch_main && touch_chip)
						p[x0] = tmpc;
				}
			}
			msk.replaceColor(bnd, tmpc, mainc);
			msk.replaceColor(bnd, chipc, mainc);
			
			mptc.bnd = bnd;
			msk.rescanParticle(mptc, mainc);
			msk.paintParticle(mptc, vc);
		}
	}
	return cnt;
}

int Assembler3d::find_conflicts()
{
	for (AsmOrigin& orig : origins) {
		if (!orig.valid) continue;
		estimate_zbounds(orig);
		orig.conflicts.clear();
	}
	if (origins.size() < 2) return 0;
	
	for (int idx=0; size_t(idx)<origins.size()-1; idx++) {
		AsmOrigin& orig = origins[idx];
		if (!orig.valid) continue;
		for (int oidx=idx+1; size_t(oidx)<origins.size(); oidx++) {
			AsmOrigin& other = origins[oidx];
			if (!other.valid) continue;
			if (orig.zmax < other.zmin || other.zmax < orig.zmin) continue;
			int minarea = orig.area < other.area ? orig.area : other.area;
			int ovl = orig.overlay_area(other);
			if (ovl <= int(minarea * conf_ovl)) continue;
			orig.conflicts.insert(oidx);
			other.conflicts.insert(idx);
		}
	}
	
	int cnt = 0;
	for (AsmOrigin& orig : origins) {
		if (orig.valid && orig.conflicts.size() > 0) ++cnt;
	}
	return cnt;
}

int Assembler3d::resolve_nonconflicts()
{
	int cnt = 0;
	for (AsmOrigin& orig : origins) {
		if (!orig.valid) continue;
		if (orig.conflicts.size() > 0) continue;
		++cnt;
		assemble(orig);
	}
	return cnt;
}

/*
static void print_conflicts(std::vector<int>& conflicted) {
	for (int idx : conflicted) {
		std::cout << " " << idx;
	}
	std::cout << std::endl;
}
static void show_conflicts(AsmOrigin& orig) {
	if (!orig.valid || orig.conflicts.empty())
		return;
	std::cout << "idx=" << orig.idx << " :";
	for (int i : orig.conflicts) {
		std::cout << " " << i;
	}
	std::cout << std::endl;
}
*/

int Assembler3d::resolve_conflicts()
{
	if (origins.size() == 0) return 0;
	
// Pass 1
	for (int idx=0; size_t(idx)<origins.size(); idx++) {
		AsmOrigin& first = origins[idx];
		if (!first.valid || first.conflicts.empty()) continue;
		
		std::vector<int> conflicted;
		conflicted.push_back(idx);
		for (int i : first.conflicts)
			conflicted.push_back(i);
		
		if (conflicted.size() > 3) {
			resolve_m(conflicted);
		}
		
		if (conflicted.size() == 2) {
			AsmOrigin& second = origins[conflicted[1]];
			if (second.conflicts.size() <= 1 && resolve_2(first, second)) {
				continue;
			}
		} else if (conflicted.size() == 3) {
			AsmOrigin& second = origins[conflicted[1]];
			AsmOrigin& third = origins[conflicted[2]];
			if (second.conflicts.size() <= 2 && third.conflicts.size() <= 2 &&
					resolve_3(first, second, third)) {
				continue;
			}
		}
		
	}
	
// Pass 2
	for (int idx=0; size_t(idx)<origins.size(); idx++) {
		AsmOrigin& first = origins[idx];
		if (!first.valid || first.conflicts.empty()) continue;
		
		std::vector<int> conflicted;
		conflicted.push_back(idx);
		
		for (int oidx : first.conflicts) {
			AsmOrigin& other = origins[oidx];
			if (other.valid)
				conflicted.push_back(oidx);
			else
				first.conflicts.erase(oidx);
		}
		
		if (conflicted.size() == 2) {
			AsmOrigin& second = origins[conflicted[1]];
			if (resolve_2(first, second)) {
				continue;
			}
		} else if (conflicted.size() == 3) {
			AsmOrigin& second = origins[conflicted[1]];
			AsmOrigin& third = origins[conflicted[2]];
			if (resolve_3(first, second, third)) {
				continue;
			}
		}
		
		// Give up
		assemble(first);
		first.conflicts.clear();
	}
	
	return conflict_count();
}

static bool find_zrange(std::vector<int> & zmap, int *pzmin, int *pzmax, int maxgap=2)
{
	int zmin=-1, zmax=-1, dz=-1;
	int in_range = 0;
	int z0=-1, z1=-1;
	for (int z=0; size_t(z)<zmap.size(); z++) {
		if (zmap[z] >= 0) {
			in_range = maxgap;
			if (z0 < 0) z0 = z;
			z1 = z;
			continue;
		}
		if (in_range == 0) continue;
		--in_range;
		if (in_range > 0) continue;
		int _dz = z1 - z0;
		if (_dz >= dz) {
			zmin = z0;
			zmax = z1;
			dz = _dz;
		}
		z0 = z1 = -1;
	}
	if (z0 >= 0 && (z1 - z0) >= dz) {
		zmin = z0;
		zmax = z1;
		dz = z1 - z0;
	}
	if (dz >= 0) {
		*pzmin = zmin;
		*pzmax = zmax;
		return true;
	}
	return false;
}

int Assembler3d::bread_crumbs()
{
	Raster8 msk(w, h, (unsigned char *)(dat.buf));
	msk.fill(0);
	for (int z=0; z<d; z++) {
		for (AsmSlice& ptc : all_slices[z]) {
			if (ptc.valid() && !ptc.taken())
				msk.paintParticleADD(ptc, 1);
		}
	}
	for (long long lx=0; lx<msk.len; lx++) {
		if (msk.buf[lx] > 8) msk.buf[lx] = 8;
	}
	msk.fillBorder(0xE0, 2);
	std::vector<Slice> crumbs;
	
	unsigned char tmp_fg = 0xA0;
	unsigned char tmp_bk = 0xA1;
	unsigned char cur_fg, next_fg;
	for (int lvl=8; lvl>3; lvl--) {
		cur_fg = (unsigned char)(lvl);
		next_fg = cur_fg - 1;
		
		Particle tmp_ptc;
		for (Slice &ptc : crumbs) {
			for (HSeg &hs : ptc.fill) {
				unsigned char *p = msk.scanLine(hs.y);
				for (int x=hs.xl; x<=hs.xr; x++) {
					for (int j=1; j<HOOD_SIZE_NEUMANN; j++) {
						int x0 = x + hood_pts[j].dx;
						int y0 = hs.y + hood_pts[j].dy;
						if (msk.value(x0, y0) == cur_fg) {
							msk.detectParticle(tmp_ptc, x0, y0, tmp_bk);
							msk.paintParticle(tmp_ptc, next_fg);
						}
					}
				}
			}
		}
		
		for (int y0=1; y0<h-1; y0++) {
			unsigned char *p = msk.scanLine(y0);
			for (int x0=1; x0<w-1; x0++) {
				if (p[x0] != cur_fg) continue;
				
				crumbs.resize(crumbs.size()+1);
				Slice &ptc = crumbs[crumbs.size()-1];
				ptc.area = msk.detectParticle(ptc, x0, y0, tmp_fg);
				if (ptc.area < 100) {
					msk.paintParticle(ptc, next_fg);
					crumbs.resize(crumbs.size()-1);
				}
			}
		}

		msk.replaceColor(cur_fg, next_fg);
	}
	msk.replaceColor(next_fg, 0x40);
	for (int pass=0; pass<4; pass++) {
		for (Slice &ptc : crumbs) {
			msk.expandParticle(ptc, tmp_fg, tmp_bk, 0x40,
				(pass & 1) ? HOOD_SIZE_NEUMANN : HOOD_SIZE_MOORE);
		}
	}
	
	int cnt = 0;
	
	std::vector<int> zmap(d);
	for (Slice &crumb : crumbs) {
		crumb.area = fill_area(crumb.fill);
		crumb.bnd = fill_boundary(crumb.fill);
		int crumb_ovl = int(crumb.area * conf_ovl);
		for (int z=0; z<d; z++) {
			zmap[z] = -1;
			std::vector<AsmSlice> & slices = all_slices[z];
			for (int idx=0; size_t(idx)<slices.size(); idx++) {
				AsmSlice& ptc = slices[idx];
				if (!ptc.valid() || ptc.taken()) continue;
				int ovl = crumb.overlay_area(ptc);
				if (ovl > crumb_ovl) {
					zmap[z] = idx;
					break;
				}
			}
		}
		
		int zmin, zmax;
		if (!find_zrange(zmap, &zmin, &zmax)) continue;
		if (zmax-zmin < 3) continue;
		
		int oidx = int(origins.size());
		origins.resize(origins.size() + 1);
		AsmOrigin &orig = origins[oidx];
		orig.init(d, oidx);
		for (int z=zmin; z<=zmax; z++) {
			int idx = zmap[z];
			if (idx < 0) continue;
			AsmSlice& ptc = all_slices[z][idx];
			ptc.orig_idx = oidx;
			orig.zmap[z] = idx;
			Raster8 msk = mstack.getPlane(z);
			msk.paintParticle(ptc, vc);
		}
		estimate_zbounds(orig);
		assemble(orig);
		++cnt;
	}

	
	return cnt;
}

void Assembler3d::toParticle3D(AsmOrigin &orig, Particle3D &pt)
{
	pt.fills.resize(orig.zmap.size());
	pt.x0 = pt.y0 = pt.z0 = 0;

	Boundary3D bnd(0,0,0,0,0,0);
	bool is_first = true;
	for (int z=0; size_t(z)<orig.zmap.size(); z++) {
		std::vector<HSeg> &tfill = pt.fills[z];
		tfill.resize(0);
		int idx = orig.zmap[z];
		if (idx < 0) {
			continue;
		}
		AsmSlice &ptc = all_slices[z][idx];
		if (ptc.fill.empty()) continue;
		bnd.zmax = z;
		tfill.resize(ptc.fill.size());
		for (size_t i=0; i<ptc.fill.size(); i++) {
			HSeg &hs = ptc.fill[i];
			if (is_first) {
				is_first = false;
				pt.x0 = hs.xl;
				pt.y0 = hs.y;
				pt.z0 = z;
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
			tfill[i] = hs;
		}
		tfill.shrink_to_fit();
	}
	pt.bnd = bnd;
}


void Assembler3d::invalidate_orig(AsmOrigin &orig)
{
	for (int z=0; z<d; z++) {
		int idx = orig.zmap[z];
		if (idx < 0) continue;
		orig.zmap[z] = -1;
		AsmSlice &ptc = all_slices[z][idx];
		ptc.orig_idx = -1;
		Raster8 msk = mstack.getPlane(z);
		msk.paintParticle(ptc, nvc);
	}
	orig.invalidate();
}

int Assembler3d::find_full_match(AsmSlice &ptc, int z)
{
	if (ptc.z == z) return -1;
	int ptc_ovl = int(ptc.area * match_ovl);
	std::vector<AsmSlice> & slices = all_slices[z];
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		AsmSlice &other = slices[idx];
		if (!other.valid()) continue;
		int ovl = ptc.overlay_area(other);
		if (ovl < ptc_ovl) continue;
		if (ovl >= int(other.area * match_ovl))
			return idx;
	}
	return -1;
}

int Assembler3d::calc_average(AsmOrigin &orig)
{
	std::vector<AsmSlice *> pslices;
	int zmin=-1, zmax=-1;
	for (int z=0; z<d; z++) {
		int idx = orig.zmap[z];
		if (idx < 0) continue;
		pslices.push_back(&(all_slices[z][idx]));
		if (zmin < 0) zmin = z;
		zmax = z;
	}
	if (pslices.size() < 2) {
		invalidate_orig(orig);
		return 0;
	}
	orig.zmin = zmin;
	orig.zmax = zmax;
	Boundary bnd = pslices[0]->bnd;
	for (int i=1; size_t(i)<pslices.size(); i++)
		bnd.combo(pslices[i]->bnd);
	dat.fill(bnd, 0);
	for (AsmSlice *pptc : pslices)
		dat.paintParticleADD(*pptc, 1);
	orig.bnd = bnd;
	unsigned short tc = pslices.size() < 5 ? 2 : 3;
	orig.area = dat.rescanParticleThresh(orig, tc);
	return int(pslices.size());
}

int Assembler3d::find_chips(std::vector<int> &ichips, AsmOrigin &orig, int z)
{
	ichips.clear();
	if (calc_average(orig) == 0)
		return 0;
	std::vector<AsmSlice> & slices = all_slices[z];
	Raster8 msk = mstack.getPlane(z);
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		AsmSlice &ptc = slices[idx];
		if (!ptc.valid() || ptc.taken()) continue;
		int ovl = orig.overlay_area(ptc);
		if (ovl >= int(ptc.area * chip_ovl)) {
			ichips.push_back(idx);
		}
	}
	return int(ichips.size());
}

int Assembler3d::estimate_zbounds(AsmOrigin &orig)
{
	int nvs = calc_average(orig);
	if (nvs == 0) return 0;
	int zmin = orig.zmin;
	if (zmin > 0) {
		for (--zmin; zmin>=0; zmin--) {
			Raster8 msk = mstack.getPlane(zmin);
			int ba = msk.countColor(orig, bkg);
			if (ba * 2 > orig.area) break;
		}
	}
	if (zmin < 0) zmin = 0;
	orig.zmin = zmin;
	int zmax = orig.zmax;
	if (zmax < d-1) {
		for (++zmax; zmax<d; zmax++) {
			Raster8 msk = mstack.getPlane(zmax);
			int ba = msk.countColor(orig, bkg);
			if (ba * 2 > orig.area) break;
		}
	}
	if (zmax >= d) zmax = d - 1;
	orig.zmax = zmax;
	return nvs;
}

int Assembler3d::assemble_slice(AsmOrigin &orig, int z)
{
	unsigned char tmp_bk = tmpc + 1;
	unsigned char tmp_fg = tmpc + 2;
	unsigned char tmp_new = tmpc + 3;

	std::vector<AsmSlice> & slices = all_slices[z];
	Raster8 msk = mstack.getPlane(z);
	
	Boundary bnd = orig.bnd;
	std::vector<int> bkg_idxes;
	int njoined = 0;
	int cidx = orig.zmap[z];
	if (cidx >= 0) {
		AsmSlice& ptc = slices[cidx];
		bnd.combo(ptc.bnd);
		msk.paintParticle(ptc, tmp_fg);
		++njoined;
	}
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		AsmSlice& ptc = slices[idx];
		if (!ptc.valid() || ptc.taken()) continue;
		int ovl = orig.overlay_area(ptc);
		if (ovl <= int(ptc.area * conf_ovl)) continue;
		if ((ovl >= int(asm_ovl * ptc.area)) ||
			(ovl >= int(0.5 * ptc.area) && ptc.area*4 < orig.area)) {
			bnd.combo(ptc.bnd);
			msk.paintParticle(ptc, tmp_fg);
			++njoined;
			if (cidx < 0) {
				ptc.orig_idx = orig.idx;
				orig.zmap[z] = cidx = idx;
			} else {
				ptc.invalidate();
			}
			continue;
		}
	
		msk.paintParticle(ptc, tmpc);
		++njoined;
		bkg_idxes.push_back(idx);
	
		msk.paintParticleInto(orig, tmp_fg, tmpc);
		msk.paintParticleInto(ptc, tmp_bk, tmpc);
		msk.expandBordersInto(bnd, tmp_fg, tmp_bk, tmpc);
		msk.expandBordersInto(bnd, tmp_fg, tmp_bk, tmpc);
	}
	
	if (cidx < 0 && msk.countColor(orig, tmp_fg) > 0) {
		orig.zmap[z] = cidx = int(slices.size());
		slices.resize(slices.size()+1);
		AsmSlice& ptc = slices[cidx];
		ptc.z = z;
		ptc.orig_idx = orig.idx;
	}
	if (cidx >= 0) {
		AsmSlice& ptc = slices[cidx];
		ptc.bnd = bnd;
		ptc.area = msk.rescanParticle(ptc, tmp_fg);
		if (njoined > 1) {
			msk.removeFalseWalls(ptc, vc, bkg, tmpc+4);
		} else {
			msk.paintParticle(ptc, vc);
		}
	}
	for (int i : bkg_idxes) {
		AsmSlice& ptc = slices[i];
		msk.paintParticleInto(ptc, tmpc, tmp_bk);
		ptc.area = msk.rescanParticle(ptc, tmpc);
		msk.paintParticle(ptc, nvc);
	}
	
	return cidx;
}

void Assembler3d::assemble(AsmOrigin &orig)
{
	int zlo = -1, zhi = -1;
	for (int z=0; z<d; z++) {
		if (orig.zmap[z] < 0) continue;
		zhi = z;
		if (zlo < 0) zlo = z;
	}
	if (zlo < 0) {
		invalidate_orig(orig);
		return;
	}
	for (int z=zlo; z<=zhi; z++) {
		int idx = assemble_slice(orig, z);
	}
	int zmin = orig.zmin - 2;
	if (zmin < 0) zmin = 0;
	int zmax = orig.zmax + 2;
	if (zmax >= d) zmax = d - 1;
	int cutfac = 10;
	for (int z=zlo-1; z>=zmin; z--) {
		int idx = assemble_slice(orig, z);
		if (idx < 0) break;
		orig.zmin = z;
		AsmSlice& ptc = all_slices[z][idx];
		if (ptc.area * cutfac < orig.area) break;
	}
	for (int z=zhi+1; z<=zmax; z++) {
		int idx = assemble_slice(orig, z);
		if (idx < 0) break;
		orig.zmax = z;
		AsmSlice& ptc = all_slices[z][idx];
		if (ptc.area * cutfac < orig.area) break;
	}
}

bool Assembler3d::resolve_2(AsmOrigin& first, AsmOrigin& second)
{
	int ovl;
	if (first.zmin > second.zmax || second.zmin > first.zmax)
		goto resolved2;
	
	if (can_merge_2(first, second)) {
		merge_2(first, second);
		goto resolved1;
	}
	
	ovl = first.overlay_area(second);
	if (ovl < int(first.area * conf_ovl)) {
		separate_2(second, first);
	} else if (ovl < int(second.area * conf_ovl)) {
		separate_2(first, second);
	} else {
		decouple_2(first, second);
	}
	
resolved2:
	second.conflicts.erase(first.idx);
	assemble(second);
resolved1:
	first.conflicts.erase(second.idx);
	assemble(first);
	return true;
}

bool Assembler3d::can_merge_2(AsmOrigin& first, AsmOrigin& second)
{
	int minovl = int(first.area * merge_ovl);
	int minovl2 = int(second.area * merge_ovl);
	if (minovl2 < minovl) minovl = minovl2;
	int ovl = first.overlay_area(second);
	if (ovl < minovl)
		return false;
	
	int zlo, zhi;
	int z11 = first.zfirst();
	int z12 = first.zlast() + 1;
	int z21 = second.zfirst();
	int z22 = second.zlast() + 1;
	if (z12 < z21) {
		zlo = z12;
		zhi = z21;
	} else if (z22 < z11) {
		zlo = z22;
		zhi = z11;
	} else
		return true;
	
	for (int z=zlo; z<zhi; z++) {
		Raster8 msk = mstack.getPlane(z);
		ovl = first.area - msk.countColor(first, bkg);
		if (ovl < minovl)
			return false;
		ovl = second.area - msk.countColor(second, bkg);
		if (ovl < minovl)
			return false;
	}
	return true;
}

void Assembler3d::merge_2(AsmOrigin& first, AsmOrigin& second)
{
	for (int z=0; z<d; z++) {
		int idx2 = second.zmap[z];
		if (idx2 < 0) continue;
		AsmSlice& ptc2 = all_slices[z][idx2];
		int idx1 = first.zmap[z];
		if (idx1 < 0) {
			first.zmap[z] = idx2;
			ptc2.orig_idx = first.idx;
		} else {
			AsmSlice& ptc1 = all_slices[z][idx1];
			Raster8 msk = mstack.getPlane(z);
			ptc1.area = msk.mergeParticles(ptc1, ptc2, vc, bkg, tmpc);
			ptc2.invalidate();
		}
		second.zmap[z] = -1;
	}
	second.invalidate();
	estimate_zbounds(first);
}

void Assembler3d::separate_slice_2(AsmOrigin& first, AsmOrigin& second, int z)
{
	std::vector<AsmSlice> & slices = all_slices[z];
	int idx2 = second.zmap[z];
	int ovl = first.overlay_area(slices[idx2]);
	if (ovl < int(first.area * conf_ovl)) return;

	Raster8 msk = mstack.getPlane(z);
	slices[idx2].orig_idx = -1;
	msk.paintParticle(slices[idx2], nvc);
	int idx1 = assemble_slice(first, z);
	if (idx1 >= 0) {
		msk.paintParticle(slices[idx1], vc);
	}
	if (slices[idx2].valid()) {
		slices[idx2].orig_idx = second.idx;
		msk.paintParticle(slices[idx2], vc);
	} else {
		second.zmap[z] = -1;
	}
}

void Assembler3d::separate_2(AsmOrigin& first, AsmOrigin& second)
{
	for (int z=0; z<d; z++) {
		if (second.zmap[z] >= 0)
			separate_slice_2(first, second, z);
	}
}

void Assembler3d::decouple_2(AsmOrigin& first, AsmOrigin& second)
{
	int zmid1 = (first.zfirst() + first.zlast()) / 2;
	int zmid2 = (second.zfirst() + second.zlast()) / 2;
	if (zmid1 > zmid2) {
		decouple_2(second, first);
		return;
	}
	int zmin = first.zmin;
	if (zmin < second.zmin) zmin = second.zmin;
	int zmax = first.zmax;
	if (zmax > second.zmax) zmax = second.zmax;
	if (zmin > zmax) return;

	bool first_over = false;
	for (int z=zmin; z<=zmax; z++) {
		Raster8 msk = mstack.getPlane(z);
		if (!first_over) {
			int bka = msk.countColor(first, bkg);
			double occ1 = double(first.area - bka) / first.area;
			bka = msk.countColor(second, bkg);
			double occ2 = double(second.area - bka) / second.area;
			if (occ1+0.3 < occ2) {
				first_over = true;
			} else {
				if (second.zmap[z] >= 0)
					separate_slice_2(first, second, z);
				continue;
			}
		}
		if (first.zmap[z] >= 0)
			separate_slice_2(second, first, z);
	}
}

bool Assembler3d::resolve_3(AsmOrigin& first, AsmOrigin& second, AsmOrigin &third)
{
	if (second.area > first.area && second.area > third.area) {
		return resolve_3(second, first, third);
	} else if (third.area > first.area && third.area > second.area) {
		return resolve_3(third, first, second);
	}
	
	int oa12 = first.overlay_area(second);
	int oa13 = first.overlay_area(third);
	int oa23 = second.overlay_area(third);
	if ( ((oa12+oa13) >= int(first.area * match_ovl)) &&
		 (oa23 < int(second.area * conf_ovl)) &&
		 (oa23 < int(third.area * conf_ovl)) ) {
		separate_3(first, second, third);
		return true;
	}

	if (can_merge_2(first, second)) {
		first.conflicts.erase(second.idx);
		merge_2(first, second);
		return resolve_2(first, third);
	} else if (can_merge_2(first, third)) {
		first.conflicts.erase(third.idx);
		merge_2(first, third);
		return resolve_2(first, second);
	}
	
	if (oa12 < int(second.area * conf_ovl)) {
		separate_2(second, first);
	} else {
		decouple_2(first, second);
	}
	if (oa13 < int(third.area * conf_ovl)) {
		separate_2(third, first);
	} else {
		decouple_2(first, third);
	}

	bool rc = resolve_2(second, third);
	second.conflicts.erase(first.idx);
	first.conflicts.erase(second.idx);
	third.conflicts.erase(first.idx);
	first.conflicts.erase(third.idx);
	assemble(first);
	return rc;
}

void Assembler3d::separate_3(AsmOrigin& first, AsmOrigin& second, AsmOrigin &third)
{
	unsigned char tmp_fg1 = tmpc + 1;
	unsigned char tmp_fg2 = tmpc + 2;
	unsigned char tmp_fg3 = tmpc + 3;
	
	for (int z=0; z<d; z++) {
		int idx1 = first.zmap[z];
		if (idx1 < 0) continue;
		
		std::vector<AsmSlice> & slices = all_slices[z];
		Raster8 msk = mstack.getPlane(z);
		slices[idx1].orig_idx = -1;
		first.zmap[z] = -1;
		Boundary bnd = slices[idx1].bnd;
		msk.paintParticle(slices[idx1], tmp_fg1);
		
		int idx2 = second.zmap[z];
		if (idx2 < 0) {
			idx2 = int(slices.size());
			slices.resize(slices.size()+1);
			second.zmap[z] = idx2;
			slices[idx2].z = z;
			slices[idx2].orig_idx = second.idx;
		} else {
			bnd.combo(slices[idx2].bnd);
			msk.paintParticle(slices[idx2], tmp_fg2);
		}
		
		int idx3 = third.zmap[z];
		if (idx3 < 0) {
			idx3 = int(slices.size());
			slices.resize(slices.size()+1);
			third.zmap[z] = idx3;
			slices[idx3].z = z;
			slices[idx3].orig_idx = third.idx;
		} else {
			bnd.combo(slices[idx3].bnd);
			msk.paintParticle(slices[idx3], tmp_fg3);
		}
		
		msk.paintParticleInto(second, tmp_fg2, tmp_fg1);
		bnd.combo(second.bnd);
		msk.paintParticleInto(third, tmp_fg3, tmp_fg1);
		bnd.combo(third.bnd);
		
		for (int k=0; k<4; k++) {
			msk.expandBordersInto(bnd, tmp_fg2, tmp_fg1, tmpc);
			msk.expandBordersInto(bnd, tmp_fg3, tmp_fg1, tmpc);
		}

		AsmSlice & ptc1 = slices[idx1];
		ptc1.bnd = bnd;
		ptc1.area = msk.rescanParticle(ptc1, tmp_fg1);
		if (ptc1.area > 0)
			msk.paintParticle(ptc1, nvc);
		else
			ptc1.invalidate();

		AsmSlice & ptc2 = slices[idx2];
		ptc2.bnd = bnd;
		ptc2.area = msk.rescanParticle(ptc2, tmp_fg2);
		msk.paintParticle(ptc2, 0xA0);

		AsmSlice & ptc3 = slices[idx3];
		ptc3.bnd = bnd;
		ptc3.area = msk.rescanParticle(ptc3, tmp_fg3);
		msk.paintParticle(ptc3, 0xC0);
	}
	
	second.conflicts.erase(first.idx);
	second.conflicts.erase(third.idx);
	third.conflicts.erase(first.idx);
	third.conflicts.erase(second.idx);
	
	first.invalidate();
	
	estimate_zbounds(second);
	assemble(second);
	estimate_zbounds(third);
	assemble(third);
}

void Assembler3d::resolve_m(std::vector<int>& conflicted)
{
	std::vector<AsmOrigin *> corigs;
	for (int i : conflicted) {
		corigs.push_back(&(origins[i]));
	}
	
	for (int i=0; size_t(i)<corigs.size()-1; i++) {
		AsmOrigin & first = *(corigs[i]);
		if (!first.valid) continue;
		for (int j=i+1; size_t(j)<corigs.size(); j++) {
			AsmOrigin & second = *(corigs[j]);
			if (!second.valid) continue;
			if (can_merge_2(first, second)) {
				merge_2(first, second);
				for (AsmOrigin *porig : corigs)
					porig->conflicts.erase(second.idx);
			}
		}
	}
	
	conflicted.clear();
	for (AsmOrigin *porig : corigs) {
		if (!porig->valid) continue;
		if (porig->conflicts.empty())
			assemble(*porig);
		else
			conflicted.push_back(porig->idx);
	}
}

//----------- ActinAssembler3d -----------

int ActinAssembler3d::read_z01_data(const char *z01_csvfile)
{
	z01slices.clear();

	CsvCellDataReader rdr(z01_csvfile);
	if (rdr.eof || !rdr.is_2d) return -1;

	int cur_id = -1;
	AsmSlice *pptc = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			z01slices.resize(z01slices.size()+1);
			pptc = & z01slices[z01slices.size()-1];
			pptc->z = 0;
			cur_id = id;
		}
		if (!pptc) continue;
		pptc->fill.push_back(HSeg(y, xl, xr));
	}
	
	for (AsmSlice& ptc : z01slices) {
		ptc.fill.shrink_to_fit();
		ptc.area = ptc.update_from_fill();
	}
	return cur_id;
}

int ActinAssembler3d::z01_slice_linkage()
{
	if (z01slices.empty()) return 0;
	
	std::vector<std::vector<SliceId>> z01map;
	z01map.resize(z01slices.size());
	
	for (int z=d-1; z>=0; z--) {
		std::vector<AsmSlice> & slices = all_slices[z];
		for (int z01_idx=0; size_t(z01_idx)<z01slices.size(); z01_idx++) {
			AsmSlice & ptcz = z01slices[z01_idx];
			int minovl = int(ptcz.area * match_ovl);
			for (int slice_idx=0; size_t(slice_idx)<slices.size(); slice_idx++) {
				AsmSlice & ptc = slices[slice_idx];
				if (!ptc.valid() || ptc.taken()) continue;
				int ovl = ptcz.overlay_area(ptc);
				if (ovl<minovl || ovl <int(ptc.area * match_ovl)) continue;
				z01map[z01_idx].push_back(SliceId(z, slice_idx));
			}
		}
	}
	
	int cnt = 0;
	for (int z01_idx=0; size_t(z01_idx)<z01slices.size(); z01_idx++) {
		std::vector<SliceId>& zslices = z01map[z01_idx];
		if (zslices.size() < 2) continue;
		for (size_t j=1; j<zslices.size(); j++) {
			if (zslices[j].z - zslices[j-1].z > 5) {
				zslices.resize(j);
				break;
			}
		}
		if (zslices.size() < 2) continue;
		AsmSlice & ptcz = z01slices[z01_idx];
		
		int orig_idx = int(origins.size());
		origins.resize(origins.size()+1);
		AsmOrigin & orig = origins[orig_idx];
		orig.init(d, orig_idx);
		
		int zkey = -1;
		double ovl_rat = 0.;
		for (SliceId & si : zslices) {
			AsmSlice & ptc = all_slices[si.z][si.idx];
			if (!ptc.valid() || ptc.taken()) continue;
			int ovl = ptcz.overlay_area(ptc);
			double rat = double(ovl)*2. / (ptcz.area + ptc.area);
			if (zkey < 0 || ovl_rat < rat) {
				zkey = si.z;
				ovl_rat = rat;
			}
			orig.zmap[si.z] = si.idx;
			ptc.orig_idx = orig_idx;
		}
		if (zkey < 0) {
			// Very strange if this ever happens...
			invalidate_orig(orig);
			origins.resize(origins.size()-1);
			continue;
		}
		orig.zkey = zkey;
		orig.zidx = z01_idx;
		++cnt;
	}
	
	for (int z=0; z<d; z++) {
		Raster8 msk = mstack.getPlane(z);
		for (AsmOrigin& orig : origins) {
			if (orig.zmap[z] < 0) continue;
			AsmSlice & ptc = all_slices[z][orig.zmap[z]];
			msk.paintParticle(ptc, vc);
		}
	}
	
	return cnt;
}

int ActinAssembler3d::read_dna_data(const char *dna_csvfile)
{
	nuclei.clear();
	CsvCellDataReader rdr(dna_csvfile);
	if (rdr.eof || !rdr.is_3d) return -1;

	int cur_id = -1;
	Nucleus *pnuc = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (z<0 || z>=d || y<0 || y>=h || xl<0 || xr>=w || xl>xr) continue;
		if (id != cur_id) {
			nuclei.resize(nuclei.size()+1);
			pnuc = & nuclei[nuclei.size()-1];
			pnuc->fills.resize(d);
			cur_id = id;
		}
		if (!pnuc) continue;
		pnuc->fills[z].push_back(HSeg(y, xl, xr));
	}
	
	for (int idx=0; size_t(idx)<nuclei.size(); idx++) {
		Nucleus &nuc = nuclei[idx];
		for (int z=0; z<d; z++) nuc.fills[z].shrink_to_fit();
		nuc.fills.shrink_to_fit();
		nuc.update_from_fill();
	}
	return cur_id;
}

NucKey ActinAssembler3d::key_for_slice(int z, int slice_idx)
{
	NucKey nk;
	AsmSlice& ptc = all_slices[z][slice_idx];
	for (int nuc_idx=0; size_t(nuc_idx)<nuclei.size(); nuc_idx++) {
		Nucleus & nuc = nuclei[nuc_idx];
		int ovl = nuc.overlay_area2d(ptc, z);
		if (ovl == 0) continue;
		if (ovl > int(fill_area(nuc.fills[z]) * nuc_ovl))
			nk.add_nuc_idx(nuc_idx);
	}
	
	return nk;
}

void ActinAssembler3d::split_slicenuc_2(int sidx)
{
	unsigned char tmp_fg1 = tmpc + 1;
	unsigned char tmp_fg2 = tmpc + 2;

	NucKey nk0 = sliceNucMap[sidx].nk;
	int nuc_idx1 = nk0.idxs[0];
	int nuc_idx2 = nk0.idxs[1];
	NucKey nk1, nk2;
	nk1.add_nuc_idx(nuc_idx1);
	nk2.add_nuc_idx(nuc_idx2);
	int tidx1 = find_or_create_slicenuc(nk1);
	int tidx2 = find_or_create_slicenuc(nk2);
	SliceNucMatch& snm0 = sliceNucMap[sidx];
	SliceNucMatch& tnm1 = sliceNucMap[tidx1];
	SliceNucMatch& tnm2 = sliceNucMap[tidx2];
	Nucleus& nuc1 = nuclei[nuc_idx1];
	Nucleus& nuc2 = nuclei[nuc_idx2];
	for (SliceId& si : snm0.sliceMap) {
		int z = si.z;
		if (nuc1.fills[z].empty() || nuc2.fills[z].empty()) continue;
		std::vector<AsmSlice> & slices = all_slices[z];
		int slice_idx1 = si.idx;
		int slice_idx2 = int(slices.size());
		slices.resize(slices.size() + 1);
		AsmSlice & ptc1 = slices[slice_idx1];
		AsmSlice & ptc2 = slices[slice_idx2];
		ptc2.z = z;
		double xc1, yc1, xc2, yc2;
		fill_centroid(nuc1.fills[z], &xc1, &yc1);
		fill_centroid(nuc2.fills[z], &xc2, &yc2);
		Raster8 msk = mstack.getPlane(z);
		Boundary bnd = ptc1.bnd;
		msk.paintParticle(ptc1, tmpc);
		msk.paintParticleFillInto(nuc1.fills[z], tmp_fg1, tmpc);
		msk.paintParticleFillInto(nuc2.fills[z], tmp_fg2, tmpc);
		for (int y=bnd.ymin; y<=bnd.ymax; y++) {
			unsigned char *p = msk.scanLine(y);
			for (int x=bnd.xmin; x<=bnd.xmax; x++) {
				if (p[x] != tmpc) continue;
				double dx = xc1 - x;
				double dy = yc1 - y;
				double d1 = dx*dx + dy*dy;
				dx = xc2 - x;
				dy = yc2 - y;
				double d2 = dx*dx + dy*dy;
				p[x] = (d1 < d2) ? tmp_fg1 : tmp_fg2;
			}
		}
		ptc1.bnd = bnd;
		msk.rescanParticle(ptc1, tmp_fg1);
		ptc2.bnd = bnd;
		msk.rescanParticle(ptc2, tmp_fg2);
		
		tnm1.add_slice(z, slice_idx1);
		tnm2.add_slice(z, slice_idx2);
		
		msk.paintParticle(ptc1, nvc);
		msk.paintParticle(ptc2, nvc);
	}
	snm0.invalidate();
}

int ActinAssembler3d::nuc_slice_linkage()
{
	sliceNucMap.clear();
	for (int z=d-1; z>=0; z--) {
		std::vector<AsmSlice> & slices = all_slices[z];
		for (int slice_idx=0; size_t(slice_idx)<slices.size(); slice_idx++) {
			NucKey nk = key_for_slice(z, slice_idx);
			if (nk.len < 1 || nk.len > 4) continue;
			int i = find_or_create_slicenuc(nk);
			sliceNucMap[i].add_slice(z, slice_idx);
		}
	}
	
	std::vector<int> to_invalidate;
	for (int l=4; l>=2; l--) {
		for (int sidx=0; size_t(sidx)<sliceNucMap.size(); sidx++) {
			SliceNucMatch & snm = sliceNucMap[sidx];
			if (snm.nk.len != l || snm.len() < 3) continue;
			for (int tidx=0; size_t(tidx)<sliceNucMap.size(); tidx++) {
				if (tidx == sidx) continue;
				SliceNucMatch & tnm = sliceNucMap[tidx];
				if (!snm.conflicts(tnm)) continue;
				if (l > 2) {
					if (tnm.nk.len >= 2)
						to_invalidate.push_back(tidx);
					to_invalidate.push_back(sidx);
					continue;
				} else if (l == 2) {
					if (tnm.nk.len > 2) {
						to_invalidate.push_back(tidx);
					} else if (tnm.nk.len == 2) {
						if (tnm.len() <= snm.len())
							to_invalidate.push_back(tidx);
						else
							to_invalidate.push_back(sidx);
					} else if (tnm.nk.len == 1) {
						if (tnm.len() < 3)
							snm.merge(tnm);
						else
							to_invalidate.push_back(sidx);
					}
				}
			}
		}
	}
	
	for (int sidx : to_invalidate) {
		SliceNucMatch & snm = sliceNucMap[sidx];
		if (snm.len() == 0) continue;
		if (snm.nk.len == 2 && snm.len() > 4) {
			split_slicenuc_2(sidx);
			continue;
		}
		sliceNucMap[sidx].invalidate();
	}
	
/* DEBUG
	int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
	for (SliceNucMatch & snm : sliceNucMap) {
		if (snm.len() < 3) continue;
		int n = snm.nk.len;
		if (n==1) ++cnt1;
		else if (n==2) ++cnt2;
		else if (n==3) ++cnt3;
		else if (n==4) ++cnt4;
	}
std::cout << "cnt1=" << cnt1 << " cnt2=" << cnt2 << " cnt3=" << cnt3 << " cnt4=" << cnt4 << std::endl;
*/
	
	for (SliceNucMatch & snm : sliceNucMap) {
		if (snm.len() < 3) continue;
		
		int orig_idx = -1;
		for (SliceId& si : snm.sliceMap) {
			AsmSlice& ptc = all_slices[si.z][si.idx];
			if (!ptc.valid()) continue;
			if (ptc.orig_idx >= 0) {
				orig_idx = ptc.orig_idx;
				break;
			}
		}
		
		if (orig_idx < 0) {
			orig_idx = int(origins.size());
			origins.resize(origins.size()+1);
			origins[orig_idx].init(d, orig_idx);
		}
		AsmOrigin & orig = origins[orig_idx];
		
		for (SliceId& si : snm.sliceMap) {
			int z = si.z;
			AsmSlice& ptc = all_slices[z][si.idx];
			if (!ptc.valid() || ptc.taken()) continue;
			Raster8 msk = mstack.getPlane(z);
			if (orig.zmap[z] >= 0) {
				AsmSlice& optc = all_slices[z][orig.zmap[z]];
				optc.area = msk.mergeParticles(optc, ptc, vc, bkg, tmpc);
				ptc.invalidate();
			} else {
				orig.zmap[z] = si.idx;
				ptc.orig_idx = orig_idx;
				msk.paintParticle(ptc, vc);
			}
		}
	}
	
	return int(origins.size());
}

int ActinAssembler3d::initial_linkage()
{
	Assembler3d::initial_linkage();
	int rc = 0;
	for (AsmOrigin& orig : origins) {
		if (!orig.valid) continue;
		if (orig.zlast() - orig.zfirst() < 4) {
			invalidate_orig(orig);
		} else ++rc;
	}
	return rc;
}

int ActinAssembler3d::nuc_cell_linkage()
{
	nucCellMap.clear();
	for (AsmOrigin & orig : origins) {
		if (orig.valid) estimate_zbounds(orig);
	}
	return int(nucCellMap.size());
}

int ActinAssembler3d::resolve_nuclear_conflicts()
{
	int cnt = 0;
	for (AsmOrigin & orig : origins)
		if (orig.valid) ++cnt;
	
	for (int nuc_idx=0; size_t(nuc_idx)<nuclei.size(); nuc_idx++) {
		// Nucleus &nuc = nuclei[nuc_idx];
		std::vector<int> cellids = find_cells_bynuc(nuc_idx);
		if (cellids.size() < 2) continue;
		AsmOrigin & first = origins[cellids[0]];
		if (!first.valid) {
			invalidate_nuccell(cellids[0]);
			continue;
		}
		for (int i=nuc_idx+1; size_t(i)<cellids.size(); i++) {
			AsmOrigin & second = origins[cellids[i]];
			if (second.valid)
				merge_2(first, second);
			invalidate_nuccell(cellids[i]);
		}
	}
	
	for (AsmOrigin & orig : origins)
		if (orig.valid) --cnt;
	return cnt;
}

int ActinAssembler3d::resolve_conflicts()
{
	int cnt = 0;
	
	for (int idx1=0; size_t(idx1)<origins.size()-1; idx1++) {
		AsmOrigin & first = origins[idx1];
		if (!first.valid || count_nuclei_bycell(idx1)>0) continue;
		for (int idx2=idx1+1; size_t(idx2)<origins.size(); idx2++) {
			AsmOrigin & second = origins[idx2];
			if (!second.valid || count_nuclei_bycell(idx2)>0) continue;
			int minarea = first.area < second.area ? first.area : second.area;
			int ovl = first.overlay_area(second);
			if (ovl <= int(minarea * conf_ovl)) continue;
			if (can_merge_2(first, second)) {
				merge_2(first, second);
				invalidate_nuccell(idx2);
				++cnt;
			}
		}
	}
	
	return cnt;
}

int ActinAssembler3d::fill_missing_slices()
{
	int cnt = 0;
	
	for (int pass=0; pass<2; pass++) {
		for (int orig_idx=0; size_t(orig_idx)<origins.size(); orig_idx++) {
			AsmOrigin & orig = origins[orig_idx];
			if (!orig.valid) continue;
			std::vector<int> nucids = find_nuclei_bycell(orig_idx);
			if (pass==0 && nucids.size() != 1) continue;
			if (pass==1 && nucids.size() == 1) continue;
			
			std::vector<AsmSlice> seeds;
			for (int z=orig.zmin; z<=orig.zmax; z++) {
				if (orig.zmap[z] >= 0) continue;
				seeds.resize(seeds.size()+1);
				AsmSlice & seed = seeds[seeds.size()-1];
				calculate_seed(seed, orig, z, nucids);
			}
			if (seeds.empty()) continue;
			++cnt;
			
			int zmin = orig.zmin;
			int zmax = orig.zmax;
			for (AsmSlice & seed : seeds) {
				int idx = slice_from_seed(seed, orig);
				int z = seed.z;
				if (z < zmin || z > zmax) continue;
				AsmSlice& ptc = all_slices[z][idx];
				double circ = fill_circularity(ptc.fill);
				if (circ < min_slice_circ) {
					if (z < orig.zfirst()) zmin = z;
					else if (z > orig.zlast()) zmax = z;
					orig.zmap[z] = -1;
					ptc.orig_idx = -1;
					Raster8 msk = mstack.getPlane(z);
					msk.paintParticle(ptc, nvc);
				}
			}
		}
	}
	
	return cnt;
}

int ActinAssembler3d::estimate_zbounds(AsmOrigin &orig)
{
	int nvs = calc_average(orig);
	if (nvs == 0) return 0;
	
	for (int nuc_idx=0; size_t(nuc_idx)<nuclei.size(); nuc_idx++) {
		Nucleus &nuc = nuclei[nuc_idx];
		if (!nuc.bnd.intersects2d(orig.bnd)) continue;
		long long ovol = overlay_volume(nuc, orig);
		double ovl_fract = double(ovol) / nuc.vol;
		if (ovl_fract < nuc_vol_min_fract) continue; 
		update_nuc_cell_map(orig.idx, nuc_idx, ovl_fract);
		
		if (orig.zmin > nuc.z_lo) orig.zmin = nuc.z_lo;
		if (orig.zmax < nuc.z_hi) orig.zmax = nuc.z_hi;
	}
	
	return nvs;
}

long long ActinAssembler3d::overlay_volume(Nucleus &nuc, AsmOrigin &orig)
{
	long long ovol = 0;
	for (int z=0; z<d; z++) {
		if (nuc.fills[z].empty() || orig.zmap[z] < 0) continue;
		AsmSlice &ptc = all_slices[z][orig.zmap[z]];
		ovol += ptc.overlay_area(nuc.fills[z]);
	}
	return ovol;
}

void ActinAssembler3d::update_nuc_cell_map(int cell_idx, int nuc_idx, double ovl_fract)
{
	for (NucCellMatch ncm : nucCellMap) {
		if (ncm.cell_idx == cell_idx && ncm.nuc_idx == nuc_idx) {
			ncm.ovl_fract = ovl_fract;
			ncm.valid = true;
			return;
		}
	}
	nucCellMap.push_back(NucCellMatch(cell_idx, nuc_idx, ovl_fract));
}

void ActinAssembler3d::calculate_seed(AsmSlice &seed, AsmOrigin &orig, int z, std::vector<int> &nucids)
{
	seed.z = z;
	seed.orig_idx = orig.idx;
	
	if ((orig.zkey >= 0) && (orig.zidx >= 0) && (z >= orig.zkey-7) && (z <= orig.zkey+5)) {
		// Use Z01 slice as template seed if available and close enough along z-axis
		seed.copyfrom(z01slices[orig.zidx]);
		seed.area = fill_area(seed.fill);
		return;
	}

	Boundary bnd = orig.bnd;
	for (int nuc_idx : nucids) {
		Nucleus &nuc = nuclei[nuc_idx];
		if (nuc.fills[z].empty()) continue;
		Boundary nbnd = fill_boundary(nuc.fills[z]);
		bnd.combo(nbnd);
	}
	std::vector<AsmSlice *> nb_slices;
	for (int dz=1; dz<5; dz++) {
		int zz = z + dz;
		if (zz <= orig.zmax && orig.zmap[zz] >= 0) {
			AsmSlice & ptc = all_slices[zz][orig.zmap[zz]];
			bnd.combo(ptc.bnd);
			nb_slices.push_back(&ptc);
		}
		zz = z - dz;
		if (zz >= orig.zmin && orig.zmap[zz] >= 0) {
			AsmSlice & ptc = all_slices[zz][orig.zmap[zz]];
			bnd.combo(ptc.bnd);
			nb_slices.push_back(&ptc);
		}
		if (nb_slices.size() >= 4) break;
	}
	
	int tc = int(nb_slices.size() / 2) + 1;
	
	dat.fill(bnd, 0);
	for (int nuc_idx : nucids) {
		dat.paintParticleFill(nuclei[nuc_idx].fills[z], 2);
	}
	dat.paintParticleADD(orig, 1);
	for (AsmSlice * pptc : nb_slices)
		dat.paintParticleADD(*pptc, 1);

	seed.bnd = bnd;
	seed.area = dat.rescanParticleThresh(seed, (unsigned short)(tc));
}

int ActinAssembler3d::slice_from_seed(AsmSlice& seed, AsmOrigin& orig)
{
	int z = seed.z;
	int confl_area = int(seed.area * conf_ovl);
	
	unsigned char tmp_bk = tmpc + 1;
	unsigned char tmp_fg = tmpc + 2;
	std::vector<AsmSlice> & slices = all_slices[z];
	Raster8 msk = mstack.getPlane(z);
	Boundary bnd = seed.bnd;
	
	for (int idx=0; size_t(idx)<slices.size(); idx++) {
		AsmSlice & ptc = slices[idx];
		if (!ptc.valid() || ptc.taken()) continue;
		int ovl = seed.overlay_area(ptc);
		if (ovl > int(ptc.area * seed_ovl)) {
			bnd.combo(ptc.bnd);
			msk.paintParticle(ptc, tmp_fg);
			ptc.invalidate();
		} else if (ovl > confl_area) {
			msk.paintParticle(ptc, tmp_bk);
			msk.paintParticleInto(seed, tmp_fg, tmp_bk);
			msk.rescanParticle(ptc, tmp_bk);
			msk.paintParticle(ptc, nvc);
		}
	}
	
	int cidx = orig.zmap[z];
	if (cidx >= 0) {
		msk.paintParticle(slices[cidx], tmp_fg);
		bnd.combo(slices[cidx].bnd);
	} else {
		cidx = int(slices.size());
		slices.resize(slices.size() + 1);
		AsmSlice & ptc = slices[cidx];
		ptc.z = z;
		ptc.orig_idx = orig.idx;
		orig.zmap[z] = cidx;
	}
	msk.paintParticleInto(seed, tmp_fg, bkg);
	AsmSlice & nptc = slices[cidx];
	nptc.bnd = bnd;
	msk.rescanParticle(nptc, tmp_fg);

	msk.paintParticle(nptc, vc);
	
	return cidx;
}


//----------------------------- TopBottomDetector --------------------------------

int TopBottomDetector::read_dna_data(const char *dna_csvfile)
{
	nuclei.clear();
	CsvCellDataReader rdr(dna_csvfile);
	if (rdr.eof) return -1;

	int cur_id = -1;
	Nucleus *pnuc = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (z<0 || z>=d || y<0 || y>=h || xl<0 || xr>=w) continue;
		if (id != cur_id) {
			nuclei.resize(nuclei.size()+1);
			pnuc = & nuclei[nuclei.size()-1];
			pnuc->fills.resize(d);
			cur_id = id;
		}
		if (!pnuc) continue;
		pnuc->fills[z].push_back(HSeg(y, xl, xr));
	}
	
	for (int idx=0; size_t(idx)<nuclei.size(); idx++) {
		Nucleus &nuc = nuclei[idx];
		nuc.idx = idx;
		for (int z=0; z<d; z++) nuc.fills[z].shrink_to_fit();
		nuc.fills.shrink_to_fit();
		nuc.update_from_fill(0.);
	}
	return cur_id;
}

int TopBottomDetector::read_actin_data(const char *actin_csvfile)
{
	cells.clear();
	CsvCellDataReader rdr(actin_csvfile);
	if (rdr.eof) return -1;

	int cur_id = -1;
	Cell *pcell = NULL;
	int id, z, y, xl, xr;
	while(rdr.read_hs(&id, &z, &y, &xl, &xr) > 0) {
		if (z<0 || z>=d || y<0 || y>=h || xl<0 || xr>=w) continue;
		if (id != cur_id) {
			cells.resize(cells.size()+1);
			pcell = & cells[cells.size()-1];
			pcell->fills.resize(d);
			cur_id = id;
		}
		if (!pcell) continue;
		pcell->fills[z].push_back(HSeg(y, xl, xr));
	}
	
	for (int idx=0; size_t(idx)<cells.size(); idx++) {
		Cell &cell = cells[idx];
		cell.idx = idx;
		for (int z=0; z<d; z++) cell.fills[z].shrink_to_fit();
		cell.fills.shrink_to_fit();
		cell.update_from_fill();
	}
	return cur_id;
}

void TopBottomDetector::map_cells_to_nuclei()
{
	nucCellMap.clear();
//	int n0=0, n1=0, n2=0, nn=0;
	for (int cell_idx=0; size_t(cell_idx)<cells.size(); cell_idx++) {
		Cell &cell = cells[cell_idx];
		cell.nnucs = 0;
		for (int nuc_idx=0; size_t(nuc_idx)<nuclei.size(); nuc_idx++) {
			Nucleus &nuc = nuclei[nuc_idx];
			if (!cell.bnd.intersects(nuc.bnd)) continue;
			long long ovl_vol = cell.overlay_volume(nuc);
			double ovl_fract = double(ovl_vol) / nuc.vol;
			if (ovl_fract >= min_nuc_fract) {
				nucCellMap.push_back(NucCellMatch(cell_idx, nuc_idx, ovl_fract));
				++cell.nnucs;
			}
		}
		// if (cell.nnucs == 0) ++n0;
		// else if (cell.nnucs == 1) ++n1;
		// else if (cell.nnucs == 2) ++n2;
		// else if (cell.nnucs > 2) ++nn;
	}
// std::cout << "cells with no nuclei: " << n0 << std::endl;
// std::cout << "cells with 1 nucleus: " << n1 << std::endl;
// std::cout << "cells with 2 nuclei: " << n2 << std::endl;
// std::cout << "cells with >2 nuclei: " << nn << std::endl;
}

void TopBottomDetector::actin_cell_limits()
{
	for (int cell_id=0; size_t(cell_id)<cells.size(); cell_id++) {
		Cell &cell = cells[cell_id];
		std::vector<Nucleus *> cnucs = find_cell_nuclei(cell_id);
		Slice inner, outer;
		int z_lo, z_hi, z0, z1;
		int minoutv;
		bool first;
		int maxdz = (cell.bnd.zmax - cell.bnd.zmin + 2) / 3;
		
		// Cell top (Actin)
		z_hi = cell.bnd.zmax;
		z_lo = z_hi - 2;
		if (z_lo < cell.bnd.zmin) z_lo = cell.bnd.zmin;
		combo_slice(outer, inner, cell, z_lo, z_hi);
		z0 = cell.bnd.zmax;
		for (Nucleus *pnuc : cnucs) {
			if (z0 < pnuc->bnd.zmax) z0 = pnuc->bnd.zmax;
		}
		first = true;
		minoutv = 5;
		z1 = z0;
		for (int z=z0; z<mstack.d; z++) {
			Raster8 msk = mstack.getPlane(z);
			int outv = msk.countColor(outer, 0x80) + 2*msk.countColor(outer, 0xFF);
			if (first) {
				minoutv = 5 + outv/4;
				first = false;
			} else
				if (outv < minoutv) break;
			int inv = msk.countColor(inner, 0x80) + 2*msk.countColor(inner, 0xFF);
			z1 = z;
			if ((z1 - z0) >= maxdz) break;
		}
		cell.z_hi = z1;
		
		// Cell bottom (Actin)
		z_lo = cell.bnd.zmin;
		z_hi = z_lo + 2;
		if (z_hi > cell.bnd.zmax) z_hi = cell.bnd.zmax;
		combo_slice(outer, inner, cell, z_lo, z_hi);
		z0 = cell.bnd.zmin;
		for (Nucleus *pnuc : cnucs) {
			if (z0 > pnuc->bnd.zmin) z0 = pnuc->bnd.zmin;
		}
		first = true;
		minoutv = 5;
		z1 = z0;
		for (int z=z0; z>=0; z--) {
			Raster8 msk = mstack.getPlane(z);
			int outv = msk.countColor(outer, 0x80) + 2*msk.countColor(outer, 0xFF);
			if (first) {
				minoutv = 5 + outv/4;
				first = false;
			} else
				if (outv < minoutv) break;
			int inv = msk.countColor(inner, 0x80) + 2*msk.countColor(inner, 0xFF);
			z1 = z;
			if ((z0 - z1) >= maxdz) break;
		}
		cell.z_lo = z1;
	}
}

void TopBottomDetector::detect_top_bottom()
{
	int z_lo, z_hi;
	int z0, z1;

	mstack.fill(0);
	for (Cell &cell : cells) {
		mstack.paintParticle(cell, 0xFF);
	}

	for (int cell_id=0; size_t(cell_id)<cells.size(); cell_id++) {
		Cell &cell = cells[cell_id];
		std::vector<Nucleus *> cnucs = find_cell_nuclei(cell_id);
		Slice inner, outer;
		
		z0 = cell.bnd.zmax + 1;
		z1 = cell.z_hi;
		if (z0 <= z1) {
			// Extrapolate cell top
			z_hi = cell.bnd.zmax;
			z_lo = z_hi - 2;
			if (z_lo < cell.bnd.zmin) z_lo = cell.bnd.zmin;
			combo_slice(outer, inner, cell, z_lo, z_hi);
			
			int zrad = z1 - z0 + 1;
			if (zrad < minzrad) zrad = minzrad;
			int iz = 0;
			for (int z=z0; z<=z1; z++, iz++) {
				Raster8 msk = mstack.getPlane(z);
				double gamma = 1. - pow(double(iz)/zrad, 2.);
				int ar = int(outer.area * gamma);
				if (iz == 0) {
					msk.paintParticleInto(outer, 0x80, 0);
					inner.bnd = outer.bnd;
				} else {
					msk.paintParticleInto(inner, 0x80, 0);
				}
				for (Nucleus *pnuc : cnucs) {
					msk.paintParticleFillInto(pnuc->fills[z], 0x80, 0);
				}
				msk.rescanParticle(inner, 0x80);
				if (inner.area < 50 || fill_circularity(inner.fill) < min_slice_circ) {
					msk.paintParticle(inner, 0);
					break;
				}
				int nbsz = HOOD_SIZE_NEUMANN;
				int iter = 0;
				while (inner.area > ar ) {
					msk.shrinkParticle(inner, 0x80, 0x50, 0x51, nbsz);
					inner.area = fill_area(inner.fill);
					msk.replaceColor(outer.bnd, 0x51, 0);
					nbsz = (nbsz==HOOD_SIZE_NEUMANN) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
					++ iter;
					if (iter > 100) {
						std::cout << "too many shrink iterations!" << std::endl;
						break;
					}
				}
				msk.paintParticle(inner, 0xA0);
				cell.fills[z] = inner.fill;
				cell.bnd.zmax = z;
			}
			
			// double gamma = 1. - pow(double(ii)/delta, 2.);
		}
		
		z0 = cell.bnd.zmin - 1;
		z1 = cell.z_lo;
		if (z0 >= z1) {
			// Extrapolate cell bottom
			z_lo = cell.bnd.zmin;
			z_hi = z_lo + 2;
			if (z_hi > cell.bnd.zmax) z_hi = cell.bnd.zmax;
			combo_slice(outer, inner, cell, z_lo, z_hi);

			int zrad = z0 - z1 + 1;
			if (zrad < minzrad) zrad = minzrad;
			int iz = 0;
			for (int z=z0; z>=z1; z--, iz++) {
				Raster8 msk = mstack.getPlane(z);
				double gamma = 1. - pow(double(iz)/zrad, 2.);
				int ar = int(outer.area * gamma);
				if (iz == 0) {
					msk.paintParticleInto(outer, 0x80, 0);
					inner.bnd = outer.bnd;
				} else {
					msk.paintParticleInto(inner, 0x80, 0);
				}
				for (Nucleus *pnuc : cnucs) {
					msk.paintParticleFillInto(pnuc->fills[z], 0x80, 0);
				}
				msk.rescanParticle(inner, 0x80);
				if (inner.area < 50 || fill_circularity(inner.fill) < min_slice_circ) {
					msk.paintParticle(inner, 0);
					break;
				}
				int nbsz = HOOD_SIZE_NEUMANN;
				int iter = 0;
				while (inner.area > ar ) {
					msk.shrinkParticle(inner, 0x80, 0x50, 0x51, nbsz);
					inner.area = fill_area(inner.fill);
					msk.replaceColor(outer.bnd, 0x51, 0);
					nbsz = (nbsz==HOOD_SIZE_NEUMANN) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
					++ iter;
					if (iter > 100) {
						std::cout << "too many shrink iterations!" << std::endl;
						break;
					}
				}
				msk.paintParticle(inner, 0xA0);
				cell.fills[z] = inner.fill;
				cell.bnd.zmax = z;
			}
		}
	}
}


void TopBottomDetector::combo_slice(Slice &outer, Slice &inner, Cell &cell, int z_lo, int z_hi)
{
	Boundary bnd = cell.bnd.boundary2d();
	dat.fill(bnd, 0);
	unsigned char zt = (unsigned char)(z_hi - z_lo + 1);
	if (zt > 2) zt = 2;
	for (int z=z_lo; z<=z_hi; z++) {
		dat.paintParticleFillADD(cell.fills[z], 1);
	}
	for (int y=bnd.ymin; y<=bnd.ymax; y++) {
		unsigned char *p = dat.scanLine(y);
		for (int x=bnd.xmin; x<=bnd.xmax; x++) {
			p[x] = (p[x] >= zt) ? 0xFF : 0;
		}
	}
	outer.bnd = bnd;
	outer.area = dat.rescanParticle(outer, 0xFF);
	
	bnd = outer.bnd;
	for (int pass=0; pass<4; pass++) {
		int nbsz = (pass & 1) ? HOOD_SIZE_MOORE : HOOD_SIZE_NEUMANN;
		for (int y0=bnd.ymin; y0<=bnd.ymax; y0++) {
			unsigned char *p = dat.scanLine(y0);
			for (int x0=bnd.xmin; x0<=bnd.xmax; x0++) {
				if (p[x0] != 0xFF) continue;
				for (int j=1; j<nbsz; j++) {
					if (dat.value(x0+hood_pts[j].dx, y0+hood_pts[j].dy) == 0) {
						p[x0] = 0x80;
						break;
					}
				}
			}
		}
		dat.replaceColor(0x80, 0);
	}
	
	inner.bnd = bnd;
	inner.area = dat.rescanParticle(inner, 0xFF);
	
	dat.paintParticle(inner, 0);
}

//----------------------------- SegmentationComparator --------------------------------

void SegmentationComparator::load_cells()
{
	cells.clear();
	msk.fillBorder(bordc, 1);
	for (int y0=1; y0<msk.h-1; y0++) {
		unsigned char *p = msk.scanLine(y0);
		for (int x0=1; x0<w-1; x0++) {
			if (p[x0] != fg) continue;
			cells.resize(cells.size()+1);
			Slice &cell = cells[cells.size()-1];
			cell.area = msk.detectParticle(cell, x0, y0, tmpc);
		}
	}
}

void SegmentationComparator::load_rs_cells(const char *csvfile)
{
	rs_cells.clear();
	rmsk.fillBorder(bordc, 1);
	Boundary mbnd = rmsk.getBoundary();

	CsvReader rdr(csvfile);
	if (!rdr.next()) return;
// std::cout << "Reading: " << csvfile << std::endl;
	
	CsvRow headers = rdr.GetRow();
	rdr.setHeaders(headers);

	int xs_idx = rdr.header_index("XStart");
	int ys_idx = rdr.header_index("YStart");
	while (rdr.next()) {
		int x0 = rdr.GetInt(xs_idx, -1);
		int y0 = rdr.GetInt(ys_idx, -1);
// std::cout << x0 << "," << y0 << std::endl;
		if (!mbnd.IsInside(x0, y0)) continue;
		rs_cells.resize(rs_cells.size()+1);
		Slice &cell = rs_cells[rs_cells.size()-1];
		cell.area = rmsk.detectParticle(cell, x0, y0, tmpc);
	}
}

void SegmentationComparator::match_cells()
{
	matchmap.clear();
	for (int idx=0; size_t(idx)<cells.size(); idx++) {
		Slice &cell = cells[idx];
		int minovl = int(movl * cell.area);
		for (int rs_idx=0; size_t(rs_idx)<rs_cells.size(); rs_idx++) {
			int rs_ovl = cell.overlay_area(rs_cells[rs_idx]);
			if (rs_ovl >= minovl) {
				matchmap.push_back(SliceMatch(idx, rs_idx, rs_ovl));
				break;
			}
		}
	}
}

int SegmentationComparator::find_matches(std::vector<Slice *> & matches, int rs_idx)
{
	int res = 0;
	matches.clear();
	for (SliceMatch &sm : matchmap) {
		if (sm.rs_idx != rs_idx) continue;
		res += sm.rs_ovl;
		matches.push_back(& cells[sm.idx]);
	}
	return res;
}

