#ifndef assembly3d_h
#define assembly3d_h

#include "raster.h"
#include "csv.h"
#include <unordered_set>

struct AsmSlice : public Slice
{
	int z;
	int orig_idx;
	
	//
	AsmSlice() : Slice(), z(-1), orig_idx(-1) {}
	
	bool valid() { return z >= 0; }
	bool taken() { return orig_idx >= 0; }
	
	void invalidate() {
		z = -1;
		fill.clear();
		bnd = Boundary(0,0,0,0);
	}
	
	//
	static void detect_all(std::vector<AsmSlice> & slices, Raster8 & msk,
		unsigned char oldc, unsigned char newc, int z=-1)
	{
		for (int y0=1; y0<msk.h-1; y0++) {
			unsigned char *p = msk.scanLine(y0);
			for (int x0=1; x0<msk.w-1; x0++) {
				if (p[x0] != oldc) continue;
				slices.resize(slices.size()+1);
				AsmSlice &ptc = slices[slices.size()-1];
				ptc.area = msk.detectParticle(ptc, x0, y0, newc);
				ptc.z = z;
			}
		}
	}
};

struct AsmOrigin : public Slice
{
	int idx;
	std::vector<int> zmap;
	bool valid;
	int zmin, zmax;
	std::unordered_set<int> conflicts;
	int zkey;
	int zidx;
	//
	AsmOrigin() : Slice(), valid(true) {}
	//
	void init(int d, int _idx) {
		idx = _idx;
		zmap.resize(d);
		for (int i=0; size_t(i)<zmap.size(); i++)
			zmap[i] = -1;
		zkey = zidx = -1;
	}
	void invalidate() {
		for (int i=0; size_t(i)<zmap.size(); i++)
			zmap[i] = -1;
		conflicts.clear();
		fill.clear();
		area = 0;
		valid = false;
	}
	int count() {
		int cnt = 0;
		for (int idx : zmap)
			if (idx >= 0) ++cnt;
		return cnt;
	}
	int zfirst() {
		for (int z=0; size_t(z)<zmap.size(); z++) {
			if (zmap[z] >= 0) return z;
		}
		return -1;
	}
	int zlast() {
		for (int z = int(zmap.size()-1); z>=0; z--) {
			if (zmap[z] >= 0) return z;
		}
		return -1;
	}
};

class Assembler3d
{
public:
	//
	unsigned char bkg = 0;
	unsigned char fgd = 0xFF;
	unsigned char vc = 0xFF;
	unsigned char nvc = 0x80;
	unsigned char tmpc = 0x50;	// temp. color range: tmpc..tmpc+0xF
	//
	double match_ovl = 0.75;
	double chip_ovl = 0.75;
	double conf_ovl = 0.2;
	double asm_ovl = 0.85;
	double merge_ovl = 0.85;

	//
	int w, h, d;
	Raster3D mstack;
	Raster16 dat;
	std::vector<std::vector<AsmSlice>> all_slices;
	std::vector<AsmOrigin> origins;
	
	Assembler3d(int _w, int _h, int _d, unsigned char *mask3d) :
			w(_w), h(_h), d(_d),
			mstack(_w, _h, _d, mask3d), dat(_w, _h, NULL) {
		dat.fill(0);
		all_slices.resize(_d);
	}
	
	int conflict_count() {
		int cnt = 0;
		for (AsmOrigin & orig : origins) {
			if (orig.valid && !orig.conflicts.empty())
				++cnt;
		}
		return cnt;
	}

	int detect_slices();
	int initial_linkage();
	int absorb_chips();
	int find_conflicts();
	int resolve_nonconflicts();
	int resolve_conflicts();
	int bread_crumbs();
	void toParticle3D(AsmOrigin &orig, Particle3D &pt);
	
	//
	virtual int estimate_zbounds(AsmOrigin &orig);
	//
	void invalidate_orig(AsmOrigin &orig);
	int find_full_match(AsmSlice &ptc, int z);
	int calc_average(AsmOrigin &orig);
	int find_chips(std::vector<int> &ichips, AsmOrigin &orig, int z);
	int assemble_slice(AsmOrigin &orig, int z);
	void assemble(AsmOrigin &orig);
	bool resolve_2(AsmOrigin& first, AsmOrigin& second);
	bool can_merge_2(AsmOrigin& first, AsmOrigin& second);
	void merge_2(AsmOrigin& first, AsmOrigin& second);
	void separate_slice_2(AsmOrigin& first, AsmOrigin& second, int z);
	void separate_2(AsmOrigin& first, AsmOrigin& second);
	void decouple_2(AsmOrigin& first, AsmOrigin& second);
	bool resolve_3(AsmOrigin& first, AsmOrigin& second, AsmOrigin &third);
	void separate_3(AsmOrigin& first, AsmOrigin& second, AsmOrigin &third);
	void resolve_m(std::vector<int>& conflicted);

};


struct NucKey
{
	int idxs[5];
	int len;
	NucKey() {
		idxs[0] = idxs[1] = idxs[2] = idxs[3] = idxs[4] = -1;
		len = 0;
	}
	bool add_nuc_idx(int idx) {
		if (len >= 5) return false;
		idxs[len++] = idx;
		return true;
	}
	bool has_nuc_idx(int nuc_idx) {
		for (int i=0; i<len; i++)
			if (idxs[i] == nuc_idx) return true;
		return false;
	}
	bool empty() { return len == 0; }
	bool equals(NucKey& other) {
		if (len != other.len)
			return false;
		for (int i=0; i<len; i++)
			if (idxs[i] != other.idxs[i])
				return false;
		return true;
	}
	bool operator == (NucKey& other) { return equals(other); }
	bool operator != (NucKey& other) { return !equals(other); }
};

struct SliceId
{
	int z, idx;
	SliceId() {}
	SliceId(int _z, int _idx) : z(_z), idx(_idx) {}
};

struct SliceNucMatch
{
	NucKey nk;
	std::vector<SliceId> sliceMap;
	//
	SliceNucMatch(NucKey _nk) : nk(_nk) {}
	void add_slice(int z, int slice_idx) {
		sliceMap.push_back(SliceId(z, slice_idx));
	}
	int len() { return int(sliceMap.size()); }
	bool conflicts(SliceNucMatch& other) {
		if (other.len() == 0) return false;
		for (int j=0; j<nk.len; j++) {
			if (other.nk.has_nuc_idx(nk.idxs[j])) return true;
		}
		return false;
	}
	void invalidate() { sliceMap.clear(); }
	void merge(SliceNucMatch& other) {
		for (SliceId& si : other.sliceMap) {
			add_slice(si.z, si.idx);
		}
		other.invalidate();
	}
};

struct NucCellMatch
{
	int cell_idx;
	int nuc_idx;
	double ovl_fract;
	bool valid;
	NucCellMatch() { valid=true; }
	NucCellMatch(int _cell_idx, int _nuc_idx, double _ovl_fract) :
		cell_idx(_cell_idx), nuc_idx(_nuc_idx), ovl_fract(_ovl_fract) { valid=true; }
};

class ActinAssembler3d : public Assembler3d
{
public:
	double nuc_ovl = 0.7;
	double seed_ovl = 0.75;
	double nuc_vol_min_fract = 0.15;
	double min_slice_circ = 0.25;
	//
	std::vector<Nucleus> nuclei;
	std::vector<AsmSlice> z01slices;
	//
	std::vector<SliceNucMatch> sliceNucMap;
	std::vector<NucCellMatch> nucCellMap;
	//
	ActinAssembler3d(int _w, int _h, int _d, unsigned char *mask3d) :
			Assembler3d(_w, _h, _d, mask3d) {
		chip_ovl = 0.65;
	}
	//
	int read_z01_data(const char *z01_csvfile);
	int z01_slice_linkage();
	int read_dna_data(const char *dna_csvfile);
	int nuc_slice_linkage();
	int initial_linkage();
	int nuc_cell_linkage();
	int resolve_nuclear_conflicts();
	int resolve_conflicts();
	int fill_missing_slices();
	
	//
	virtual int estimate_zbounds(AsmOrigin &orig) override;
	
protected:
	//
	NucKey key_for_slice(int z, int slice_idx);
	int find_or_create_slicenuc(NucKey nk) {
		for (int i=0; size_t(i)<sliceNucMap.size(); i++) {
			SliceNucMatch& snm = sliceNucMap[i];
			if (snm.nk == nk) return i;
		}
		sliceNucMap.push_back(SliceNucMatch(nk));
		return int(sliceNucMap.size() - 1);
	}
	void split_slicenuc_2(int sidx);
	//
	std::vector<int> find_cells_bynuc(int nuc_idx) {
		std::vector<int> res;
		for (NucCellMatch &ncm : nucCellMap) {
			if (ncm.valid && ncm.nuc_idx == nuc_idx)
				res.push_back(ncm.cell_idx);
		}
		return res;
	}
	std::vector<int> find_nuclei_bycell(int cell_idx) {
		std::vector<int> res;
		for (NucCellMatch &ncm : nucCellMap) {
			if (ncm.valid && ncm.cell_idx == cell_idx)
				res.push_back(ncm.nuc_idx);
		}
		return res;
	}
	int count_nuclei_bycell(int cell_idx) {
		int cnt = 0;
		for (NucCellMatch &ncm : nucCellMap) {
			if (ncm.valid && ncm.cell_idx == cell_idx)
				++cnt;
		}
		return cnt;
	}
	void invalidate_nuccell(int cell_idx) {
		for (NucCellMatch &ncm : nucCellMap) {
			if (ncm.cell_idx == cell_idx)
				ncm.valid = false;
		}
	}
	//
	// int split_slice(Slice &seed, int z, int slice_idx);
	long long overlay_volume(Nucleus &nuc, AsmOrigin &orig);
	void update_nuc_cell_map(int cell_idx, int nuc_idx, double ovl_fract);
	void calculate_seed(AsmSlice &seed, AsmOrigin &orig, int z, std::vector<int> &nucids);
	int slice_from_seed(AsmSlice& seed, AsmOrigin& orig);
};

class TopBottomDetector
{
public:
	double min_nuc_fract = 0.4;
	int minzrad = 2;
	double min_slice_circ = 0.25;
	//
	int w, h, d;
	Raster3D mstack;
	Raster8 dat;
	std::vector<Nucleus> nuclei;
	std::vector<Cell> cells;
	std::vector<NucCellMatch> nucCellMap;
	//
	TopBottomDetector(int _w, int _h, int _d, unsigned char *mask3d) :
			w(_w), h(_h), d(_d),
			mstack(_w, _h, _d, mask3d), dat(_w, _h, NULL) {
	}
	//
	int read_dna_data(const char *dna_csvfile);
	int read_actin_data(const char *actin_csvfile);
	void map_cells_to_nuclei();
	void actin_cell_limits();
	void detect_top_bottom();
	//
protected:
	std::vector<Nucleus *> find_cell_nuclei(int cell_idx) {
		std::vector<Nucleus *> res;
		for (NucCellMatch &nm : nucCellMap) {
			if (nm.valid && nm.cell_idx == cell_idx)
				res.push_back(&nuclei[nm.nuc_idx]);
		}
		return res;
	}
	void combo_slice(Slice &outer, Slice &inner, Cell &cell, int z_lo, int z_hi);
};


//----- Experimental stuff (compare to REShAPE, etc.) -----

struct SliceMatch
{
	int idx;
	int rs_idx;
	int rs_ovl;
	SliceMatch() {}
	SliceMatch(int _idx, int _rs_idx, int _rs_ovl) :
		idx(_idx), rs_idx(_rs_idx), rs_ovl(_rs_ovl) {}
};

class SegmentationComparator
{
public:
	double movl = 0.8;
	unsigned char bordc = 0x10;
	unsigned char tmpc = 0x50;
	unsigned char bk = 0;
	unsigned char fg = 0xFF;
	//
	int w, h;
	Raster8 msk;
	Raster8 rmsk;
	//
	std::vector<Slice> cells;
	std::vector<Slice> rs_cells;
	std::vector<SliceMatch> matchmap;
	//
	SegmentationComparator(int _w, int _h, unsigned char *mask, unsigned char *rs_mask) :
		w(_w), h(_h), msk(_w, _h, mask), rmsk(_w, _h, rs_mask) {
	}
	//
	bool has_rs_match(int idx) {
		for (SliceMatch &sm : matchmap) {
			if (sm.idx == idx) return true;
		}
		return false;
	}
	//
	void load_cells();
	void load_rs_cells(const char *csvfile);
	void match_cells();
	int find_matches(std::vector<Slice *> & matches, int rs_idx);
};


#endif
